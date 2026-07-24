#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""Shared bank-pool engine for T3/T4/T5.

The expensive step (scalar propagation to the object-plane speckle intensity) is
OBJECT-INDEPENDENT, so we precompute the block-summed speckle kernel per code per
bank ONCE.  A bucket for ANY 64x64 object is then a cheap tensor contraction:
    flux_i(bank) = <block_sum(S_i(bank)), x_64> .
Shot noise (1e4 pe / bucket) is applied at bucket formation on the flux.
Complementary exposures (+/-) are kept separate so the SIGNED bucket b = B^+ - B^-
matches the sealed detector convention.
"""
import numpy as np
import torch
import wave_twin as wt

DEV = wt.DEV; c = wt._cx; NM = wt.NMACRO; HALF = wt.DMD_PIX // 2
CHUNK = 8

def _block_sum(Sfull):
    a = c - HALF
    reg = Sfull[..., a:a+wt.DMD_PIX, a:a+wt.DMD_PIX]
    sh = reg.shape[:-2]
    reg = reg.reshape(*sh, NM, wt.MACRO_PX, NM, wt.MACRO_PX)
    return reg.sum(dim=(-3, -1))

def code_fields_at_diffuser(M, z1, seed=10):
    """Complementary code fields propagated to the diffuser plane: Ep,Em (M,N,N)."""
    S = wt.signed_codes(M=M, seed=seed)
    Ep = torch.stack([wt.propagate(wt.dmd_field((1.0+S[i])*0.5), z1) for i in range(M)])
    Em = torch.stack([wt.propagate(wt.dmd_field((1.0-S[i])*0.5), z1) for i in range(M)])
    return Ep, Em, S

def speckle_pool(Ep, Em, z2, sigma_phi, l_c_um, n_bank, seed=100, step_px=300):
    """Precompute block-summed object-grid speckle kernels Gp,Gm (n_bank,M,64,64)
    for the + and - complementary exposures, one independent diffuser per bank."""
    M = Ep.shape[0]
    scr = wt.make_screen(l_c_um, sigma_phi, seed=seed); big = scr.shape[0]
    offs = wt.bank_offsets(n_bank, big, step_px, seed=seed+1)
    Gp = torch.empty(n_bank, M, NM, NM, device=DEV, dtype=wt.RDT)
    Gm = torch.empty(n_bank, M, NM, NM, device=DEV, dtype=wt.RDT)
    for t in range(n_bank):
        ph = wt.screen_crop(scr, int(offs[t,0]), int(offs[t,1]))
        for s in range(0, M, CHUNK):
            ep = Ep[s:s+CHUNK]*ph[None]; em = Em[s:s+CHUNK]*ph[None]
            if z2 > 0:
                ep = wt.propagate(ep, z2); em = wt.propagate(em, z2)
            Gp[t, s:s+CHUNK] = _block_sum(ep.real**2 + ep.imag**2)
            Gm[t, s:s+CHUNK] = _block_sum(em.real**2 + em.imag**2)
    return Gp, Gm

def signed_buckets(Gp, Gm, x64, phot=wt.PHOT, shot=True, gen=None, scp=None):
    """Signed buckets b (n_bank,M) = B^+ - B^- for object x64, with Poisson shot at
    `phot` mean pe/bucket.  scp (pe per unit flux) is shared across H0/H1 if given."""
    x = torch.as_tensor(x64, device=DEV, dtype=wt.RDT)
    fp = torch.einsum('bmij,ij->bm', Gp, x)
    fm = torch.einsum('bmij,ij->bm', Gm, x)
    if not shot:
        return fp - fm, None
    if scp is None:
        scp = phot / (0.5*(fp+fm)).mean().clamp(min=1e-30)
    npos = torch.poisson((fp*scp).clamp(min=0), generator=gen)
    nneg = torch.poisson((fm*scp).clamp(min=0), generator=gen)
    return (npos - nneg)/scp, scp
