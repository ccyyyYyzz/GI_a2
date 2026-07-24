#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
WAVE-OPTICS DIGITAL TWIN of the ROUND63 beyond-band bucket-optics bench.

Scalar angular-spectrum (band-limited) propagation, float64 / complex128, GPU.
This is the wave-physics realization of the thin-screen statistical model in
docs/ROUND63_GPT_ROUND39_RULING_RAW.md (eqs 1.1-1.4) and the developed-speckle
endpoint in SCRAMBLE_EXT/DERIVATION.md.  It feeds the bench / collaboration
phase only -- QUARANTINED from the PRL manuscript (R43 green light).  No git.

Optical chain (per exposure):
    DMD code amplitude  --prop(z1)-->  thin random phase screen (diffuser)
    --prop(z2)-->  object plane  --|.|^2 * x(r)-->  bucket (one scalar) + shot.

Design choices (stated + verified in check_sampling()):
  * lambda = 532 nm.
  * Grid N = 2048, pixel dx = 4.0 um  => window 8.192 mm.
  * Band-limited angular-spectrum method (Matsushima-Shimobaba); the plain-ASM
    alias-free distance is z_max = N*dx^2/lambda = 61.6 mm, comfortably beyond
    every propagation leg here (z1 <= 20, z2 <= 20 mm).  We ADDITIONALLY apply
    the analytic band limit so residual aliasing is bounded for all distances.
  * DMD: 10.8 um micromirrors, fill factor 0.92.  The 64x64 logical code is
    rendered on macro-pixels of 32 um (= 8 sim px = ~3 micromirrors); the code
    band cut-off k_p=5 is enforced spectrally (band_modes) exactly as in the
    statistical model, and the macro-pixel HOLD + micromirror aperture supply
    the real pixelation envelope (the thing T1 probes for wall leakage).

Author: wave-twin build for ROUND63 bench transfer.  Honest ledger throughout.
"""
import numpy as np
import torch

# ----------------------------------------------------------------- device / dtype
# WAVE_DTYPE=64 (default, precision-critical T1) or 32 (fast complex64 for the T2-T5
# Monte Carlo).  In fast mode FIELD FFTs run in complex64 (FP32, ~4-8x faster on the
# 4060's weak FP64 path); the angular-spectrum transfer function is still BUILT in
# float64 then cast, so its phase stays accurate.  Energy conservation verified in both.
import os as _os
DEV = "cuda" if torch.cuda.is_available() else "cpu"
_FAST = _os.environ.get("WAVE_DTYPE", "64") == "32"
CDT = torch.complex64 if _FAST else torch.complex128
RDT = torch.float32 if _FAST else torch.float64

# ----------------------------------------------------------------- optical ledger
LAMBDA = 0.532e-6            # m  (532 nm)
N = 2048                    # grid
DX = 4.0e-6                 # m  pixel pitch
WIN = N * DX                # 8.192 mm window
MICRON = 1e-6

# DMD geometry
MIRROR_PITCH = 10.8e-6      # m  micromirror pitch
FILL = 0.92                 # areal fill factor
MACRO_PX = 8               # sim pixels per macro-pixel  (8*4um = 32 um ~ 3 mirrors)
NMACRO = 64                # 64x64 logical code / scene grid (matches sealed n=64)
DMD_PIX = NMACRO * MACRO_PX # 512 px active DMD region
KP = 5                     # code band cut-off (Chebyshev), matches sealed KP

# photon ledger (matches sealed_common: 1e4 detected photoelectrons / bucket)
PHOT = 1e4

_cx = N // 2
# spatial-frequency grid (cycles / m); ALWAYS float64 for transfer-function accuracy
_fx = torch.fft.fftfreq(N, d=DX).to(device=DEV, dtype=torch.float64)   # (N,)
FX, FY = torch.meshgrid(_fx, _fx, indexing="ij")              # (N,N) cyc/m float64


# ============================================================ propagation
def _bandlimited_H(z):
    """Band-limited angular-spectrum transfer function for distance z (m).

    H = exp(i k z sqrt(1 - (l fx)^2 - (l fy)^2)) for propagating waves, with the
    Matsushima-Shimobaba local-frequency band limit that zeroes components whose
    quadratic-phase gradient exceeds the sampling Nyquist -> alias-free at any z.
    """
    l = LAMBDA
    fx2 = (l * FX) ** 2
    fy2 = (l * FY) ** 2
    arg = 1.0 - fx2 - fy2
    prop = arg > 0.0                                  # propagating (not evanescent)
    sq = torch.sqrt(torch.clamp(arg, min=0.0))        # float64
    phase = (2.0 * np.pi / l) * z * sq                # float64 phase
    H = torch.complex(torch.cos(phase), torch.sin(phase))   # unit-modulus, float64
    # local fringe frequency of H: |d phi / d fx| / (2 pi) <= 1/(2 dfx) Nyquist
    dfx = float(_fx[1] - _fx[0])
    flim_x = 1.0 / (l * np.sqrt((2.0 * dfx * z) ** 2 + 1.0))
    band = (FX.abs() <= flim_x) & (FY.abs() <= flim_x) & prop
    H = H * band.to(torch.complex128)
    return H.to(CDT)                                   # cast to field dtype (c64 fast / c128)


_H_CACHE = {}


def propagate(u, z):
    """Band-limited angular-spectrum propagation of complex field u by z (m)."""
    if abs(z) < 1e-12:
        return u
    key = round(z, 12)
    H = _H_CACHE.get(key)
    if H is None:
        H = _bandlimited_H(z)
        _H_CACHE[key] = H
    U = torch.fft.fft2(u)
    return torch.fft.ifft2(U * H)


# ============================================================ DMD code rendering
def _macro_to_full(code64):
    """Upsample a 64x64 macro array to the full-grid DMD region (macro-pixel HOLD),
    embedded centered in the NxN window with an opaque surround.  Adds the finite
    micromirror APERTURE envelope by low-pass with the 10.8 um mirror footprint."""
    c = torch.as_tensor(code64, device=DEV, dtype=RDT)
    # macro-pixel hold: repeat each cell MACRO_PX times
    held = c.repeat_interleave(MACRO_PX, 0).repeat_interleave(MACRO_PX, 1)  # (512,512)
    full = torch.zeros(N, N, device=DEV, dtype=RDT)
    a = _cx - DMD_PIX // 2
    full[a:a + DMD_PIX, a:a + DMD_PIX] = held
    return full


def dmd_field(code64):
    """Complex field just AFTER the DMD for a real amplitude code (band-limited
    to k_p already).  Includes fill-factor throughput and the micromirror aperture
    envelope (a smooth low-pass ~ FT of a 10.8um*sqrt(FILL) square) applied in
    Fourier so sub-macro pixelation diffraction is represented without needing to
    resolve the ~0.9um mirror gap on the grid."""
    amp = _macro_to_full(code64)
    # micromirror aperture envelope: sinc from a square of side w = sqrt(FILL)*pitch
    w = np.sqrt(FILL) * MIRROR_PITCH
    env = torch.sinc(w * FX) * torch.sinc(w * FY)         # separable sinc (unitless)
    F = torch.fft.fft2(amp.to(CDT))
    field = torch.fft.ifft2(F * env.to(CDT))
    return field * np.sqrt(FILL)                          # areal throughput


# ============================================================ codes / scenes (n=64)
_XX, _YY = np.meshgrid(np.arange(NMACRO), np.arange(NMACRO), indexing="ij")


def band_modes(k_lo, k_hi):
    """Orthonormal REAL Fourier modes, Chebyshev freq max(|kx|,|ky|) in [k_lo,k_hi].
    Identical convention to sealed_common.band_modes (64-grid)."""
    cols = []
    for kx in range(-k_hi, k_hi + 1):
        for ky in range(-k_hi, k_hi + 1):
            if k_lo <= max(abs(kx), abs(ky)) <= k_hi:
                ph = 2 * np.pi * (kx * _XX + ky * _YY) / NMACRO
                cols.append(np.cos(ph).ravel())
                cols.append(np.sin(ph).ravel())
    Q, R = np.linalg.qr(np.array(cols).T)
    return Q[:, np.abs(np.diag(R)) > 1e-8]


def signed_codes(M=64, seed=10):
    """M signed band-limited codes in [-1,1], band 1<=|k|<=k_p (DC-free, zero-mean).
    Matches sealed_common.signed_codes up to M."""
    rp = np.random.default_rng(seed)
    Uc = band_modes(1, KP)
    P = rp.standard_normal((M, Uc.shape[1])) @ Uc.T
    P /= np.abs(P).max(axis=1, keepdims=True)
    return P.reshape(M, NMACRO, NMACRO)


def witness_scene(seed=4):
    """In-band content + deliberate beyond-band structure -> nonneg image in [0,1].
    Matches sealed_common.witness_scene."""
    rs = np.random.default_rng(seed)
    U_in = band_modes(0, KP)
    a_in = U_in @ rs.standard_normal(U_in.shape[1])
    Ube = band_modes(KP + 1, 9)
    a_be = Ube @ (2.0 * rs.choice([-1, 1], Ube.shape[1]))
    x = a_in + a_be
    x -= x.min()
    x = x / x.max()
    return x.reshape(NMACRO, NMACRO)


def beyond_band_delta(eps, xnorm, k_lo=KP + 1, k_hi=9, seed=7):
    """Sub-DMD-resolution beyond-band change delta with ||delta|| = eps*||x||.
    Realized on 64-grid Fourier modes in the beyond-band annulus [k_lo,k_hi]."""
    rng = np.random.default_rng(seed)
    Ube = band_modes(k_lo, k_hi)
    c = rng.standard_normal(Ube.shape[1])
    d = Ube @ c
    d = d / np.linalg.norm(d) * eps * xnorm
    return d.reshape(NMACRO, NMACRO)


def single_freq_scene(k, seed=0):
    """A single beyond/near-band cosine grating at Chebyshev freq k (for T3 aperture)."""
    ph = 2 * np.pi * (k * _XX) / NMACRO
    g = np.cos(ph)
    return (g / np.linalg.norm(g)).reshape(NMACRO, NMACRO)


def object_full(x64):
    """Embed a 64x64 object (reflectance/transmittance in [0,1]) on the full grid,
    aligned with the DMD footprint; opaque (0) outside."""
    return _macro_to_full(x64)


# ============================================================ diffuser (thin screen)
def make_screen(l_c_um, sigma_phi, big=None, seed=0):
    """Precompute a LARGE Gaussian-correlated random phase screen (radians).

    Gaussian surface autocorrelation of correlation length l_c (um); roughness set
    by sigma_phi (rad RMS).  developed regime: sigma_phi >= 2*pi; weak/thin-
    multiplicative regime: sigma_phi ~ 0.3 rad.  Returns a (big,big) tensor; a
    per-bank lateral SHIFT crops an NxN window = the bench 'rotation' protocol.
    """
    if big is None:
        big = N + 1024                        # room for lateral shifts
    l_c = l_c_um * MICRON
    sig_px = (l_c / DX)                        # correlation length in pixels
    rng = np.random.default_rng(seed)
    white = rng.standard_normal((big, big))
    # Gaussian filter in Fourier: corr length l_c  <=> spectral width
    fb = np.fft.fftfreq(big)                  # cycles/pixel
    FXb, FYb = np.meshgrid(fb, fb, indexing="ij")
    # Gaussian autocov exp(-r^2/(2 s^2)) with s=sig_px -> spectrum exp(-2 pi^2 s^2 f^2)
    Hk = np.exp(-2.0 * np.pi ** 2 * sig_px ** 2 * (FXb ** 2 + FYb ** 2))
    filt = np.fft.ifft2(np.fft.fft2(white) * np.sqrt(Hk)).real
    filt = (filt - filt.mean()) / (filt.std() + 1e-12)
    phase = filt * sigma_phi                  # target RMS phase
    return torch.tensor(phase, device=DEV, dtype=RDT)


def screen_crop(screen, ox, oy):
    """Crop an NxN window from the big screen at offset (ox,oy) and return exp(i phi)."""
    sub = screen[ox:ox + N, oy:oy + N]
    return torch.exp(1j * sub.to(CDT))


def bank_offsets(n_banks, big, step_px, seed=0):
    """Lateral shift offsets for n_banks banks (decorrelated between banks).
    step_px >> correlation length -> independent banks (matches sealed protocol)."""
    rng = np.random.default_rng(seed)
    maxoff = big - N
    offs = rng.integers(0, maxoff, size=(n_banks, 2))
    return offs


# ============================================================ forward: buckets
def object_plane_intensity(code_field_at_diffuser, screen_phase_crop, z2):
    """Given the code field already at the diffuser plane, apply the diffuser phase,
    propagate z2 to the object, and return |field|^2 on the full grid."""
    e = code_field_at_diffuser * screen_phase_crop
    if z2 > 0:
        e = propagate(e, z2)
    return (e.real ** 2 + e.imag ** 2)


def buckets_for_bank(code_fields_diff, screen_phase_crop, z2, x_full, phot=PHOT,
                     scp=None, shot=True, gen=None):
    """Vectorized buckets for a stack of code fields under ONE diffuser realization.

    code_fields_diff: (M,N,N) complex, code fields already propagated to the diffuser.
    Returns signed buckets b (M,) = sum_r |E|^2 x(r), with optional Poisson shot at
    `phot` mean photoelectrons per bucket (scp = photoelectrons per unit flux)."""
    e = code_fields_diff * screen_phase_crop[None]
    if z2 > 0:
        e = propagate(e, z2)
    S = e.real ** 2 + e.imag ** 2                        # (M,N,N) intensity
    flux = (S * x_full[None]).sum(dim=(-2, -1))          # (M,) bucket flux
    if shot:
        if scp is None:
            scp = phot / flux.mean().clamp(min=1e-30)
        lam = (flux * scp).clamp(min=0.0)
        counts = torch.poisson(lam, generator=gen)
        return counts / scp
    return flux


# ============================================================ sampling verification
def check_sampling():
    """Return a dict of sampling / anti-aliasing diagnostics for the chosen grid."""
    dfx = float(_fx[1] - _fx[0])
    out = dict(
        lambda_nm=LAMBDA * 1e9, N=N, dx_um=DX * 1e6, window_mm=WIN * 1e3,
        macro_px=MACRO_PX, macro_um=MACRO_PX * DX * 1e6, dmd_active_mm=DMD_PIX * DX * 1e3,
        mirror_pitch_um=MIRROR_PITCH * 1e6, fill=FILL, kp=KP,
        asm_alias_free_mm=(N * DX ** 2 / LAMBDA) * 1e3,       # plain-ASM critical distance
        nyquist_freq_cyc_per_mm=(0.5 / DX) / 1e3,
        max_scene_freq_cyc_per_mm=(NMACRO / 2) / (NMACRO * MACRO_PX * DX) / 1e3,
    )
    return out


if __name__ == "__main__":
    import json
    s = check_sampling()
    print(json.dumps(s, indent=2))
    print("DEV", DEV)
