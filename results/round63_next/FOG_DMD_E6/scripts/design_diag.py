# Design diagnostics: Fourier bank rank/overlap + A2 schedule redundancy structure.
import numpy as np, fog_e6 as fe, json
n=32; N=1024; M=96; d=48; T=128; S=96
x=fe.make_scene(n); U=fe.dct_basis(d,n)
rngp=np.random.default_rng(777); P_bin=fe.make_binary_bank(M,N,rngp)
P_four,meta=fe.make_fourier_overlap_bank(M,n,spacing=4,band=4.0)
out={}
for name,P in [('binary',P_bin),('fourier',P_four)]:
    Pi,rangeP=fe.projectors(P)
    s=np.linalg.svd(P,compute_uv=False)
    out[name]=dict(rank=int((s>1e-9*s[0]).sum()), cond=float(s[0]/s[max((s>1e-9*s[0]).sum()-1,0)]),
                   null_frac_scene=fe.null_fraction(x,rangeP),
                   rowmean=float(P.mean()), rowstd=float(P.std()))
out['fourier_meta']=meta
# A2 schedule redundancy: distinct patterns per block, coverage
rng=np.random.default_rng(3)
idx,flat=fe.slot_idx(T,S,M,rng)
distinct=[len(np.unique(idx[t])) for t in range(T)]
counts=np.bincount(flat,minlength=M)
# medium-observation subspace variation: fraction of pattern PAIRS co-measured within a block
comeasure=np.zeros((M,M))
for t in range(T):
    u=np.unique(idx[t])
    comeasure[np.ix_(u,u)]+=1
offdiag=comeasure[~np.eye(M,dtype=bool)]
out['A2_schedule']=dict(distinct_per_block_mean=float(np.mean(distinct)),
                        distinct_per_block_min=int(np.min(distinct)),
                        distinct_per_block_max=int(np.max(distinct)),
                        pattern_count_min=int(counts.min()), pattern_count_max=int(counts.max()),
                        pattern_count_mean=float(counts.mean()),
                        pair_comeasure_mean=float(offdiag.mean()),
                        pair_comeasure_frac_zero=float((offdiag==0).mean()))
# certificate contrast A0 vs A2 (same scene)
idx0=fe.cartesian_idx(T,M)
out['cert_A0']=fe.carrier_certificate(P_bin,idx0.reshape(-1),x)
out['cert_A2']=fe.carrier_certificate(P_bin,flat,x)
print(json.dumps(out,indent=1))
