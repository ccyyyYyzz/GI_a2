# IT5 — DPSS/Slepian concentration window (internal divergence finer #4)

**Verdict: PASS — a DPSS render+object window holds the beyond-2kp witness leak
≤2.4e-4 uniformly in z1. The hardware concentration route lives** (unlike the naive
Tukey window, which IT4b killed). `IT5_dpss_concentration.py` + `.json`. Deterministic
float64, M=12, ~10 s on one L4.

## Witness leak (shells k∈[12,14], beyond 2kp=10) vs z1
| window | z1=0 | z1=1 | z1=5 | z1=10 | z1=20 |
|---|---|---|---|---|---|
| none | 5.3e-3 | 9.9e-3 | 3.9e-2 | 6.1e-2 | 8.3e-2 |
| Tukey (render+obj) | 4.0e-2 | 3.9e-2 | 3.9e-2 | 3.8e-2 | 3.3e-2 |
| **DPSS (render+obj)** | **2.4e-4** | 2.3e-4 | 2.3e-4 | 2.3e-4 | 2.1e-4 |

- **P_dpss_concentration_lives: PASS** — DPSS holds the leak **≤2.4e-4 ≪ 1e-3 uniformly**
  across the full z1 sweep (349× below the no-window pedestal, which grows to 8.3e-2 at
  z1=20mm).
- **Tukey confirms its kill** — the naive separable Tukey window sits at ~4e-2 (worse
  than no-window at z1=0), the rim effect IT4b flagged. Only the concentration-optimal
  (DPSS/Slepian) window works.

## Combined route (DPSS + single-z1 SVD subspace)
The two window/code routes do **NOT** stack: applying the DPSS object window to the
single-z1 (z1=0) SVD leak operator **shrinks** the certified-blind subspace
(d@1e-4: 80 → 24; d@1e-5: 80 → 9). Like the Tukey window in IT4b, the concentration
window trades spectral-leak suppression for code-subspace dimension. The two are
alternative wall-restoration mechanisms, used separately: DPSS for the beyond-2kp
witness leak, SVD-blind codes for the mean-blind subspace.

## Reading
The DPSS/Slepian concentration window is the surviving window route after the Tukey
rim-effect kill: it drives the mean-channel leak into the beyond-2kp witness annulus to
≤2.4e-4 at every propagation distance — well inside the 1e-3 bench spec. It provides an
exact-enough detector-level guard for the pupil-hardened wall (KT1 P2b's open item).
Nothing here touches the Letter or sealed artifacts.
