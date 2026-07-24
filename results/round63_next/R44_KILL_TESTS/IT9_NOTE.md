# IT9 — Carrier-cross exactness, grid-independence confirm (internal finer #2 residual)

**Verdict: CONFIRMED — the WAVE_TWIN_REPORT §1 carrier-cross correction can be frozen
with confidence.** `IT9_carrier_cross_grid_independence.py` + `.json`. float64, M=32,
~14 s on one L4.

## Measured
| grid | carrier-cross rel-err (max) | beyond-2kp leak (mean ± std) |
|---|---|---|
| N=2048 | **5.05e-16** | 9.581e-3 ± 1.32e-3 |
| N=4096 | **4.71e-16** | 9.581e-3 ± 1.32e-3 |

- **Carrier-cross exact (both grids):** `|E⁺|²−|E⁻|² = Re[C·conj(S)]` holds to machine
  zero (~5e-16) — the `|S|²` intensity-autocorrelation term **cancels exactly** in the
  complementary difference. There is **no 2k_p intensity-autocorrelation edge**; the
  beyond-band leak is a carrier–code heterodyne cross term, as both divergence engines
  derived.
- **Grid-independent:** the beyond-2k_p carrier-pedestal leak is **9.581e-3 at both
  N=2048 and N=4096 (0.0% difference)** — it is a real physical pedestal, not a
  2048-grid aliasing artifact.

## Reading
The report correction (no 2k_p edge; leak = carrier-cross pedestal) is exact and
grid-independent to M=32 and a doubled 4096² window. The WAVE_TWIN_REPORT §1 mechanism
sentence and the annulus-profile numbers can be frozen with confidence. Nothing here
touches the Letter or sealed artifacts.
