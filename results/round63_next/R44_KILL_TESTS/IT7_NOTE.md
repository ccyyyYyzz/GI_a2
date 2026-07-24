# IT7 — Contact-arm coincidence veto (internal divergence stranger #6; D3/D5 FA repair)

**Verdict: KILL by power cost — the hardware identity is real and the medium false
alarm IS repaired (17.8×), but a beyond-band scene change is too weak in the contact
channel to survive the AND coincidence (only 24% scene power retained). Fall back to the
quarantined semiparametric/cross-fit statistical repair.** `IT7_contact_coincidence_veto.py`
+ `.json`. M=16, 700 banks, T_eff=300, N_mc=2000, ~6 min on one L4.

## Design
- **Contact arm (z2_eff=0):** the diffuser phase cancels in |E·e^{iφ}|², so the contact
  bucket is **exactly** medium-independent (hardware identity) and sees a scene change only
  through the mean channel (Δb_c = J_mean·δ). Statistic = matched mean filter, T shots.
- **z2=5 mm arm:** a **generic covariance-deviation** detector ‖V⁻¹ᐟ²(Ĉ−C₀)V⁻¹ᐟ²‖²_F —
  fires on ANY covariance change (scene OR medium) = the FA-prone primary.
- **2D coincidence** thresholded on H0 (each arm at its 95th percentile).

## Measured
| quantity | value |
|---|---|
| contact AUC(H0 vs **medium**) | **0.5009** (exactly medium-blind — hardware identity ✓) |
| contact AUC(H0 vs scene) | 0.741 (sees the beyond-band change only moderately) |
| z2=5 AUC(H0 vs medium) | 0.999 (medium strongly fires the primary = the FA source) |
| z2=5 AUC(H0 vs scene) | 0.993 |
| single-arm z2=5 power (scene) | 0.971 |
| **single-arm z2=5 FA (medium)** | **0.9965** (the sealed campaign's failed gate) |
| coincidence FA (medium) | **0.056** (17.8× repair) |
| coincidence power (scene) | 0.233 |
| **power ratio coinc/single** | **0.24** (76% of scene power lost) |

## Reading (honest KILL with strong diagnostics)
Two of the three decision criteria PASS: the contact arm is exactly medium-blind
(AUC 0.5009 ≈ 0.5, the flip-mirror hardware identity holds), and the coincidence cuts
the medium false alarm from **0.9965 (single-arm) to 0.056 (17.8× repair)**, essentially
at the 0.05 target. **The third fails**: the coincidence retains only **24%** of single-arm
scene power (≪ the 90% floor), because a beyond-band object change leaks only weakly into
the contact/mean channel (contact scene AUC 0.74), so demanding both arms fire discards
most true detections. This is the coordinator's explicit KILL branch ("power cost >20% …
scene leak too weak to witness → fall back to the quarantined semiparametric/cross-fit
repair"). The mechanism is fundamental to *beyond-band* changes (they barely leak into the
mean channel by design — that is the wall); the veto would work for in-band/amplitude
changes but not for the campaign's beyond-band witness. Nothing here touches the Letter or
sealed artifacts.
