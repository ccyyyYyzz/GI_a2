# M1 verdicts (unblinded 2026-07-19, tag m1-freeze @ 6f00932)

RIDGE_OPERATING_PASS = True
  median dQ = 1.867 dB (bar >= 1.0)
  LB2.5 = 0.120 (bar > 0)
  n_pos = 19/24 (bar >= 18)

RIDGE_SPEED_PASS = False
  median S = 0.276 (bar >= 3)
  S LB2.5 = 0.172 (bar > 1)
  n_S>1 = 1/24 (bar >= 18)

cert branch = FALLBACK_DESCRIPTIVE
  descriptive certified fraction = 0.0
  status distribution: {'NUMERICAL_UNRESOLVED': 181, 'COUNTEREXAMPLE': 299}
  by anchor: {'nu200_b0.05': {'NUMERICAL_UNRESOLVED': 43, 'COUNTEREXAMPLE': 77}, 'nu200_b0.6': {'COUNTEREXAMPLE': 67, 'NUMERICAL_UNRESOLVED': 53}, 'nu2000_b0.05': {'COUNTEREXAMPLE': 84, 'NUMERICAL_UNRESOLVED': 36}, 'nu2000_b0.6': {'COUNTEREXAMPLE': 71, 'NUMERICAL_UNRESOLVED': 49}}

per-family median dQ: {'chirp': 1.4046299999999983, 'contour': -4.17524, 'glyph': 3.5824799999999986, 'maze': 6.368990000000002, 'microtexture': 1.4872999999999994, 'spokes': 0.71488}
