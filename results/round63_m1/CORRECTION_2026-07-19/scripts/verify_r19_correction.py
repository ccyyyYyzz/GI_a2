"""Small deterministic checks for the additive R19 correction bundle."""
import json
from pathlib import Path


def pava(values):
    blocks = []
    for value in values:
        blocks.append([float(value), 1])
        while len(blocks) >= 2 and blocks[-2][0] > blocks[-1][0]:
            left, right = blocks[-2:]
            weight = left[1] + right[1]
            blocks[-2:] = [[(left[0] * left[1] + right[0] * right[1]) / weight, weight]]
    return [mean for mean, _ in blocks for _ in range(1)]  # block means for test vector


def main():
    root = Path(__file__).resolve().parents[1]
    doc = json.loads((root / "M1_VERDICTS_SPEC_CORRECTED_R19.json").read_text())
    assert pava([3, 2, 1]) == [2.0]
    assert doc["verdicts"]["RIDGE_OPERATING_PASS"] is True
    assert doc["verdicts"]["RIDGE_SPEED_PASS"] is True
    assert doc["verdicts"]["posthoc_nu_rho_sensitivity_verdict"] is None
    cert = doc["full_stack_certificate_descriptive"]
    assert cert["branch"] == "FALLBACK_DESCRIPTIVE"
    assert cert["n_cert_cells"] == 480
    assert cert["n_certified"] == 0
    assert cert["n_counterexample"] == 299
    assert cert["n_numerical_unresolved"] == 181
    print("R19 correction checks: PASS")


if __name__ == "__main__":
    main()
