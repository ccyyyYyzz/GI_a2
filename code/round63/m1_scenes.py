"""M1 fresh confirmatory scene set (spec docs/ROUND63_METHOD_SPEC_M1.md §3).

24 confirmatory instances: the SAME six frozen detail32 family generators,
NEW seed namespace s_{f,r} = 633000 + 100 f + r (R10 Q7 §3); six dev-only
instances 632900 + f. Purely additive: reuses detail24's _build32 pipeline
(PNG truth of record, generator-frozen ROIs, SHA256 manifest) with a new
cache root data/r63_images_m1/. No frozen table or cache is touched.
No acceptance filter of any kind is applied (R10: no OED/C_u/Gamma/PSNR/
CNR/visual filtering).
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import detail24 as d24  # noqa: E402

M1_CONF_BASE = 633000
M1_DEV_BASE = 632900
M1_IMG_ROOT = os.path.join(d24.DATA, "r63_images_m1")
M1_DEV_IMG_ROOT = os.path.join(d24.DATA, "r63_images_m1_dev")

_FAMS = [fam for _, fam, _ in d24._conf32_table()[::4]]  # ordered six families


def _m1_conf_table():
    out = []
    for f, fam in enumerate(_FAMS):
        for r in range(4):
            out.append(("m1_%s_%d" % (fam, r), fam,
                        M1_CONF_BASE + 100 * f + r))
    return out


def _m1_dev_table():
    return [("m1_dev_%s" % fam, fam, M1_DEV_BASE + f)
            for f, fam in enumerate(_FAMS)]


def build_m1_set(side=32):
    return d24._build32(_m1_conf_table(), M1_IMG_ROOT, side,
                        ruling="M1 spec §3 (R10 Q7 fresh scenes)")


def build_m1_dev_set(side=32):
    return d24._build32(_m1_dev_table(), M1_DEV_IMG_ROOT, side,
                        ruling="M1 spec §3 dev instances")


if __name__ == "__main__":
    conf = build_m1_set()
    dev = build_m1_dev_set()
    print("[m1_scenes] conf:", len(conf), "images ->", M1_IMG_ROOT)
    print("[m1_scenes] dev :", len(dev), "images ->", M1_DEV_IMG_ROOT)
    print("[m1_scenes] families:", _FAMS)
