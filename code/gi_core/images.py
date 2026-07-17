"""Ground-truth image loading — spec §1. Canonical: data/*.png (64x64, 8-bit).
Truth is sum-normalized (sum x = 1) before use; Phase A uses bicubic 16x16
downsamples of the same images.
"""
import os

import numpy as np

from . import config as C

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.path.join(ROOT, "data")


def load_truths(side):
    """Returns dict name -> flattened float64 truth with sum = 1."""
    from imageio.v2 import imread
    from skimage.transform import resize

    out = {}
    for name in C.IMAGES:
        img = imread(os.path.join(DATA, "%s.png" % name)).astype(np.float64) / 255.0
        assert img.shape == (64, 64)
        if side != 64:
            img = np.clip(resize(img, (side, side), order=3, anti_aliasing=True,
                                 mode="reflect"), 0.0, 1.0)
        x = img.ravel()
        s = x.sum()
        assert s > 0
        out[name] = x / s
    return out
