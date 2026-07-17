"""Fetch the 8 Stage-0 test images per ROUND59 spec §1.

Preferred source: STL-10 test split, classes {cat, deer, dog, horse} and
{airplane, car, ship, truck}, first image of each class in test-split order.
Grayscale (skimage rgb2gray, Rec.709 luma), bicubic resize 96->64
(skimage.transform.resize order=3, anti_aliasing), saved as 8-bit PNG
(canonical ground truth) + sha256.txt.

Fallback (only if download infeasible): skimage.data built-ins, declared
per-image in data/PROVENANCE.md.
"""
import hashlib
import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
RAW = os.path.join(DATA, "stl10_raw")

# spec listing order
WANTED = ["cat", "deer", "dog", "horse", "airplane", "car", "ship", "truck"]


def to_gray64(img_hwc_uint8):
    from skimage.color import rgb2gray
    from skimage.transform import resize

    g = rgb2gray(img_hwc_uint8)  # float64 in [0,1]
    g64 = resize(g, (64, 64), order=3, anti_aliasing=True, mode="reflect")
    g64 = np.clip(g64, 0.0, 1.0)
    return np.round(g64 * 255.0).astype(np.uint8)


def main():
    os.makedirs(DATA, exist_ok=True)
    if os.path.exists(os.path.join(DATA, "sha256.txt")):
        print("data already fetched (sha256.txt exists) - not overwriting; "
              "delete data/sha256.txt to refetch", flush=True)
        return 0
    provenance = {"source": None, "images": []}
    try:
        import torchvision

        ds = torchvision.datasets.STL10(root=RAW, split="test", download=True)
        classes = ds.classes  # ['airplane','bird','car','cat',...]
        idx_of = {c: classes.index(c) for c in WANTED}
        first = {}
        labels = np.asarray(ds.labels)
        for name in WANTED:
            pos = int(np.nonzero(labels == idx_of[name])[0][0])
            first[name] = pos
        provenance["source"] = "STL-10 test split (torchvision %s)" % torchvision.__version__
        imgs = {}
        for name in WANTED:
            pil_img, _ = ds[first[name]]
            arr = np.asarray(pil_img)  # HWC uint8 96x96x3
            imgs[name] = to_gray64(arr)
            provenance["images"].append(
                {"name": name, "stl10_test_index": first[name], "class_id": idx_of[name]}
            )
    except Exception as e:  # download infeasible -> declared fallback
        print("STL-10 unavailable (%r); falling back to skimage.data" % (e,), flush=True)
        from skimage import data as skdata
        from skimage.color import rgb2gray
        from skimage.transform import resize

        pool = [
            ("cat", skdata.chelsea()),
            ("deer", skdata.coffee()),
            ("dog", skdata.astronaut()),
            ("horse", np.dstack([255 - skdata.horse().astype(np.uint8) * 255] * 3)),
            ("airplane", skdata.camera()),
            ("car", skdata.brick()),
            ("ship", skdata.coins()),
            ("truck", skdata.moon()),
        ]
        provenance["source"] = "skimage.data fallback (STL-10 download infeasible)"
        imgs = {}
        for name, arr in pool:
            if arr.ndim == 2:
                g = arr.astype(np.float64) / 255.0
            else:
                g = rgb2gray(arr)
            g64 = resize(g, (64, 64), order=3, anti_aliasing=True, mode="reflect")
            imgs[name] = np.round(np.clip(g64, 0, 1) * 255.0).astype(np.uint8)
            provenance["images"].append({"name": name, "substitute": True})

    from imageio.v2 import imwrite

    sha_lines = []
    for name in WANTED:
        p = os.path.join(DATA, "%s.png" % name)
        imwrite(p, imgs[name])
        with open(p, "rb") as f:
            h = hashlib.sha256(f.read()).hexdigest()
        sha_lines.append("%s  %s.png" % (h, name))
        print("saved", p, h, flush=True)
    with open(os.path.join(DATA, "sha256.txt"), "w") as f:
        f.write("\n".join(sha_lines) + "\n")
    with open(os.path.join(DATA, "PROVENANCE.json"), "w") as f:
        json.dump(provenance, f, indent=2)
    print("DONE data_fetch", flush=True)


if __name__ == "__main__":
    sys.exit(main())
