"""Shared publication style + source-table loaders for the R40 flagship figures.

Palette: Okabe-Ito colorblind-safe. Fonts >= 7 pt at single-column (8.6 cm) width.
Every figure reads ONLY from paper_aperture/source_tables/ (RULE ZERO).
"""
import csv, json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(__file__)
ST   = os.path.join(HERE, "..", "source_tables")
OUT  = HERE

# Okabe-Ito colorblind-safe palette
C = {
    "green":     "#009E73",   # certified / living / M1
    "red":       "#D55E00",   # killed / harmful (vermillion)
    "amber":     "#E69F00",   # pocket / borderline
    "blue":      "#0072B2",   # ridge / support / high-flux
    "sky":       "#56B4E9",   # secondary curves
    "gray":      "#9A9A9A",   # neutral / support fill edge
    "grayfill":  "#DCDCDC",   # pale formal-support fill
    "darkfill":  "#3A6E8F",   # Fisher-usable fill (muted blue)
    "ink":       "#222222",   # text / axes
    "line":      "#555555",
}

# cm -> inch
def cm(x): return x / 2.54

def apply_style():
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
        "font.size": 8,
        "axes.titlesize": 8.5,
        "axes.labelsize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "axes.linewidth": 0.7,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "axes.edgecolor": C["ink"],
        "text.color": C["ink"],
        "axes.labelcolor": C["ink"],
        "xtick.color": C["ink"],
        "ytick.color": C["ink"],
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "pdf.fonttype": 42,   # embed TrueType (editable text)
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })

def load_json(name):
    with open(os.path.join(ST, name), encoding="utf-8") as f:
        return json.load(f)

def load_csv(name):
    """Return (comment_lines, list-of-dict rows). Lines starting with '#' are comments."""
    comments, data_lines = [], []
    with open(os.path.join(ST, name), encoding="utf-8") as f:
        for ln in f:
            # tolerate csv-quoted comment lines (leading '"' before '#')
            if ln.lstrip().lstrip('"').startswith("#"):
                comments.append(ln.rstrip("\n"))
            else:
                data_lines.append(ln)
    rdr = csv.DictReader(data_lines)
    return comments, list(rdr)

def save(fig, stem):
    """Save both vector PDF and PNG preview."""
    pdf = os.path.join(OUT, stem + ".pdf")
    png = os.path.join(OUT, stem + ".png")
    fig.savefig(pdf, bbox_inches="tight", pad_inches=0.02)
    fig.savefig(png, bbox_inches="tight", pad_inches=0.02, dpi=200)
    print("  saved", os.path.basename(pdf), "+", os.path.basename(png))
