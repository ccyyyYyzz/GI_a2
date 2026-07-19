"""Shared helpers for the Optics-Express SUPPLEMENT table generators
(code/round63/figs/supp_s*.py).

Every supp_s*.py script reads a frozen campaign artifact, prints the numbers it
extracts to stdout (auditable), and writes a self-contained LaTeX table fragment
to paper/supp_tables/ that paper/supplement.tex \\input's.

PROVENANCE: values are read verbatim from the committed artifacts; no imputation
is performed. Empty / 'nan' cells are dropped from every statistic (never
imputed). Each fragment carries a provenance comment naming its source file and
the git commit that last touched it.

Run (cwd = repo root D:\\GI_another), e.g.:
    D:\\Anacondar\\anaconda3\\python.exe code\\round63\\figs\\supp_s4_exact.py
"""
import os
import csv
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))          # code/round63/figs
REPO = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))  # repo root
OUT_DIR = os.path.join(REPO, "paper", "supp_tables")


def to_float(s):
    """Parse a cell to float; '' / 'nan' / 'inf' / '--' -> None (dropped, never imputed)."""
    if s is None:
        return None
    s = str(s).strip()
    if s in ("", "--", "None") or s.lower() in ("nan", "inf", "-inf", "+inf"):
        return None
    try:
        return float(s)
    except ValueError:
        return None


def load_csv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def git_commit(path):
    """Short hash of the commit that last touched `path` (abs or repo-relative)."""
    try:
        h = subprocess.check_output(
            ["git", "log", "-1", "--format=%h", "--", path],
            cwd=REPO, text=True, stderr=subprocess.DEVNULL).strip()
        return h if h else "UNTRACKED"
    except Exception:
        return "UNKNOWN"


def tex_escape(s):
    """Escape LaTeX specials in free text / identifiers."""
    s = str(s)
    out = []
    for ch in s:
        if ch == "\\":
            out.append(r"\textbackslash{}")
        elif ch in "_%&#${}":
            out.append("\\" + ch)
        elif ch == "~":
            out.append(r"\textasciitilde{}")
        elif ch == "^":
            out.append(r"\textasciicircum{}")
        else:
            out.append(ch)
    return "".join(out)


def tt(s):
    """Monospace identifier with escaped specials."""
    return r"\texttt{%s}" % tex_escape(s)


def fnum(x, nd=2, plus=False):
    """Format a float to nd decimals; None -> en-dash. plus=True keeps + sign."""
    if x is None:
        return "--"
    fmt = ("%+." + str(nd) + "f") if plus else ("%." + str(nd) + "f")
    return fmt % x


def md_tables(text):
    """Parse every pipe-delimited GitHub-markdown table in `text`.

    Returns a list of (header, rows) where header is a list of column names and
    rows is a list of lists of cell strings. A table is a maximal run of lines
    beginning with '|', whose second line is a '---' separator.
    """
    tables = []
    lines = text.splitlines()
    i = 0
    n = len(lines)
    while i < n:
        ln = lines[i].strip()
        if ln.startswith("|") and i + 1 < n and set(lines[i + 1].strip()) <= set("|-: "):
            # header at i, separator at i+1, data from i+2
            header = [c.strip() for c in ln.strip("|").split("|")]
            rows = []
            j = i + 2
            while j < n and lines[j].strip().startswith("|"):
                cells = [c.strip() for c in lines[j].strip().strip("|").split("|")]
                rows.append(cells)
                j += 1
            tables.append((header, rows))
            i = j
        else:
            i += 1
    return tables


def find_table(text, must_have):
    """Return the first (header, rows) whose header contains all `must_have` names."""
    for header, rows in md_tables(text):
        hs = set(header)
        if all(m in hs for m in must_have):
            return header, rows
    raise SystemExit("no markdown table with columns %s" % (must_have,))


def write_fragment(name, lines):
    os.makedirs(OUT_DIR, exist_ok=True)
    path = os.path.join(OUT_DIR, name)
    with open(path, "w", newline="\n") as f:
        f.write("\n".join(lines) + "\n")
    print("[supp] wrote %s" % path, flush=True)
    return path
