"""Generate fig/blosum62.svg (coloured BLOSUM62 matrix used in docs/sequence_alignment.md).

Matrix values: https://ftp.ncbi.nih.gov/blast/matrices/BLOSUM62 (20 standard amino acids).
Amino acids are ordered by side-chain class, same as in fig/aa_codes.svg.

Usage: python3 scripts/make_blosum62.py
PNG version: rsvg-convert -z 2 fig/blosum62.svg -o fig/blosum62.png
"""
from pathlib import Path

BLOSUM62 = """
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4
"""

lines = BLOSUM62.strip().splitlines()
cols = lines[0].split()
M = {}
for line in lines[1:]:
    f = line.split()
    for c, v in zip(cols, f[1:]):
        M[f[0], c] = int(v)
assert all(M[a, b] == M[b, a] for a in cols for b in cols)

# order and colours by side-chain class (matches fig/aa_codes.svg)
groups = [
    ("negative", "#1f5f99", "DE"),
    ("positive", "#a3282c", "RKH"),
    ("polar", "#8a6400", "NQSTY"),
    ("nonpolar", "#3c6e1f", "AGVLIPFMWC"),
]
order = "".join(g[2] for g in groups)
label_color = {aa: col for _, col, aas in groups for aa in aas}


def cell_color(v):
    """Diverging scale: red for negative, white for 0, blue for positive."""
    if v == 0:
        return "#f7f7f7"
    if v < 0:
        t = min(-v, 4) / 4
        a, b = (247, 247, 247), (214, 96, 77)
    else:
        t = min(v, 11) / 11
        t = t ** 0.6
        a, b = (247, 247, 247), (33, 102, 172)
    r, g, bl = (round(a[i] + (b[i] - a[i]) * t) for i in range(3))
    return f"#{r:02x}{g:02x}{bl:02x}"


C = 40            # cell size
X0, Y0 = 70, 90   # top-left of the matrix
N = len(order)
W = X0 + N * C + 30
H = Y0 + N * C + 90
FONT = "font-family='DejaVu Sans, Arial, Helvetica, sans-serif'"
out = []
add = out.append

add(f"<text x='{X0}' y='30' {FONT} font-size='20' font-weight='bold' fill='#222'>BLOSUM62 substitution matrix</text>")

for i, aa in enumerate(order):
    col = label_color[aa]
    add(f"<text x='{X0 + i*C + C/2}' y='{Y0-12}' {FONT} font-size='18' font-weight='bold' fill='{col}' text-anchor='middle'>{aa}</text>")
    add(f"<text x='{X0-16}' y='{Y0 + i*C + C/2 + 6}' {FONT} font-size='18' font-weight='bold' fill='{col}' text-anchor='middle'>{aa}</text>")

for i, a in enumerate(order):
    for j, b in enumerate(order):
        v = M[a, b]
        x, y = X0 + j * C, Y0 + i * C
        add(f"<rect x='{x}' y='{y}' width='{C}' height='{C}' fill='{cell_color(v)}' stroke='#ffffff' stroke-width='1'/>")
        dark = (v >= 5) or (v <= -4)
        fill = "#ffffff" if dark else "#222"
        weight = "bold" if i == j else "normal"
        add(f"<text x='{x + C/2}' y='{y + C/2 + 6}' {FONT} font-size='16' font-weight='{weight}' fill='{fill}' text-anchor='middle'>{v}</text>")

# group separators
pos = 0
for _, col, aas in groups[:-1]:
    pos += len(aas)
    p = pos * C
    add(f"<line x1='{X0 + p}' x2='{X0 + p}' y1='{Y0}' y2='{Y0 + N*C}' stroke='#555' stroke-width='2'/>")
    add(f"<line x1='{X0}' x2='{X0 + N*C}' y1='{Y0 + p}' y2='{Y0 + p}' stroke='#555' stroke-width='2'/>")
add(f"<rect x='{X0}' y='{Y0}' width='{N*C}' height='{N*C}' fill='none' stroke='#555' stroke-width='2'/>")

# legend
ly = Y0 + N * C + 30
add(f"<text x='{X0}' y='{ly + 17}' {FONT} font-size='15' fill='#444'>score:</text>")
lx = X0 + 60
for v in range(-4, 12):
    add(f"<rect x='{lx}' y='{ly}' width='30' height='24' fill='{cell_color(v)}' stroke='#ccc' stroke-width='0.5'/>")
    add(f"<text x='{lx + 15}' y='{ly + 40}' {FONT} font-size='12' fill='#444' text-anchor='middle'>{v}</text>")
    lx += 30
add(f"<text x='{lx + 14}' y='{ly + 12}' {FONT} font-size='13' fill='#444'>positive: substitution seen</text>")
add(f"<text x='{lx + 14}' y='{ly + 28}' {FONT} font-size='13' fill='#444'>more often than by chance</text>")

svg = (f"<svg xmlns='http://www.w3.org/2000/svg' width='{W}' height='{H}' viewBox='0 0 {W} {H}'>\n"
       f"<rect width='{W}' height='{H}' fill='#ffffff'/>\n" + "\n".join(out) + "\n</svg>\n")
(Path(__file__).resolve().parent.parent / "fig" / "blosum62.svg").write_text(svg)
