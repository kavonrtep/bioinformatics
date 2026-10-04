"""Generate fig/blosum62.svg (coloured BLOSUM62 matrix used in docs/sequence_alignment.md).

Matrix values: https://ftp.ncbi.nih.gov/blast/matrices/BLOSUM62 (20 standard amino acids).
Lower triangle, amino acids ordered and coloured by physicochemical group
(same scheme as the BLOSUM62/PAM250 figures in the lecture slides).

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

# groups, order and colours (same as lecture slides)
GROUPS = [("C", "#C9A400", "sulfur"), ("STPAG", "#548235", "small"),
          ("NDEQ", "#C55A11", "acidic / amide"), ("HRK", "#1F4E79", "basic"),
          ("MILV", "#7EA6E0", "hydrophobic"), ("FYW", "#8E6BB0", "aromatic")]
ORDER = "".join(g for g, _, _ in GROUPS)
COL = {a: c for g, c, _ in GROUPS for a in g}
TEXT, GREY, TITLE = "#1A1A1A", "#595959", "#1F4E79"

C = 40            # cell size
X0, Y0 = 50, 90   # top-left of the matrix
N = len(ORDER)
LX = X0 + N * C + 40   # legend x
W = LX + 330
H = Y0 + N * C + 20
FONT = "font-family='Liberation Sans, Arial, Helvetica, sans-serif'"
out = []
add = out.append

add(f"<text x='{X0}' y='40' {FONT} font-size='26' font-weight='bold' fill='{TITLE}'>BLOSUM62 substitution scores</text>")

for i, a in enumerate(ORDER):
    add(f"<text x='{X0 - 20}' y='{Y0 + i*C + C/2 + 7}' {FONT} font-size='20' font-weight='bold' fill='{COL[a]}' text-anchor='middle'>{a}</text>")
    add(f"<text x='{X0 + i*C + C/2}' y='{Y0 - 13}' {FONT} font-size='20' font-weight='bold' fill='{COL[a]}' text-anchor='middle'>{a}</text>")
    for j, b in enumerate(ORDER[:i + 1]):
        v = M[a, b]
        same = COL[a] == COL[b]
        face = COL[a] if same else ("#F2F2F2" if v <= 0 else "#DDE7F0")
        opacity = " fill-opacity='0.9'" if same else ""
        x, y = X0 + j * C, Y0 + i * C
        add(f"<rect x='{x}' y='{y}' width='{C}' height='{C}' fill='{face}'{opacity} stroke='#ffffff' stroke-width='1.5'/>")
        fill = "#ffffff" if same else (TEXT if v > 0 else GREY)
        weight = "bold" if v > 0 else "normal"
        add(f"<text x='{x + C/2}' y='{y + C/2 + 6}' {FONT} font-size='16' font-weight='{weight}' fill='{fill}' text-anchor='middle'>{v}</text>")

# legend
ly = Y0 + 10
for g, col, lab in GROUPS:
    add(f"<rect x='{LX}' y='{ly}' width='30' height='30' fill='{col}'/>")
    add(f"<text x='{LX + 42}' y='{ly + 22}' {FONT} font-size='19' fill='{TEXT}'>{lab}  ({' '.join(g)})</text>")
    ly += 50
ly += 30
for k, line in enumerate(["coloured cell = same group", "light blue = positive score", "across groups"]):
    add(f"<text x='{LX}' y='{ly + k*24}' {FONT} font-size='17' fill='{GREY}'>{line}</text>")

svg = (f"<svg xmlns='http://www.w3.org/2000/svg' width='{W}' height='{H}' viewBox='0 0 {W} {H}'>\n"
       f"<rect width='{W}' height='{H}' fill='#ffffff'/>\n" + "\n".join(out) + "\n</svg>\n")
(Path(__file__).resolve().parent.parent / "fig" / "blosum62.svg").write_text(svg)
