"""Generate fig/aa_codes.svg (amino acid code table used in docs/sequence_alignment.md).

Usage: python3 scripts/make_aa_codes.py
PNG version: rsvg-convert -z 2 fig/aa_codes.svg -o fig/aa_codes.png
"""
from pathlib import Path

groups_left = [
    ("negatively charged", "#cfe5f7", "#1f5f99", [
        ("Aspartic acid", "Asp", "D"), ("Glutamic acid", "Glu", "E")]),
    ("positively charged", "#f8d3d3", "#a3282c", [
        ("Arginine", "Arg", "R"), ("Lysine", "Lys", "K"), ("Histidine", "His", "H")]),
    ("uncharged polar", "#fbefc4", "#8a6400", [
        ("Asparagine", "Asn", "N"), ("Glutamine", "Gln", "Q"), ("Serine", "Ser", "S"),
        ("Threonine", "Thr", "T"), ("Tyrosine", "Tyr", "Y")]),
]
groups_right = [
    ("nonpolar", "#dcefd0", "#3c6e1f", [
        ("Alanine", "Ala", "A"), ("Glycine", "Gly", "G"), ("Valine", "Val", "V"),
        ("Leucine", "Leu", "L"), ("Isoleucine", "Ile", "I"), ("Proline", "Pro", "P"),
        ("Phenylalanine", "Phe", "F"), ("Methionine", "Met", "M"),
        ("Tryptophan", "Trp", "W"), ("Cysteine", "Cys", "C")]),
]
special = [("B", "Asx", "Asp or Asn"), ("Z", "Glx", "Glu or Gln"),
           ("X", "Xaa", "any amino acid"), ("U", "Sec", "selenocysteine"),
           ("O", "Pyl", "pyrrolysine"), ("*", "", "stop codon"), ("-", "", "gap")]

W, ROW, TOP = 1000, 34, 70
COLW = 470
X0 = [20, 510]
FONT = "font-family='DejaVu Sans, Arial, Helvetica, sans-serif'"
ink, muted, accent = "#222", "#666", "#444"
out = []
add = out.append


def column(x, groups, polar_title):
    add(f"<text x='{x+12}' y='{TOP-14}' {FONT} font-size='15' font-weight='bold' fill='{muted}'>AMINO ACID</text>")
    add(f"<text x='{x+200}' y='{TOP-14}' {FONT} font-size='15' font-weight='bold' fill='{muted}' text-anchor='middle'>CODE</text>")
    add(f"<text x='{x+272}' y='{TOP-14}' {FONT} font-size='15' font-weight='bold' fill='{muted}'>SIDE CHAIN</text>")
    y = TOP
    for label, bg, fg, rows in groups:
        h = ROW * len(rows)
        add(f"<rect x='{x}' y='{y}' width='{COLW}' height='{h}' rx='6' fill='{bg}'/>")
        for i, (name, three, one) in enumerate(rows):
            ty = y + i * ROW + 23
            if i:
                add(f"<line x1='{x+8}' x2='{x+COLW-8}' y1='{y+i*ROW}' y2='{y+i*ROW}' stroke='#ffffff' stroke-width='1.5'/>")
            add(f"<text x='{x+12}' y='{ty}' {FONT} font-size='17' fill='{ink}'>{name}</text>")
            add(f"<text x='{x+165}' y='{ty}' {FONT} font-size='17' fill='{ink}'>{three}</text>")
            add(f"<text x='{x+240}' y='{ty}' {FONT} font-size='19' font-weight='bold' fill='{fg}' text-anchor='middle'>{one}</text>")
            add(f"<text x='{x+272}' y='{ty}' {FONT} font-size='16' fill='{fg}'>{label}</text>")
        y += h + 4
    # bracket with group title
    by = TOP + 10 * ROW + 8 + 22
    add(f"<path d='M{x+2},{by-12} v12 h{COLW-4} v-12' fill='none' stroke='{accent}' stroke-width='1.5'/>")
    tw = 11.5 * len(polar_title) + 30
    add(f"<rect x='{x+COLW/2-tw/2}' y='{by-10}' width='{tw}' height='20' fill='#ffffff'/>")
    add(f"<text x='{x+COLW/2}' y='{by+6}' {FONT} font-size='16' font-weight='bold' fill='{accent}' text-anchor='middle' letter-spacing='1'>{polar_title}</text>")
    return by


column(X0[0], groups_left, "POLAR AMINO ACIDS")
by = column(X0[1], groups_right, "NONPOLAR AMINO ACIDS")

# special / ambiguity codes
sy = by + 40
add(f"<text x='20' y='{sy}' {FONT} font-size='15' font-weight='bold' fill='{muted}'>OTHER CODES USED IN SEQUENCES AND ALIGNMENTS</text>")
sy += 12
per_row = 4
cw = (W - 40) / per_row
for i, (one, three, desc) in enumerate(special):
    cx = 20 + (i % per_row) * cw
    cy = sy + (i // per_row) * ROW
    add(f"<rect x='{cx}' y='{cy}' width='{cw-8}' height='{ROW-6}' rx='5' fill='#eeeeee'/>")
    add(f"<text x='{cx+18}' y='{cy+20}' {FONT} font-size='18' font-weight='bold' fill='{ink}' text-anchor='middle'>{one}</text>")
    t = f"{three} – {desc}" if three else desc
    add(f"<text x='{cx+38}' y='{cy+20}' {FONT} font-size='15' fill='{ink}'>{t}</text>")
H = sy + 2 * ROW + 10

svg = (f"<svg xmlns='http://www.w3.org/2000/svg' width='{W}' height='{H}' viewBox='0 0 {W} {H}'>\n"
       f"<rect width='{W}' height='{H}' fill='#ffffff'/>\n" + "\n".join(out) + "\n</svg>\n")
(Path(__file__).resolve().parent.parent / "fig" / "aa_codes.svg").write_text(svg)
