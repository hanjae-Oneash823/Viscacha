"""Post-process the six remaining raw 3D renders into presentation-ready figures.

Same treatment as annotate_bace1_figures.py: bare PyMOL geometry gets cropped,
then a title, legend, optional callouts, and a footnote are added as a 2D
annotation pass. See that script for why labels are added this way rather
than as native PyMOL labels.
"""

from __future__ import annotations

from PIL import Image, ImageDraw, ImageFont

SCRATCH = "/tmp/claude-1813/-home-welcome3-Viscacha-pipeline/7f2faea9-593c-4a2e-ae4d-ecc12e47472b/scratchpad"
EXPANDED = "/home/welcome3/Viscacha_pipeline/outputs/docking_campaign/figures/expanded_campaign"
BC = "/home/welcome3/Viscacha_pipeline/outputs/docking_campaign/figures/bc_candidates"

FONT_BOLD = "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"
FONT_REG = "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"

INK = (35, 38, 42)
MUTED = (110, 118, 128)
RULE = (210, 214, 219)


def autocrop(im: Image.Image, pad: int = 30) -> tuple[Image.Image, int, int]:
    gray = im.convert("L")
    bbox = Image.eval(gray, lambda p: 0 if p > 250 else 255).getbbox()
    if bbox is None:
        return im, 0, 0
    left = max(bbox[0] - pad, 0)
    top = max(bbox[1] - pad, 0)
    right = min(bbox[2] + pad, im.width)
    bottom = min(bbox[3] + pad, im.height)
    return im.crop((left, top, right, bottom)), left, top


def draw_callout(draw, anchor, text_pos, text, font, dot_color, align=None):
    ax, ay = anchor
    tx, ty = text_pos
    bbox = draw.multiline_textbbox((0, 0), text, font=font, spacing=8)
    text_w = bbox[2] - bbox[0]
    text_h = bbox[3] - bbox[1]
    anchor_right = (tx > ax) if align is None else (align == "right")
    tx_draw = tx + 14 if anchor_right else tx - 14 - text_w
    ty_draw = ty - text_h / 2

    draw.line([(ax, ay), (tx, ty)], fill=MUTED, width=3)
    r = 9
    draw.ellipse([ax - r, ay - r, ax + r, ay + r], fill=dot_color, outline=(255, 255, 255), width=2)
    draw.rectangle(
        [tx_draw - 12, ty_draw - 10, tx_draw + text_w + 12, ty_draw + text_h + 10],
        fill=(255, 255, 255), outline=RULE, width=2,
    )
    draw.multiline_text((tx_draw, ty_draw), text, fill=INK, font=font, spacing=8)


def build_figure(raw_path, title, subtitle, legend, footer, out_path, callouts=None, canvas_w=2000):
    im = Image.open(raw_path).convert("RGB")
    im, off_x, off_y = autocrop(im, pad=30)

    scale = canvas_w / im.width
    canvas_h = int(im.height * scale)
    im = im.resize((canvas_w, canvas_h), Image.LANCZOS)

    title_h, legend_h, footer_h = 130, 90, 70
    total_h = title_h + canvas_h + legend_h + footer_h
    page = Image.new("RGB", (canvas_w, total_h), (255, 255, 255))
    draw = ImageDraw.Draw(page)

    f_title = ImageFont.truetype(FONT_BOLD, 40)
    f_sub = ImageFont.truetype(FONT_REG, 24)
    f_label = ImageFont.truetype(FONT_BOLD, 25)
    f_legend = ImageFont.truetype(FONT_REG, 24)
    f_footer = ImageFont.truetype(FONT_REG, 22)

    draw.text((40, 28), title, fill=INK, font=f_title)
    draw.text((40, 82), subtitle, fill=MUTED, font=f_sub)
    draw.line([(0, title_h - 1), (canvas_w, title_h - 1)], fill=RULE, width=2)

    page.paste(im, (0, title_h))

    def to_canvas(pt_raw):
        x = (pt_raw[0] - off_x) * scale
        y = (pt_raw[1] - off_y) * scale + title_h
        return x, y

    for callout in callouts or []:
        raw_anchor, text_pos, text, color = callout[:4]
        align = callout[4] if len(callout) > 4 else None
        draw_callout(draw, to_canvas(raw_anchor), (text_pos[0], text_pos[1] + title_h), text, f_label, color, align=align)

    legend_y0 = title_h + canvas_h
    draw.line([(0, legend_y0), (canvas_w, legend_y0)], fill=RULE, width=2)
    lx = 40
    ly = legend_y0 + legend_h // 2
    for color, label in legend:
        draw.rectangle([lx, ly - 14, lx + 28, ly + 14], fill=color, outline=(255, 255, 255))
        lx += 40
        draw.text((lx, ly - 13), label, fill=INK, font=f_legend)
        lx += int(draw.textlength(label, font=f_legend)) + 36

    footer_y0 = legend_y0 + legend_h
    draw.line([(0, footer_y0), (canvas_w, footer_y0)], fill=RULE, width=2)
    draw.text((40, footer_y0 + 20), footer, fill=MUTED, font=f_footer)

    page.save(out_path)
    print("wrote", out_path, page.size)


build_figure(
    raw_path=f"{SCRATCH}/3D_BACE1_variant_pocket_overlay_raw.png",
    title="Both BACE1 deletion isoforms displace verubecestat from the crystal pose",
    subtitle="Canonical redocking (blue) matches the crystal pose (gray); BACE1-476 (teal) and BACE1-457 (orange) land elsewhere",
    legend=[
        ((60, 60, 60), "Crystal pose (5HU1)"),
        ((30, 90, 200), "Canonical redocked"),
        ((0, 150, 140), "BACE1-476"),
        ((230, 130, 30), "BACE1-457"),
    ],
    footer="Mean top-pose RMSD to crystal (5 seeds each): canonical 0.96 A; BACE1-476 7.12 A; BACE1-457 7.01 A.",
    out_path=f"{EXPANDED}/3D_BACE1_variant_pocket_overlay.png",
)

build_figure(
    raw_path=f"{SCRATCH}/3D_CHRNA7_canonical_encenicline_raw.png",
    title="Canonical alpha7 recovers the encenicline crystal pose",
    subtitle="7EKP homopentamer; the docked pose sits at the A/B intersubunit site, matching the crystal ligand",
    legend=[
        ((160, 160, 160), "Other subunits"),
        ((90, 150, 220), "Chain A"),
        ((190, 225, 225), "Chain B"),
        ((200, 30, 180), "Docked encenicline"),
    ],
    footer="Mean top-pose RMSD to crystal: 0.297 A across 6 independent runs -- the tightest pose recovery in the campaign.",
    out_path=f"{EXPANDED}/3D_CHRNA7_canonical_encenicline.png",
    callouts=[((894, 726), (1550, 90), "Docked encenicline\n(matches crystal pose)", (200, 30, 180), "right")],
)

build_figure(
    raw_path=f"{SCRATCH}/3D_CHRFAM7A_B_face_encenicline_raw.png",
    title="A CHRFAM7A fusion at the B face displaces encenicline",
    subtitle="Same intersubunit site as canonical alpha7, with the fusion domain (violet) replacing the B-face subunit",
    legend=[
        ((160, 160, 160), "Other subunits"),
        ((90, 150, 220), "Chain A"),
        ((205, 140, 210), "CHRFAM7A fusion (chain B)"),
        ((230, 130, 30), "Docked encenicline"),
    ],
    footer="Mean pose displacement: 8.45 +/- 0.15 A across 6 runs. Topology hypothesis, not a confirmed biological assembly.",
    out_path=f"{EXPANDED}/3D_CHRFAM7A_B_face_encenicline.png",
    callouts=[((1404, 708), (1550, 90), "Docked encenicline\n(~8.5 A from crystal pose)", (230, 130, 30), "right")],
)

build_figure(
    raw_path=f"{SCRATCH}/3D_CHRNA7_topology_site_overlay_raw.png",
    title="Ligand pose depends on which face carries the CHRFAM7A fusion",
    subtitle="Pocket-lining residues for canonical (blue) and the B-face hybrid (violet) around their respective docked poses",
    legend=[
        ((60, 60, 60), "Crystal pose"),
        ((200, 30, 180), "Canonical docked pose"),
        ((230, 130, 30), "B-face hybrid docked pose"),
        ((90, 150, 220), "Canonical / hybrid pocket residues"),
    ],
    footer="Canonical and A-face poses stay near the crystal pose; the B-face hybrid is displaced ~8.5 A (see canonical_crossdock_stability).",
    out_path=f"{EXPANDED}/3D_CHRNA7_topology_site_overlay.png",
)

build_figure(
    raw_path=f"{SCRATCH}/C_CACNA1D_full_structure_overlay_raw.png",
    title="CACNA1D-214 retains the isradipine pocket but loses the distal C-terminus",
    subtitle="Canonical CaV1.3 (gray/teal shared fold) vs. the AD-associated 1625-aa isoform",
    legend=[
        ((170, 170, 170), "Canonical-only fold"),
        ((10, 110, 100), "Shared fold"),
        ((225, 75, 65), "Deleted in CACNA1D-214"),
        ((108, 43, 217), "Canonical pocket residues"),
        ((242, 177, 52), "Alternate pocket residues"),
    ],
    footer="All 14 mapped isradipine-contact residues are retained and sequence-identical; only the distal C-terminus (1606-2161) is lost.",
    out_path=f"{BC}/C_CACNA1D_full_structure_overlay.png",
)

build_figure(
    raw_path=f"{SCRATCH}/C_CACNA1D_retained_pocket_overlay_raw.png",
    title="The 14 isradipine-contact residues are retained and spatially matched",
    subtitle="Pocket residues modeled independently for canonical (purple) and CACNA1D-214 (gold), then superimposed",
    legend=[
        ((108, 43, 217), "Canonical pocket residues"),
        ((242, 177, 52), "CACNA1D-214 pocket residues"),
        ((200, 200, 200), "Backbone context"),
    ],
    footer="Local pocket C-alpha RMSD 1.49 A across all 14 residues -- the isoform difference is distal, not at the ligand site.",
    out_path=f"{BC}/C_CACNA1D_retained_pocket_overlay.png",
)
