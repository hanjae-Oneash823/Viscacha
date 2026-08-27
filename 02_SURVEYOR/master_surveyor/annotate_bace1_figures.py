"""Post-process the raw PyMOL BACE1 renders into presentation-ready figures.

PyMOL's native atom labels collide when features sit close together on
screen (see the v1 render). Cleaner approach used everywhere in structural
biology figure-making: render bare geometry in PyMOL, annotate with a proper
2D text/graphics layer afterward. This script crops the raw renders, adds a
title bar, a color legend, and callout leader lines to fixed screen anchors
measured on the rendered image.
"""

from __future__ import annotations

from PIL import Image, ImageDraw, ImageFont

SCRATCH = "/tmp/claude-1813/-home-welcome3-Viscacha-pipeline/7f2faea9-593c-4a2e-ae4d-ecc12e47472b/scratchpad"
OUT_DIR = "/home/welcome3/Viscacha_pipeline/outputs/docking_campaign/figures/bc_candidates"

FONT_BOLD = "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"
FONT_REG = "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"

INK = (35, 38, 42)
MUTED = (110, 118, 128)
RULE = (210, 214, 219)


def autocrop(im: Image.Image, pad: int = 40) -> tuple[Image.Image, int, int]:
    gray = im.convert("L")
    bbox = Image.eval(gray, lambda p: 0 if p > 250 else 255).getbbox()
    if bbox is None:
        return im, 0, 0
    left = max(bbox[0] - pad, 0)
    top = max(bbox[1] - pad, 0)
    right = min(bbox[2] + pad, im.width)
    bottom = min(bbox[3] + pad, im.height)
    return im.crop((left, top, right, bottom)), left, top


def draw_callout(draw, anchor, text_pos, text, font, dot_color, canvas_w, align=None):
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


def build_figure(raw_path, title, subtitle, callouts, legend, footer, out_path, canvas_w=2000):
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
    f_label = ImageFont.truetype(FONT_BOLD, 26)
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

    for callout in callouts:
        raw_anchor, text_pos, text, color = callout[:4]
        align = callout[4] if len(callout) > 4 else None
        draw_callout(draw, to_canvas(raw_anchor), (text_pos[0], text_pos[1] + title_h), text, f_label, color, canvas_w, align=align)

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


# ---------------------------------------------------------------------------
# Figure 1: domain-loss overview
# ---------------------------------------------------------------------------
build_figure(
    raw_path=f"{SCRATCH}/bace1_overview_raw.png",
    title="BACE1-202 deletes part of the verubecestat pocket",
    subtitle="Canonical BACE1 (gray/teal shared fold) vs. the 401-aa AD-associated isoform BACE1-202",
    callouts=[
        ((714, 552), (60, 120), "Residues 57-120: present in\ncanonical, absent in BACE1-202\n(12/28 drug-contact residues)", (224, 72, 62), "right"),
        ((1146, 774), (1550, 140), "Verubecestat\n(crystal pose, PDB 5HU1)", (244, 185, 66)),
        ((1038, 882), (1550, 320), "Asp93 - catalytic residue\nlost in BACE1-202", (176, 32, 42)),
        ((1182, 954), (1550, 460), "Asp289 - catalytic residue\nretained (= Asp189 in alt)", (27, 79, 145)),
    ],
    legend=[
        ((141, 150, 161), "Canonical-only fold"),
        ((15, 155, 142), "Shared fold (canonical + alt)"),
        ((224, 72, 62), "Deleted in BACE1-202"),
        ((244, 185, 66), "Verubecestat"),
    ],
    footer="Structures aligned to the 5HU1 coordinate frame. Deleted region shown on the canonical backbone; BACE1-202 has no residues there at all.",
    out_path=f"{OUT_DIR}/B_BACE1_isoform_overlay.png",
)

# ---------------------------------------------------------------------------
# Figure 2: active-site close-up
# ---------------------------------------------------------------------------
build_figure(
    raw_path=f"{SCRATCH}/bace1_closeup_raw.png",
    title="Same docking protocol, very different outcome",
    subtitle="Matched Vina docking, 5 seeds each: canonical recovers the crystal pose, BACE1-202 does not",
    callouts=[
        ((1392, 840), (1550, 130), "Canonical pose\n0.96 A RMSD to crystal", (244, 185, 66)),
        ((948, 912), (120, 130), "BACE1-202 pose\n10.6 A RMSD to crystal", (142, 68, 173), "right"),
        ((1242, 1056), (1550, 320), "Asp93 side chain\n(absent in BACE1-202)", (176, 32, 42)),
    ],
    legend=[
        ((244, 185, 66), "Canonical top pose"),
        ((142, 68, 173), "BACE1-202 top pose"),
        ((176, 32, 42), "Asp93 (canonical only)"),
        ((68, 68, 68), "Pose displacement"),
    ],
    footer="GNINA agrees: CNNscore 0.92 (canonical) vs 0.22 (alt); CNNaffinity 7.37 vs 5.52 -- consistent across all 5 seeds.",
    out_path=f"{OUT_DIR}/B_BACE1_pose_displacement.png",
)
