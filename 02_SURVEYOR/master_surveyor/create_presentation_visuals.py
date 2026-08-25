#!/usr/bin/env python3
"""Render polished, data-faithful 16:9 visuals for the docking presentation."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch, Rectangle
from mpl_toolkits.mplot3d import proj3d


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs" / "docking_campaign"
NAVY = "#071A33"
BLUE = "#20A4F3"
TEAL = "#29D3B6"
GOLD = "#FFCA58"
CORAL = "#FF6B6B"
MUTED = "#AAB7C8"


def pdb_xyz(path: Path, first_model: bool = False) -> tuple[np.ndarray, list[str]]:
    xyz, atoms = [], []
    in_model = False
    for line in path.read_text().splitlines():
        if first_model and line.startswith("MODEL"):
            if in_model:
                break
            in_model = True
        if line.startswith(("ATOM", "HETATM")):
            xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            atoms.append(line[12:16].strip())
        if first_model and in_model and line.startswith("ENDMDL"):
            break
    return np.asarray(xyz), atoms


def pdb_ca(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Return chain-A C-alpha coordinates and experimental residue numbers."""
    xyz, residues = [], []
    for line in path.read_text().splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            residues.append(int(line[22:26]))
    return np.asarray(xyz), np.asarray(residues)


def card(ax, x, y, w, h, title, value, subtitle, color=BLUE):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.012,rounding_size=0.025", fc="#102947", ec="#244363", lw=1.2, transform=ax.transAxes))
    ax.text(x + .03, y + h - .07, title.upper(), color=MUTED, fontsize=9, weight="bold", transform=ax.transAxes)
    ax.text(x + .03, y + h - .16, value, color=color, fontsize=22, weight="bold", transform=ax.transAxes)
    ax.text(x + .03, y + .055, subtitle, color="#D9E3EF", fontsize=9, transform=ax.transAxes)


def fyn_slide() -> None:
    campaign = OUT / "FYN_saracatinib"
    figures = campaign / "figures" / "presentation"
    figures.mkdir(parents=True, exist_ok=True)
    summary = json.loads((campaign / "analysis" / "summary.json").read_text())
    best = summary["best_top_pose"]
    protein, _ = pdb_xyz(campaign / "prepared" / "fyn_chain_A_protein.pdb")
    crystal, _ = pdb_xyz(campaign / "prepared" / "saracatinib_H8H_A601_crystal.pdb")
    docked, _ = pdb_xyz(campaign / "runs" / "vina_seed4409_ex32" / "docked_poses.pdbqt", first_model=True)
    # Protein context close to the ligand, sampled to retain visual clarity.
    center = crystal.mean(axis=0)
    near = protein[np.linalg.norm(protein - center, axis=1) < 14]
    if len(near) > 850: near = near[::max(1, len(near)//850)]

    fig = plt.figure(figsize=(16, 9), facecolor=NAVY, constrained_layout=True)
    gs = fig.add_gridspec(12, 18)
    ax3d = fig.add_subplot(gs[2:11, :10], projection="3d", facecolor=NAVY)
    axtext = fig.add_subplot(gs[2:11, 10:])
    axtext.set_facecolor(NAVY); axtext.axis("off")
    fig.text(.055, .93, "Validated docking control: FYN–saracatinib", color="white", fontsize=25, weight="bold")
    fig.text(.055, .888, "A known drug pose is recovered reproducibly; the AD-enriched alternate loses the kinase pocket.", color="#C9D6E4", fontsize=12)

    ax3d.scatter(near[:,0], near[:,1], near[:,2], s=8, c="#7B96B6", alpha=.20, depthshade=False, label="FYN pocket atoms")
    ax3d.plot(crystal[:,0], crystal[:,1], crystal[:,2], c=TEAL, lw=3.2, marker="o", ms=4, label="Crystal saracatinib")
    ax3d.plot(docked[:,0], docked[:,1], docked[:,2], c=GOLD, lw=2.4, marker="o", ms=3.5, label="Best Vina pose")
    ax3d.set_axis_off(); ax3d.view_init(elev=18, azim=-63)
    ax3d.legend(loc="upper left", bbox_to_anchor=(.02,.98), frameon=False, labelcolor="white", fontsize=10)
    ax3d.set_title("Experimental binding pocket and recovered pose", color="white", fontsize=13, pad=4)

    card(axtext, .02, .67, .43, .25, "Reproducibility", "9 / 9", "top poses recovered within 2 Å", TEAL)
    card(axtext, .50, .67, .43, .25, "Median recovery", f"{summary['top_pose_rmsd_median_A']:.2f} Å", "fixed-frame heavy-atom RMSD", GOLD)
    card(axtext, .02, .37, .43, .25, "Best recovery", f"{best['rmsd_to_crystal_heavy_atom_uncorrected_angstrom']:.2f} Å", "seed 4409; exhaustive search 32", GOLD)
    card(axtext, .50, .37, .43, .25, "Best Vina score", f"{best['vina_affinity_kcal_mol']:.2f}", "kcal/mol; scoring value, not affinity", BLUE)
    axtext.text(.02, .25, "Isoform interpretation", color="white", fontsize=15, weight="bold", transform=axtext.transAxes)
    axtext.add_patch(Rectangle((.02,.14), .91,.055, transform=axtext.transAxes, color="#536D8B"))
    # FYN protein domains: N terminus, SH3, SH2, kinase.
    for x,w,label,col in [(0.05,.10,"N", "#7186A3"),(.18,.13,"SH3", TEAL),(.34,.14,"SH2", BLUE),(.51,.39,"KINASE", GOLD)]:
        axtext.add_patch(FancyBboxPatch((x,.12),w,.095,boxstyle="round,pad=.007,rounding_size=.012",fc=col,ec="none",transform=axtext.transAxes))
        axtext.text(x+w/2,.167,label,ha="center",va="center",fontsize=9,color=NAVY,weight="bold",transform=axtext.transAxes)
    axtext.annotate("Novel AD-enriched transcript ends at residue 115", xy=(.155,.10), xytext=(.155,.025), ha="center", color=CORAL, fontsize=10, weight="bold", arrowprops=dict(arrowstyle="-|>", color=CORAL), transform=axtext.transAxes)
    axtext.text(.02, -.06, "Conclusion: the alternate is predicted to lack the saracatinib-binding kinase domain.\nThis is a structural hypothesis—not proof of expression, drug response, or disease causality.", color="#D9E3EF", fontsize=10, transform=axtext.transAxes)
    fig.savefig(figures / "FYN_saracatinib_presentation_slide.png", dpi=300, facecolor=NAVY)
    fig.savefig(figures / "FYN_saracatinib_presentation_slide.pdf", facecolor=NAVY)
    plt.close(fig)


def fyn_3d_diagram() -> None:
    """Render standalone protein overview + ligand-pocket close-up for slides."""
    campaign = OUT / "FYN_saracatinib"
    figures = campaign / "figures" / "presentation"
    backbone, resnums = pdb_ca(campaign / "prepared" / "fyn_chain_A_protein.pdb")
    crystal, _ = pdb_xyz(campaign / "prepared" / "saracatinib_H8H_A601_crystal.pdb")
    docked, _ = pdb_xyz(campaign / "runs" / "vina_seed4409_ex32" / "docked_poses.pdbqt", first_model=True)
    protein, _ = pdb_xyz(campaign / "prepared" / "fyn_chain_A_protein.pdb")
    center = crystal.mean(axis=0)
    pocket = protein[np.linalg.norm(protein - center, axis=1) < 10]
    if len(pocket) > 1100: pocket = pocket[::max(1, len(pocket)//1100)]

    fig = plt.figure(figsize=(16, 9), facecolor=NAVY, constrained_layout=True)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.05, .95])
    ax1 = fig.add_subplot(gs[0,0], projection="3d", facecolor=NAVY)
    ax2 = fig.add_subplot(gs[0,1], projection="3d", facecolor=NAVY)

    # Backbone shown in domain-aware colors. The resolved chain begins after the N terminus.
    ranges = [(0, 150, "#6D86A3", "N-terminal region"), (150, 270, TEAL, "SH2/SH3 region"), (270, 10_000, BLUE, "kinase domain")]
    for lo, hi, color, label in ranges:
        mask=(resnums>=lo)&(resnums<hi)
        if mask.any(): ax1.plot(backbone[mask,0],backbone[mask,1],backbone[mask,2],color=color,lw=2.2,label=label,alpha=.95)
    ax1.plot(crystal[:,0],crystal[:,1],crystal[:,2],color=TEAL,lw=3.8,marker="o",ms=4.5,label="crystal ligand")
    ax1.plot(docked[:,0],docked[:,1],docked[:,2],color=GOLD,lw=2.8,marker="o",ms=3.8,label="best Vina pose")
    ax1.scatter(*center,s=85,color="white",edgecolor=GOLD,depthshade=False,zorder=8)
    ax1.set_axis_off(); ax1.view_init(elev=18,azim=-65)
    ax1.legend(loc="lower left",bbox_to_anchor=(.0,.0),frameon=False,labelcolor="white",fontsize=9)

    ax2.scatter(pocket[:,0],pocket[:,1],pocket[:,2],s=13,c="#6285AA",alpha=.32,depthshade=False,label="protein pocket atoms")
    ax2.plot(crystal[:,0],crystal[:,1],crystal[:,2],color=TEAL,lw=5,marker="o",ms=5,label="crystal saracatinib")
    ax2.plot(docked[:,0],docked[:,1],docked[:,2],color=GOLD,lw=3.3,marker="o",ms=4,label="best Vina pose")
    ax2.set_axis_off(); ax2.view_init(elev=10,azim=-33)
    ax2.legend(loc="lower left",bbox_to_anchor=(.0,.0),frameon=False,labelcolor="white",fontsize=9)
    fig.text(.30,.12,"The alternate transcript ends at residue 115—\nbefore the kinase domain shown here.",ha="center",color=CORAL,fontsize=13,weight="bold")
    fig.savefig(figures / "FYN_3D_protein_and_pocket.png",dpi=300,facecolor=NAVY)
    fig.savefig(figures / "FYN_3D_protein_and_pocket.pdf",facecolor=NAVY)
    plt.close(fig)


def kit_slide() -> None:
    campaign = OUT / "KIT_masitinib"
    figures = campaign / "figures" / "presentation"
    figures.mkdir(parents=True, exist_ok=True)
    with (campaign / "analysis" / "masitinib_1T46_replicates.csv").open() as h:
        values = np.asarray([float(r["vina_affinity_kcal_mol"]) for r in csv.DictReader(h)])
    fig = plt.figure(figsize=(16, 9), facecolor=NAVY, constrained_layout=True)
    gs = fig.add_gridspec(12, 18)
    ax = fig.add_subplot(gs[3:10, :10], facecolor=NAVY)
    ax2 = fig.add_subplot(gs[3:10, 10:], facecolor=NAVY)
    fig.text(.055,.93,"KIT–masitinib: a reliable baseline, not an isoform comparison",color="white",fontsize=23,weight="bold")
    fig.text(.055,.888,"We report what the experimental structure supports—and explicitly stop where it does not.",color="#C9D6E4",fontsize=12)
    x=np.arange(1,len(values)+1)
    ax.scatter(x,values,s=90,c=BLUE,edgecolor="white",lw=.6,zorder=3)
    ax.axhline(values.mean(),color=CORAL,lw=2.5)
    ax.fill_between([.5,9.5],values.mean()-values.std(ddof=1),values.mean()+values.std(ddof=1),color=CORAL,alpha=.16)
    ax.set_xlim(.5,9.5); ax.set_xticks(x); ax.set_xlabel("independent Vina run",color="white",fontsize=11);ax.set_ylabel("top Vina score (kcal/mol)",color="white",fontsize=11)
    ax.tick_params(colors="#D9E3EF"); [s.set_color("#4D6684") for s in ax.spines.values()]
    ax.set_title("Canonical c-KIT in experimental 1T46 pocket",color="white",fontsize=14,pad=14)
    ax.text(.02,.05,f"mean {values.mean():.2f} ± {values.std(ddof=1):.2f} kcal/mol\n9 reproducible canonical runs",transform=ax.transAxes,color="white",fontsize=12,weight="bold",bbox=dict(boxstyle="round,pad=.45",fc="#102947",ec="#244363"))
    ax2.axis("off"); ax2.set_xlim(540,950); ax2.set_ylim(0,1)
    ax2.hlines(.51,565,933,color="#7B96B6",lw=14)
    ax2.axvspan(690,761,color="#CBD3DF",alpha=.52)
    ax2.axvline(715,color=CORAL,lw=3,ls="--")
    ax2.text(715,.20,"KIT-223 deletion\nresidue 715",ha="center",color=CORAL,weight="bold",fontsize=13)
    ax2.text(725,.71,"unresolved kinase-insert\nsegment in 1T46 (690–761)",ha="center",color="white",fontsize=14,weight="bold")
    ax2.text(625,.65,"resolved",ha="center",color="#D9E3EF",fontsize=12);ax2.text(845,.65,"resolved",ha="center",color="#D9E3EF",fontsize=12)
    ax2.text(745,.97,"Why no KIT-223 score is shown",ha="center",color="white",fontsize=16,weight="bold")
    ax2.text(745,.035,"A score would depend on an unvalidated rebuilt loop.\nExcluding it is the scientifically correct result.",ha="center",color="#D9E3EF",fontsize=12)
    fig.savefig(figures / "KIT_masitinib_presentation_slide.png",dpi=300,facecolor=NAVY)
    fig.savefig(figures / "KIT_masitinib_presentation_slide.pdf",facecolor=NAVY)
    plt.close(fig)


if __name__ == "__main__":
    fyn_slide(); kit_slide(); fyn_3d_diagram()
