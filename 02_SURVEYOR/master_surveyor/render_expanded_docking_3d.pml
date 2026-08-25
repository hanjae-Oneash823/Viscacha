# High-resolution 3D figures for the expanded docking campaign.
set ray_opaque_background, off
set antialias, 2
set ray_trace_mode, 1
set ray_shadows, on
set ambient, 0.38
set direct, 0.62
set specular, 0.25
set shininess, 28
set cartoon_fancy_helices, on
set cartoon_smooth_loops, on
bg_color white

# ---------------------------------------------------------------------------
# BACE1: canonical and two deletion-model pockets with top verubecestat poses.
reinitialize
bg_color white
set ray_opaque_background, on
set antialias, 2
set ray_trace_mode, 1
load /home/welcome3/Viscacha_pipeline/outputs/docking_campaign/systems/BACE1_verubecestat/prepared/bace1_5HU1_chain_A_protein.pdb, bace_can
load /home/welcome3/Viscacha_pipeline/outputs/docking_campaign/systems/BACE1_verubecestat/prepared/bace1_476_alphafold_aligned_to_5HU1.pdb, bace476
load /home/welcome3/Viscacha_pipeline/outputs/docking_campaign/systems/BACE1_verubecestat/prepared/bace1_457_alphafold_aligned_to_5HU1.pdb, bace457
load /home/welcome3/Viscacha_pipeline/outputs/docking_campaign/systems/BACE1_verubecestat/prepared/verubecestat_66F_A501_crystal.pdb, veru_xtal
load /home/welcome3/Viscacha_pipeline/outputs/docking_campaign/systems/BACE1_verubecestat/runs/canonical_obabel_seed1103_ex32/docked_poses.pdbqt, veru_can
load /home/welcome3/Viscacha_pipeline/outputs/docking_campaign/systems/BACE1_verubecestat/runs/alternate_476_seed1103_ex32/docked_poses.pdbqt, veru476
load /home/welcome3/Viscacha_pipeline/outputs/docking_campaign/systems/BACE1_verubecestat/runs/alternate_457_seed1103_ex32/docked_poses.pdbqt, veru457
frame 1
hide everything
select pocket_can, byres (bace_can within 8 of veru_xtal)
select pocket476, byres (bace476 within 8 of veru_xtal)
select pocket457, byres (bace457 within 8 of veru_xtal)
show cartoon, bace_can
color gray70, bace_can
set cartoon_transparency, 0.72, bace_can
show sticks, pocket_can
color gray55, pocket_can
show lines, pocket476
color teal, pocket476
show lines, pocket457
color salmon, pocket457
show sticks, veru_xtal or veru_can or veru476 or veru457
color gray30, veru_xtal
color marine, veru_can
color teal, veru476
color tv_orange, veru457
set stick_radius, 0.22, veru_xtal or veru_can or veru476 or veru457
set stick_radius, 0.12, pocket_can or pocket476 or pocket457
show spheres, (veru_xtal or veru_can or veru476 or veru457) and elem N+O+F
set sphere_scale, 0.22, (veru_xtal or veru_can or veru476 or veru457) and elem N+O+F
orient veru_xtal
zoom veru_xtal, 11
turn x, -12
turn y, 18
ray 2400, 1800
png /home/welcome3/Viscacha_pipeline/outputs/docking_campaign/figures/expanded_campaign/3D_BACE1_variant_pocket_overlay.png, dpi=300

# ---------------------------------------------------------------------------
# Canonical alpha7 encenicline site.
reinitialize
bg_color white
set ray_opaque_background, on
set antialias, 2
set ray_trace_mode, 1
load /home/welcome3/Viscacha_pipeline/outputs/docking_campaign/systems/CHRNA7_encenicline/prepared/chrna7_7EKP_pentamer_protein.pdb, alpha7
load /home/welcome3/Viscacha_pipeline/outputs/docking_campaign/systems/CHRNA7_encenicline/prepared/encenicline_I33_A601_crystal.pdb, ence_xtal
load /home/welcome3/Viscacha_pipeline/outputs/docking_campaign/systems/CHRNA7_encenicline/runs/canonical_obabel_seed1103_ex32/docked_poses.pdbqt, ence_can
frame 1
hide everything
show cartoon, alpha7
color gray75, alpha7
color marine, alpha7 and chain A
color teal, alpha7 and chain B
show surface, byres (alpha7 within 6 of ence_xtal)
set transparency, 0.42, byres (alpha7 within 6 of ence_xtal)
color gray85, byres (alpha7 within 6 of ence_xtal)
show sticks, ence_xtal or ence_can
color gray25, ence_xtal
color magenta, ence_can
set stick_radius, 0.24, ence_xtal or ence_can
orient alpha7
turn x, -8
turn y, 14
zoom alpha7, 0
ray 2400, 1800
png /home/welcome3/Viscacha_pipeline/outputs/docking_campaign/figures/expanded_campaign/3D_CHRNA7_canonical_encenicline.png, dpi=300

# ---------------------------------------------------------------------------
# Mixed-receptor B-face hypothesis at the same site.
reinitialize
bg_color white
set ray_opaque_background, on
set antialias, 2
set ray_trace_mode, 1
load /home/welcome3/Viscacha_pipeline/outputs/docking_campaign/systems/CHRNA7_encenicline/prepared/chrna7_chrFam7a_fusion_at_B_face.pdb, hybridB
load /home/welcome3/Viscacha_pipeline/outputs/docking_campaign/systems/CHRNA7_encenicline/prepared/encenicline_I33_A601_crystal.pdb, ence_xtal
load /home/welcome3/Viscacha_pipeline/outputs/docking_campaign/systems/CHRNA7_encenicline/runs/hybrid_B_face_seed1103_ex32/docked_poses.pdbqt, ence_hybB
frame 1
hide everything
show cartoon, hybridB
color gray75, hybridB
color marine, hybridB and chain A
color purple, hybridB and chain B
show surface, byres (hybridB within 6 of ence_xtal)
set transparency, 0.42, byres (hybridB within 6 of ence_xtal)
color gray85, byres (hybridB within 6 of ence_xtal)
show sticks, ence_xtal or ence_hybB
color gray25, ence_xtal
color tv_orange, ence_hybB
set stick_radius, 0.24, ence_xtal or ence_hybB
orient hybridB
turn x, -8
turn y, 14
zoom hybridB, 0
ray 2400, 1800
png /home/welcome3/Viscacha_pipeline/outputs/docking_campaign/figures/expanded_campaign/3D_CHRFAM7A_B_face_encenicline.png, dpi=300

# ---------------------------------------------------------------------------
# Direct site overlay for the canonical and B-face hypothesis.
reinitialize
bg_color white
set ray_opaque_background, on
set antialias, 2
set ray_trace_mode, 1
load /home/welcome3/Viscacha_pipeline/outputs/docking_campaign/systems/CHRNA7_encenicline/prepared/chrna7_7EKP_pentamer_protein.pdb, alpha7
load /home/welcome3/Viscacha_pipeline/outputs/docking_campaign/systems/CHRNA7_encenicline/prepared/chrna7_chrFam7a_fusion_at_B_face.pdb, hybridB
load /home/welcome3/Viscacha_pipeline/outputs/docking_campaign/systems/CHRNA7_encenicline/prepared/encenicline_I33_A601_crystal.pdb, ence_xtal
load /home/welcome3/Viscacha_pipeline/outputs/docking_campaign/systems/CHRNA7_encenicline/runs/canonical_obabel_seed1103_ex32/docked_poses.pdbqt, ence_can
load /home/welcome3/Viscacha_pipeline/outputs/docking_campaign/systems/CHRNA7_encenicline/runs/hybrid_B_face_seed1103_ex32/docked_poses.pdbqt, ence_hybB
frame 1
hide everything
select site_can, byres (alpha7 within 9 of ence_xtal)
select site_hyb, byres (hybridB within 9 of ence_xtal)
show cartoon, site_can or site_hyb
color marine, site_can
color purple, site_hyb
set cartoon_transparency, 0.48, site_can
set cartoon_transparency, 0.60, site_hyb
show sticks, ence_xtal or ence_can or ence_hybB
color gray25, ence_xtal
color magenta, ence_can
color tv_orange, ence_hybB
set stick_radius, 0.24, ence_xtal or ence_can or ence_hybB
orient ence_xtal
zoom ence_xtal, 13
turn x, -10
turn y, 22
ray 2400, 1800
png /home/welcome3/Viscacha_pipeline/outputs/docking_campaign/figures/expanded_campaign/3D_CHRNA7_topology_site_overlay.png, dpi=300

quit
