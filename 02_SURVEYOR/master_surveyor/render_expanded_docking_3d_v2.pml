set max_threads, 16
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
show cartoon, bace_can
color gray80, bace_can
set cartoon_transparency, 0.82, bace_can
show sticks, veru_xtal or veru_can or veru476 or veru457
color gray30, veru_xtal
color marine, veru_can
color teal, veru476
color tv_orange, veru457
set stick_radius, 0.30, veru_xtal or veru_can or veru476 or veru457
show spheres, (veru_xtal or veru_can or veru476 or veru457) and elem N+O+F
set sphere_scale, 0.24, (veru_xtal or veru_can or veru476 or veru457) and elem N+O+F
orient veru_xtal
zoom veru_xtal, 9
turn x, -12
turn y, 18
ray 2400, 1800
png /tmp/claude-1813/-home-welcome3-Viscacha-pipeline/7f2faea9-593c-4a2e-ae4d-ecc12e47472b/scratchpad/3D_BACE1_variant_pocket_overlay_raw.png, dpi=300

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
color gray80, alpha7
color skyblue, alpha7 and chain A
color palecyan, alpha7 and chain B
set cartoon_transparency, 0.35, alpha7
show sticks, ence_xtal or ence_can
color gray20, ence_xtal
color magenta, ence_can
set stick_radius, 0.42, ence_xtal or ence_can
show spheres, (ence_xtal or ence_can) and elem N+O
set sphere_scale, 0.32, (ence_xtal or ence_can) and elem N+O
orient alpha7
turn x, -8
turn y, 14
zoom alpha7, 4
ray 2400, 1800
png /tmp/claude-1813/-home-welcome3-Viscacha-pipeline/7f2faea9-593c-4a2e-ae4d-ecc12e47472b/scratchpad/3D_CHRNA7_canonical_encenicline_raw.png, dpi=300

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
color gray80, hybridB
color skyblue, hybridB and chain A
color violet, hybridB and chain B
set cartoon_transparency, 0.35, hybridB
show sticks, ence_xtal or ence_hybB
color gray20, ence_xtal
color tv_orange, ence_hybB
set stick_radius, 0.42, ence_xtal or ence_hybB
show spheres, (ence_xtal or ence_hybB) and elem N+O
set sphere_scale, 0.32, (ence_xtal or ence_hybB) and elem N+O
orient hybridB
turn x, -8
turn y, 14
zoom hybridB, 4
ray 2400, 1800
png /tmp/claude-1813/-home-welcome3-Viscacha-pipeline/7f2faea9-593c-4a2e-ae4d-ecc12e47472b/scratchpad/3D_CHRFAM7A_B_face_encenicline_raw.png, dpi=300

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
show cartoon, alpha7 or hybridB
color gray85, alpha7
color gray85, hybridB
set cartoon_transparency, 0.88, alpha7
set cartoon_transparency, 0.88, hybridB
show sticks, site_can or site_hyb
color skyblue, site_can
color violet, site_hyb
set stick_radius, 0.16, site_can or site_hyb
show sticks, ence_xtal or ence_can or ence_hybB
color gray20, ence_xtal
color magenta, ence_can
color tv_orange, ence_hybB
set stick_radius, 0.38, ence_xtal or ence_can or ence_hybB
orient ence_xtal
zoom ence_xtal or ence_can or ence_hybB, 7
turn x, -10
turn y, 22
ray 2400, 1800
png /tmp/claude-1813/-home-welcome3-Viscacha-pipeline/7f2faea9-593c-4a2e-ae4d-ecc12e47472b/scratchpad/3D_CHRNA7_topology_site_overlay_raw.png, dpi=300

quit
