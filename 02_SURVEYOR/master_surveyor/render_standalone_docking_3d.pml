# High-resolution standalone molecular renders from actual docking coordinates.

cd /home/welcome3/Viscacha_pipeline
reinitialize
set ray_opaque_background, off
set antialias, 2
set ray_trace_mode, 1
set ray_shadow, 0
set specular, 0.16
set shininess, 20
set depth_cue, 0
set orthoscopic, 1
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
set cartoon_sampling, 14
set stick_radius, 0.19
set sphere_scale, 0.24
set transparency_mode, 1
bg_color white

# FYN: whole kinase-domain context with crystallographic and docked saracatinib.
load outputs/docking_campaign/systems/FYN_saracatinib/prepared/fyn_chain_A_protein.pdb, fyn
load outputs/docking_campaign/systems/FYN_saracatinib/prepared/saracatinib_H8H_A601_crystal.sdf, crystal
load outputs/docking_campaign/figures/presentation/assets/fyn_saracatinib_pose1.sdf, docked
remove solvent
hide everything
show cartoon, fyn
color 0x8EA6B8, fyn
show sticks, crystal or docked
set stick_radius, 0.25, crystal or docked
color 0x00A896, crystal and elem C
color 0xF4A261, docked and elem C
util.cnc crystal
util.cnc docked
select fyn_pocket, byres (fyn within 4.0 of crystal)
show sticks, fyn_pocket and sidechain
set stick_radius, 0.13, fyn_pocket and sidechain
color 0xD5E0E7, fyn_pocket and sidechain and elem C
util.cnc fyn_pocket
orient fyn
turn x, 10
turn y, -25
zoom fyn, 2.5
ray 3000, 2400
png outputs/docking_campaign/figures/preliminary/3D/FYN_overview_crystal_vs_docked.png, dpi=300

# FYN: unobstructed binding-pocket overlay.
hide everything
show cartoon, fyn
color 0xAFC2CE, fyn
set cartoon_transparency, 0.20, fyn
show sticks, fyn_pocket and sidechain
set stick_radius, 0.15, fyn_pocket and sidechain
color 0xD5E0E7, fyn_pocket and sidechain and elem C
util.cnc fyn_pocket
show sticks, crystal or docked
set stick_radius, 0.27, crystal or docked
color 0x00A896, crystal and elem C
color 0xF4A261, docked and elem C
util.cnc crystal
util.cnc docked
zoom (crystal or docked or fyn_pocket), 4.0
turn x, 7
turn y, 18
ray 3000, 2400
png outputs/docking_campaign/figures/preliminary/3D/FYN_pocket_pose_overlay.png, dpi=300

# FYN: translucent local pocket surface.
hide everything
show cartoon, fyn
color 0xB8C9D4, fyn
set cartoon_transparency, 0.35, fyn
select fyn_shell, byres (fyn within 7.0 of crystal)
show surface, fyn_shell
set surface_color, 0xDCE5EA, fyn_shell
set transparency, 0.62, fyn_shell
show sticks, fyn_pocket and sidechain
set stick_radius, 0.13, fyn_pocket and sidechain
color 0xD5E0E7, fyn_pocket and sidechain and elem C
util.cnc fyn_pocket
show sticks, crystal or docked
set stick_radius, 0.27, crystal or docked
color 0x00A896, crystal and elem C
color 0xF4A261, docked and elem C
util.cnc crystal
util.cnc docked
zoom (crystal or docked or fyn_pocket), 4.5
ray 3000, 2400
png outputs/docking_campaign/figures/preliminary/3D/FYN_pocket_surface_overlay.png, dpi=300

# KIT: canonical kinase-domain context with docked masitinib.
delete all
load outputs/docking_campaign/systems/KIT_masitinib/prepared/kit_1T46_chain_A_protein.pdb, kit
load outputs/docking_campaign/figures/presentation/assets/kit_masitinib_pose1.sdf, masitinib
hide everything
show cartoon, kit
color 0x8EA6B8, kit
show sticks, masitinib
set stick_radius, 0.25, masitinib
color 0x7146C7, masitinib and elem C
util.cnc masitinib
select kit_pocket, byres (kit within 4.2 of masitinib)
show sticks, kit_pocket and sidechain
set stick_radius, 0.13, kit_pocket and sidechain
color 0xD5E0E7, kit_pocket and sidechain and elem C
util.cnc kit_pocket
orient kit
turn x, 12
turn y, -22
zoom kit, 2.5
ray 3000, 2400
png outputs/docking_campaign/figures/preliminary/3D/KIT_overview_docked_masitinib.png, dpi=300

# KIT: tight, unobstructed canonical-pocket view.
hide everything
show cartoon, kit
color 0xAFC2CE, kit
set cartoon_transparency, 0.22, kit
show sticks, kit_pocket and sidechain
set stick_radius, 0.15, kit_pocket and sidechain
color 0xD5E0E7, kit_pocket and sidechain and elem C
util.cnc kit_pocket
show sticks, masitinib
set stick_radius, 0.27, masitinib
color 0x7146C7, masitinib and elem C
util.cnc masitinib
zoom (masitinib or kit_pocket), 4.2
turn x, 8
turn y, 18
ray 3000, 2400
png outputs/docking_campaign/figures/preliminary/3D/KIT_pocket_docked_masitinib.png, dpi=300
quit
