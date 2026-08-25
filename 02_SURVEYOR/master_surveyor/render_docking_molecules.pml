# Presentation-grade ray-traced molecular assets.
# Run from the Viscacha_pipeline repository root.

cd /home/welcome3/Viscacha_pipeline
reinitialize
set ray_opaque_background, off
set antialias, 2
set ray_trace_mode, 1
set ray_shadow, 0
set specular, 0.18
set shininess, 18
set depth_cue, 0
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
set stick_radius, 0.18
set sphere_scale, 0.24
set orthoscopic, 1
bg_color white

# FYN overview
load outputs/docking_campaign/FYN_saracatinib/prepared/fyn_chain_A_protein.pdb, fyn
load outputs/docking_campaign/FYN_saracatinib/prepared/saracatinib_H8H_A601_crystal.sdf, crystal
load outputs/docking_campaign/presentation_v2/assets/fyn_saracatinib_pose1.sdf, docked
frame 1
remove solvent
hide everything
show cartoon, fyn
color 0x8EA6B8, fyn
show sticks, crystal
show sticks, docked
color 0x00A896, crystal and elem C
color 0xF4A261, docked and elem C
util.cnc crystal
util.cnc docked
select fyn_pocket, byres (fyn within 4.0 of crystal)
show sticks, fyn_pocket and sidechain
color 0xD7E3EA, fyn_pocket and sidechain and elem C
util.cnc fyn_pocket
orient fyn
turn x, 10
turn y, -25
zoom fyn, 3
ray 2200, 1700
png outputs/docking_campaign/presentation_v2/assets/fyn_overview_ray.png, dpi=300

# FYN pocket close-up with local surface
hide everything
show cartoon, fyn
color 0xB8C9D4, fyn
show surface, (fyn within 7 of crystal)
set surface_color, 0xDCE5EA, fyn
set transparency, 0.58, (fyn within 7 of crystal)
show sticks, fyn_pocket and sidechain
color 0xD7E3EA, fyn_pocket and sidechain and elem C
util.cnc fyn_pocket
show sticks, crystal
show sticks, docked
color 0x00A896, crystal and elem C
color 0xF4A261, docked and elem C
util.cnc crystal
util.cnc docked
zoom (crystal or docked), 7
turn x, 7
turn y, 18
ray 2200, 1700
png outputs/docking_campaign/presentation_v2/assets/fyn_pocket_ray.png, dpi=300

# KIT canonical docking overview
delete all
load outputs/docking_campaign/KIT_masitinib/prepared/kit_1T46_chain_A_protein.pdb, kit
load outputs/docking_campaign/presentation_v2/assets/kit_masitinib_pose1.sdf, masitinib
frame 1
hide everything
show cartoon, kit
color 0x8EA6B8, kit
select kit_pocket, byres (kit within 4.2 of masitinib)
show sticks, kit_pocket and sidechain
color 0xD7E3EA, kit_pocket and sidechain and elem C
util.cnc kit_pocket
show sticks, masitinib
color 0x7B2CBF, masitinib and elem C
util.cnc masitinib
orient kit
turn x, 12
turn y, -22
zoom kit, 3
ray 2200, 1700
png outputs/docking_campaign/presentation_v2/assets/kit_overview_ray.png, dpi=300

# KIT local pocket surface
hide everything
show cartoon, kit
color 0xB8C9D4, kit
show surface, (kit within 7 of masitinib)
set surface_color, 0xDCE5EA, kit
set transparency, 0.58, (kit within 7 of masitinib)
show sticks, kit_pocket and sidechain
show sticks, masitinib
color 0x7B2CBF, masitinib and elem C
util.cnc masitinib
zoom masitinib, 7
turn x, 8
turn y, 18
ray 2200, 1700
png outputs/docking_campaign/presentation_v2/assets/kit_pocket_ray.png, dpi=300
quit
