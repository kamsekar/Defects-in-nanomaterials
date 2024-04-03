# A script that generates a crystallographic file (.cif) for a modulated cell, consisting of several
# PtCu3 Pm-3m unit cells, separated by periodic anti-phase boundaries.
# ana.rebeka.kamsek@ki.si, 2023

# specify how many Pm-3m unit cells make up the domains on each side of periodic anti-phase boundaries
first = 6
second = 9

# specify path to save the file
file_name = "modulated_Cu3Pt_" + str(first) + "_" + str(second) + "_P4mm.cif"

# unit cell parameter of Pm-3m unit cells
a = 3.7

cells = first + second
c = a * cells
phase_name = "Cu3_Pt_" + str(cells) + "-" + str(first) + "-" + str(second) + "_Pmmm"

# consider the constituent atoms according to their relative position
Pt_first_border = []
Cu_first_middle = []
Cu_first_outer = []

Pt_second_middle = []
Cu_second_border = []
Cu_second_outer = []

# one set of constituent atoms is made up of repeating Pm-3m cells
for i in range(first):
    Pt_first_border.append("Pt" + str(i + 1) + "\tPt\t0.0\t0.0\t" + str(i / cells))
    Cu_first_middle.append("Cu" + str(i + 1) + "\tCu\t0.5\t0.5\t" + str(i / cells))
    Cu_first_outer.append("Cu" + str(first + i + 1) + "\tCu\t0.0\t0.5\t" + str((i + 0.5) / cells))

# the second set is made up of displaced repeating Pm-3m cells
for i in range(second):
    Pt_second_middle.append("Pt" + str(first + i + 1) + "\tPt\t0.5\t0.5\t" + str((first + i) / cells))
    Cu_second_border.append("Cu" + str(2 * first + i + 1) + "\tCu\t0.0\t0.0\t" + str((first + i) / cells))
    Cu_second_outer.append("Cu" + str(2 * first + second + i + 1) + "\tCu\t0.0\t0.5\t" + str((first + i + 0.5) / cells))

# create a text file and write the crystallographic information about the created phase
with open(file_name, "w") as f:
    f.write("\ndata_" + phase_name + "\n\n")
    f.write("_pd_phase_name " + "\"" + phase_name + "\"" + "\n")
    f.write("_cell_length_a\t" + str(a) + "\n")
    f.write("_cell_length_b\t" + str(a) + "\n")
    f.write("_cell_length_c\t" + str(c) + "\n")
    f.write("_cell_angle_alpha\t90" + "\n")
    f.write("_cell_angle_beta\t90" + "\n")
    f.write("_cell_angle_gamma\t90" + "\n")
    f.write("_symmetry_cell_setting\ttetragonal" + "\n")
    f.write("_symmetry_space_group_name_H-M\t\"P 4 m m\"\n\n")
    f.write("loop_\n\t_atom_site_label\n\t_atom_site_type_symbol\n"
            "\t_atom_site_fract_x\n\t_atom_site_fract_y\n\t_atom_site_fract_z\n")

    f.write("# Pt, first set of cells, border (Wyckoff 1a 4mm 0 0 z)\n")
    for i in range(first):
        f.write(Pt_first_border[i] + "\n")
    f.write("# Pt, second set of cells, middle (Wyckoff 1b 4mm 1/2 1/2 z)\n")
    for i in range(second):
        f.write(Pt_second_middle[i] + "\n")
    f.write("# Cu, first set of cells, middle (Wyckoff 1b 4mm 1/2 1/2 z)\n")
    for i in range(first):
        f.write(Cu_first_middle[i] + "\n")
    f.write("# Cu, first set of cells, outer (Wyckoff 2c mm 0 1/2 z)\n")
    for i in range(first):
        f.write(Cu_first_outer[i] + "\n")
    f.write("# Cu, second set of cells, border (Wyckoff 1a 4mm 0 0 z)\n")
    for i in range(second):
        f.write(Cu_second_border[i] + "\n")
    f.write("# Cu, second set of cells, outer (Wyckoff 2c mm 0 1/2 z)\n")
    for i in range(second):
        f.write(Cu_second_outer[i] + "\n")
f.close()
