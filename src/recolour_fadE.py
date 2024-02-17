#!/bin/env python3
# author: ph-u
# script: recolour_fadE.py
# desc: recolour fadE peptide structure
# in: "run [full_path]/src/recolour_fadE.py" -> "loadBfacts [protein_file_basename]"
# out: NA
# arg: 0
# date: 20240212

# https://pymol.org/pymol-command-ref.html

# https://stackoverflow.com/questions/10324674/parsing-a-pdb-file-in-python
# https://pymol.org/dokuwiki/doku.php?id=start
# https://pymolwiki.org/index.php/Practical_Pymol_for_Beginners
# https://pymolwiki.org/index.php/Movie_pdf
# https://doi.org/10.4310%2FSII.2015.v8.n4.a4

# https://pymolwiki.org/index.php/Simple_Scripting

#reinitialize

##### PyMol command sequence #####
# load [full_path]/data/D_1292130036_model-annotate_P1.pdb
# load [full_path]/data/D_1292130037_model-annotate_P1-coot1_waterDELETED.pdb
# run [full_path]/src/recolour_fadE.py
# loadBfacts D_1292130036_model-annotate_P1
# run [full_path]/src/recolour_fadE.py
# loadBfacts D_1292130037_model-annotate_P1-coot1_waterDELETED
# super D_1292130036_model-annotate_P1, D_1292130037_model-annotate_P1-coot1_waterDELETED
# select br. all within 2 of sel
##################################

##### Get colour reference file #####
import pandas
d00 = pandas.read_csv("[full_path]/res/02_PAO1_107_PA0506_fadE1_Cysticfibrosis.csv")
d01 = pandas.read_csv("[full_path]/res/02_PAO1_107_PA0508_fadE2_Cysticfibrosis.csv")
d10 = pandas.read_csv("[full_path]/res/02_PAO1_107_PA0506_fadE1_Others.csv")
d11 = pandas.read_csv("[full_path]/res/02_PAO1_107_PA0508_fadE2_Others.csv")
d20 = pandas.read_csv("[full_path]/res/02_PAO1_107_PA0506_fadE1_All.csv")
d21 = pandas.read_csv("[full_path]/res/02_PAO1_107_PA0508_fadE2_All.csv")

##### Reset residues colour #####
# https://www.blopig.com/blog/2020/09/pymol-colouring-proteins-by-property/
from pymol import cmd, stored, math

# https://pymolwiki.org/index.php/Load_new_B-factors
sPec, cOlumn, mAx = "rainbow", ["mEdian","dNdS"], 2.5 # "mEdian" "prot" / "dNdS" "variation"
def loadBfacts (mol, startaa=1, source=d21, visual="Y"):
    """usage: loadBfacts mol, [startaa, [source, [visual]]]"""
    obj = cmd.get_object_list(mol)[0]
    cmd.alter(mol,"b=-1.0")
    counter = int(startaa)
    for x in range(len(source)):
        cmd.alter("%s and resi %s and n. CA"%(mol,counter), "b=%s"%source[cOlumn[0]][x])
        counter += 1
    if visual=="Y":
        cmd.show_as("cartoon", mol)
# angles cgo ellipsoids *licorice* nonbonded *sticks* callback dashes *everything* *lines*
# **ribbon** *surface* **cartoon** dihedrals extent *mesh* slice volume cell *dots*
# labels nb_spheres spheres *wire*
        cmd.cartoon("automatic", mol) # arrow cylinder dumbbell oval rectangle tube automatic dash loop putty skip
        cmd.set("cartoon_putty_scale_min", 0., obj) # min(source[cOlumn])
        cmd.set("cartoon_putty_scale_max", mAx, obj) # max(source[cOlumn]) # 2.5 (dNdS), 1. (prot)
        cmd.set("cartoon_putty_transform", 0, obj)
        cmd.set("cartoon_putty_radius", .1, obj) # peptide string width
# https://pymolwiki.org/index.php/Spectrum
        cmd.spectrum("b", "%s" %sPec, "%s and n. CA " %mol) # apply colour spectrum
        cmd.ramp_new(cOlumn[1], obj, [0., mAx], sPec) # [min(source[cOlumn]), max(source[cOlumn])]        cmd.recolor()
#        cmd.set("surface_color", "dNdS", mol ) # "mesh_color" # https://pymolwiki.org/index.php/Expand_To_Surface

cmd.extend("loadBfacts", loadBfacts);
