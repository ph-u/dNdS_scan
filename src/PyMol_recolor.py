#!/bin/env python3
# author: ph-u
# script: PyMol_recolor.py
# desc: recolour a peptide structure
# in: "run /[full_path]/PyMol_recolor.py" -> "loadBfacts [protein_file_basename]"
# out: NA
# arg: 0
# date: 20240212, 20240320

##### PyMol command sequence #####
### monomer ###
# 1. load /[full/path]/AF-Q9HZR3-F1-model_v4.cif
# 2. run [full/path]/PyMol_recolor.py
# 3. loadBfacts AF-Q9HZR3-F1-model_v4
## {repeat steps 1-3}
# 4. set_name AF-Q9HZR3-F1-model_v4, PA2934.{ce}

### polymer ###
# 1. load /[full/path]/3kd2.cif
# 2. set_name 3kd2, 3kd2.{ce}
## {repeat steps 1-2}
# 3. super 3kd2.c, 3kd2.e
# 4. run [full/path]/PyMol_recolor.py
# 5. select m, (resi 25-317 and chain {ABCD} and 3kd2.{ce})
# 6. extract m{ABCD}.{ce}, m
# 7. loadBfacts m{ABCD}.{ce}
## {repeat steps 5-7 for 8 times, each 4 times run once step 4 with a modified srcNam number: each polymer has 4 subunits, 2 polymers in total}

### other useful commands ###
# bg_color white
# set grid_mode,1
# delete 3kd2
##################################

##### Get colour reference file #####
import pandas, math
pT0, pT1 = "/[full_path]/[dNdS_prefix]_", "--reCon.csv"
scaleMax, srcNam = [0.,1.2], ["Cystic-fibrosis","Environmental"]

d = pandas.read_csv(pT0 + srcNam[0] + pT1) # !!!!! CHANGE THE BRACKETED NUMBER HERE !!!!!

##### Reset residues colour #####
from pymol import cmd, stored, math

sPec, cOlumn = "rainbow", ["dNdS.mean","dNdS"]
def loadBfacts (mol, startaa=1, source=d, visual="Y"):
    """usage: loadBfacts mol, [startaa, [source, [visual]]]"""
    source.at[startaa-1,cOlumn[0]] = scaleMax[0]
    source.at[len(source)-1,cOlumn[0]] = scaleMax[1]
    obj = cmd.get_object_list(mol)[0]
    cmd.alter(mol,"b=-1.0")
    counter = int(startaa)
    xMin = min(source[cOlumn[0]])
    if math.isinf(min(source[cOlumn[0]])):
        xMax = max(source[cOlumn[0]])+2
    else:
        xMax = max(source[cOlumn[0]])

    for x in range(len(source)):
        if math.isinf(source[cOlumn[0]][x]):
            x0 = max(source[cOlumn[0]])+2
        else:
            x0 = source[cOlumn[0]][x]

        cmd.alter("%s and resi %s and n. CA"%(mol,counter), "b=%s"%x0)
        counter += 1
    if visual=="Y":
        cmd.show_as("cartoon", mol)
        cmd.cartoon("automatic", mol)
        cmd.set("cartoon_putty_scale_min", xMin, obj)
        cmd.set("cartoon_putty_scale_max", xMax, obj)
        cmd.set("cartoon_putty_transform", 0, obj)
        cmd.set("cartoon_putty_radius", .1, obj) # peptide string width
        cmd.spectrum("b", "%s" %sPec, "%s and n. CA " %mol) # apply colour spectrum
        cmd.ramp_new(cOlumn[1], obj, [xMin, xMax], sPec)

cmd.extend("loadBfacts", loadBfacts);
