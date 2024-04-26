#!/bin/env python3
# author: ph-u
# script: PA2185_recolour.py
# desc: recolour PA2185 peptide structure
# in: "run /[full_path]/PA2185_recolour.py" -> "loadBfacts [protein_file_basename]"
# out: NA
# arg: 0
# date: 20240212, 20240320

##### PyMol command sequence #####
# load [full/path]/AF-Q9I1T0-F1-model_v4.cif
# run [full/path]/PA2185_recolour.py
# loadBfacts AF-Q9I1T0-F1-model_v4
# bg_color white
# set_name dNdS, dNdS.{CF,env,oth,skin,unk}
##################################

##### Get colour reference file #####
import pandas, math
pT0, pT1 = "/[full_path]/PAO1_107_PA2185_", "--reCon.csv"
scaleMax, srcNam = 1.3, ["Cystic-fibrosis","Environmental","Other-infections","Skin","Unknown"]

d = pandas.read_csv(pT0 + srcNam[0] + pT1)

##### Reset residues colour #####
from pymol import cmd, stored, math

sPec, cOlumn = "rainbow", ["dNdS.mean","dNdS"]
def loadBfacts (mol, startaa=1, source=d, visual="Y"):
    """usage: loadBfacts mol, [startaa, [source, [visual]]]"""
    source.at[len(source)-1,"dNdS.mean"] = scaleMax
    obj = cmd.get_object_list(mol)[0]
    cmd.alter(mol,"b=-1.0")
    counter = int(startaa)
    if math.isinf(min(source[cOlumn[0]])):
        xMin = 0
        xMax = max(source[cOlumn[0]])+2
    else:
        xMin = 0 #min(source[cOlumn[0]])
        xMax = scaleMax #max(source[cOlumn[0]])

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
