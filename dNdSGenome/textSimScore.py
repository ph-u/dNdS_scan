#!/bin/env python3
# author: ph-u
# script: textSimScore.py
# desc: calculate similarity score for two text strings
# in: python3 textSimScore.py [text 1] [text 2]
# out: stdout message
# arg: 2
# date: 20240105

import sys
from difflib import SequenceMatcher

##### f: text similarity #####
def txSim(a0, a1):
    if (a0=="") | (a1==""):
        return "NA"
    else:
        return int(round(SequenceMatcher(None, a0, a1).ratio()*100))

##### Calculate text similarity score #####
print(txSim(sys.argv[1], sys.argv[2]))
