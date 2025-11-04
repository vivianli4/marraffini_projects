#!/usr/bin/env python3

#vivian li oct 30, 2025
#filter through the list of all type II CRISPR systems, if they colocalize with type III, then th

import csv, sys, re

with open ('casII_defense_systems.tsv', 'r') as casII, open ('casIIIA_only_staph_defense_systems.tsv', 'r') as casIII, open('colocalizationlist.txt', 'w') as outfile:
    read_casII=csv.DictReader(casII, delimiter= '\t')
    read_casIII=csv.DictReader(casIII, delimiter = '\t')
    accession= r"(\W+\.\W)"
    found=[]
    for sysid in read_casII['sys_id']:
        for sysid2 in read_casIII['sys_id']:
            if re.search(accession, sysid) == re.search(accession,sysid2):
                found+=sysid2.group(1)

    for gene in found:
        outfile.write(f'found+\n')

