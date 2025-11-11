#!/usr/bin/env python3

#vivian li oct 30, 2025
#filter through the list of all type II CRISPR systems, if they colocalize with type III, then th

import csv, sys, re

with open ('casII_defense_systems.tsv', 'r') as casII, open ('casIIIA_only_staph_defense_systems.tsv', 'r') as casIII, open('colocalizationlist.txt', 'w') as outfile:
    casII=csv.reader(casII, delimiter= '\t')
    print(str(casII))
    casIII=csv.reader(casIII, delimiter = '\t')
    print(casIII)
    accession= r"(\w*\.\d)(.*)"
    found=[]
    for system in casII:
        geneid=re.search(accession, system)
        for system2 in casIII:
            geneid2=re.search(accession, system)
    if geneid == geneid2:
        found+=geneid2

    for gene in found:
        outfile.write(f'{found}+\n')

