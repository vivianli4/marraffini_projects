#!/usr/bin/env python3

#vivian li oct 30, 2025
#filter through the list of all type II CRISPR systems, if they colocalize with type III, then th

import csv, sys, re

with open ('casII_defense_systems.tsv', 'r') as casII, open ('casIIIA_only_staph_defense_systems.tsv', 'r') as casIII, open('colocalizationlist.txt', 'w') as outfile:
    # read_casII=csv.DictReader(casII, delimiter= '\t')
    # read_casIII=csv.DictReader(casIII, delimiter = '\t')
    accession= r'(\w*\.\d)(.*)'
    found=[]
    for line in casII:
        for line2 in casIII:
            typeII=re.search(accession, line)
            print(typeII)
            typeIII=re.search(accession, line2)
            print(typeIII)
            if typeII == typeIII:
                found+= typeII.group(0)

    for gene in found:
        outfile.write(f'{found}+\n')

