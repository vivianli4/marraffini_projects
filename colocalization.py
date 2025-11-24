#!/usr/bin/env python3

#vivian li oct 30, 2025
#filter through the list of all type II CRISPR systems, if they colocalize with type III, then th

import csv, sys, re

with open ('casII_defense_systems.tsv', 'r') as casII, open ('casIIIA_only_staph_defense_systems.tsv', 'r') as casIII, open('colocalizationlist.txt', 'w') as outfile:
    accession= r"(\w*\.\d)(.*)"
    found=[]
    headings=next(casII)
    headings2=next(casIII)
    casII = csv.reader(casII, delimiter='\t')
    casIII = csv.reader(casIII, delimiter='\t')

    # Extract prefixes from casII
    casII_prefixes = set()
    for system in casII:
        match = re.search(accession, system[0])
        if match:
            casII_prefixes.add(match.group(1))

    # Extract prefixes from casIII and find matches
    found = []
    for system2 in casIII:
        match = re.search(accession, system2[0])
        if match:
            gene_id2 = match.group(1)
            if gene_id2 in casII_prefixes:
                found.append(gene_id2)
                print(f'Found matching prefix: {gene_id2}')

    outfile.write(f'{found} \n')