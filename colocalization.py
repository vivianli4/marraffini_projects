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
    for system in casII:
        gene_id=re.search(accession, system[0]).group(1)
        if gene_id== 'GCA_019998745.2':
            print('found it')
        for system2 in casIII:
            gene_id2=re.search(accession, system2[0]).group(1)
            if gene_id2=='GCA_019998745.2':
                print('found it here too')

        if str(gene_id) == str(gene_id2):
            found.append(gene_id2)
            print('compared the two')

    for gene in found:
        outfile.write(f'{found}+\n')

#for some reason, it's searching for the system in both files, but it's not comparing the two