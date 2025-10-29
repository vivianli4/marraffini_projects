#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "Marina Shi"
# Last Updated: 11-23-2024


# These libraries exist in standard Python, so you don't need to download anything extra. Other imported libraries may
# require downloading stuff -- in this case, I usually use Anaconda virtual environments to keep track of everything and
# prevent compatibility issues.
import os
from csv import DictReader, DictWriter

# Define global constants (should never be changed anywhere in the script) -- use all caps for these variable names
# These should be the folder where your input file resides and the folder where you want to put the output file
IN_DIR = '/Users/vivian/Documents/staph defense system colocalization'
OUT_DIR = '/Users/vivian/Documents/staph defense system colocalization'


# Define a function to filter input file for keyword and write new file containing info for all matches
def filter_write_rows(infile_path, outfile_path, filter_keyword):
    filtered_info_dict = dict()
    filtered_info_count = 0
    # Open the input file as read-only to protect original file
    # Using "with open(...) as ..." is safer practice than doing separate open() and close() statements
    with open(infile_path, 'r') as infile:
        # DictReader detects the first row in the file as fieldnames (in_reader.fieldnames) and then skips the row
        # The delimiter parameter tells DictReader that this is a tab-separated file (if you don't define the delimiter,
        # the default assumption is that the input is a comma-separated file)
        in_reader = DictReader(infile, delimiter='\t')
        for row in in_reader:
            # Check if your keyword exists under the 'sys_id' column in the current row
            if filter_keyword in row['sys_id']:
                filtered_info_count += 1
                filtered_info_dict[filtered_info_count] = dict()
                # Save all info from current row in input file
                for field in in_reader.fieldnames:
                    filtered_info_dict[filtered_info_count][field] = row[field]
    print('# of', filter_keyword, 'systems found:', filtered_info_count)

    # Create the output file and open as write-only
    # Warning!!! If you already have something in the specified output file, it will be deleted and overwritten with
    # whatever it's told to write here. If you want to add on to the end of an existing file, use 'a' instead of 'w'.
    with open(outfile_path, 'w') as outfile:
        out_writer = DictWriter(outfile, fieldnames=in_reader.fieldnames, delimiter='\t')
        out_writer.writeheader()
        for val in filtered_info_dict.values():
            # DictWriter will write row from a dictionary
            out_writer.writerow(val)


# The line below is a safety precaution that I put in all of my Python scripts (I think it's considered the standard
# in field). It means that if you have any other Python scripts that reference some specific functions/other stuff from
# this script, it won't accidentally also run whatever is below here.
if __name__ == "__main__":
    # I decided to name these variables before the function call to make it easier to change what you want in and out.
    # If you want to filter for a different defense system, you can change system_for_filter (I set it to
    # 'CAS_Class1-Subtype-III-A', but if you want type II-A instead, you can change it to 'CAS_Class2-Subtype-II-A'.
    # Make sure to change output_file_name if you don't want to overwrite the III-A file, see the warning above.)
    input_file_name = 'all_staph_defense_systems.tsv'
    output_file_name = 'casII_defense_systems.tsv'
    system_for_filter = 'CAS_Class2-Subtype-II'

    filter_write_rows(os.path.join(IN_DIR, input_file_name), os.path.join(OUT_DIR, output_file_name), system_for_filter)
