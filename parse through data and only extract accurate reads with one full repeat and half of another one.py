from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import nt_search
import operator
import mmap


def spacer_finder(deep_seq_file, name_output_file, genome_for, genome_revcomp, genome_for2, genome_revcomp2):

    handle1 = open(deep_seq_file, "r")
    recs1 = SeqIO.parse(handle1, 'fastq')

    nophage_seq = open(name_output_file, 'w')

    # open the genome txt file needed to match the spacer at the end of this code
    forward = open(genome_for)
    reverse = open(genome_revcomp)

    forward2 = open(genome_for2)
    reverse2 = open(genome_revcomp2)

    spacers_nophage = {}
    reads = 0
    # extract spacers from determined indexes file
    for r in recs1:
        reads = reads + 1
        if reads % 100000 == 0:
            print (reads)
        repeat1 = 'GTTTTAGTACTCTGTAATTTTAGGTATGAGGTAGAC'
        repeat2 = 'GTTTTAGTACTCTGTAAT'

        sequence = r.seq

        rep1 = sequence.find(repeat1)

        if (rep1 != -1):

            sequence = sequence[rep1 + 36:]
            rep2 = sequence.find(repeat2)

            if rep2 != -1:
                spacer = sequence[rep2 - 30:rep2]

                if spacer in spacers_nophage:
                    spacers_nophage[spacer] += 1

                else:
                    spacers_nophage[spacer] = 1

        else:
            rev_comp_seq = sequence.reverse_complement()
            rep1 = rev_comp_seq.find(repeat1)

            if (rep1 != -1):
                rev_comp_seq = rev_comp_seq[rep1 + 36:]
                rep2 = rev_comp_seq.find(repeat2)

                if rep2 != -1:
                    spacer = rev_comp_seq[rep2 - 30:rep2]

                    if spacer in spacers_nophage:
                        spacers_nophage[spacer] += 1

                    else:
                        spacers_nophage[spacer] = 1



    spacers_nophage = sorted(spacers_nophage.items(), key=operator.itemgetter(1), reverse=True)

    # creates a file with spacer, reads, location
    spacer_count = 0
    for spacer in spacers_nophage:

        f = mmap.mmap(forward.fileno(), 0, access=mmap.ACCESS_READ)
        r = mmap.mmap(reverse.fileno(), 0, access=mmap.ACCESS_READ)

        f2 = mmap.mmap(forward2.fileno(), 0, access=mmap.ACCESS_READ)
        r2 = mmap.mmap(reverse2.fileno(), 0, access=mmap.ACCESS_READ)

        # search in the forward orientation
        location = f.find(str(spacer[0]))
        if (location != -1):
            location = location + 31
            pamsequence = f[location-1 :location + 5]
            nophage_seq.write(str(spacer[0]) + "\t" + str(spacer[1])+ "\t" + str(location) + "\t" + str(pamsequence) + "\t" + "top_plasmid" + "\t" + "\n")


        else:
            # search in reverse orientation
            location = r.find(str(spacer[0]))
            if (location != -1):
                locationrev = (len(r) - location - 30)
                pamsequencerev = r[location+30 : location+36]

                nophage_seq.write(
                    str(spacer[0]) + "\t" + str(spacer[1]) + "\t" + str(locationrev) + "\t" + str(pamsequencerev) + "\t" + "bottom_plasmid" + "\t" + "\n")
            else:
                location = f2.find(str(spacer[0]))
            if (location != -1):
                location = location + 31
                pamsequence2 =f2[location-1: location+5]
                nophage_seq.write(
                    str(spacer[0]) + "\t" + str(spacer[1]) + "\t" + str(location) + "\t" +str(pamsequence2) + "\t"+ "top_plasmid2" + "\t" + "\n")

            else:
                location = r2.find(str(spacer[0]))
                if (location != -1):
                    locationrev = (len(r2) - location - 30)
                    pamsequencerev = r[location + 30: location + 36]
                    nophage_seq.write(str(spacer[0]) + "\t" + str(spacer[1]) + "\t" +str(locationrev) + "\t" + str(pamsequencerev)+ "\t" + "bottom_plasmid2" + "\t" + "\n")

                else:
                         nophage_seq.write(str(spacer[0]) + "\t" + str(spacer[1]) + "\t" + '0' + "\t" + "not found" + "\t" + "\n")

    nophage_seq.close()

spacer_finder("VL002_S5_L001_R1_001.fastq", "VL002_plasmidtarget_results.txt", "VL002topstrand.txt", "VL002bottomstrand.txt", "JTR399topstrand.txt", "JTR399bottomstrand.txt")

