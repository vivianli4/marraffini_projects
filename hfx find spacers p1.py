from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import nt_search
import operator
import mmap

def spacer_finder(deep_seq_file, name_output_file, genome_for, genome_revcomp, genome_for2, genome_revcomp2,
                    genome_for3, genome_revcomp3, genome_for4, genome_revcomp4, genome_for5, genome_revcomp5):

    handle1 = open(deep_seq_file, "r")
    recs1 = SeqIO.parse(handle1, 'fastq')

    nophage_seq = open(name_output_file, 'w')

    # open the genome txt file needed to match the spacer at the end of this code
    forward = open(genome_for)
    reverse = open(genome_revcomp)

    forward2 = open(genome_for2)
    reverse2 = open(genome_revcomp2)

    forward3 = open(genome_for3)
    reverse3 = open(genome_revcomp3)

    forward4 = open(genome_for4)
    reverse4 = open(genome_revcomp4)

    forward5 = open(genome_for5)
    reverse5 = open(genome_revcomp5)

    spacers_nophage = {}
    reads = 0
    # extract spacers from determined indexes file
    for r in recs1:
        reads = reads + 1
        if reads % 100000 == 0:
            print (reads)
        repeat1 = 'GTTTCAGACGAACCCTTGTGGGATTGAAGC'
        repeat2 = 'GTTTCAGACGAACCCTT'

        sequence = r.seq

        rep1 = sequence.find(repeat1)

        if (rep1 != -1):

            sequence = sequence[rep1 + 30:]
            rep2 = sequence.find(repeat2)

            if rep2 != -1:
                spacer = sequence[:rep2]

                if spacer in spacers_nophage:
                    spacers_nophage[spacer] += 1

                else:
                    spacers_nophage[spacer] = 1

        else:
            rev_comp_seq = sequence.reverse_complement()
            rep1 = rev_comp_seq.find(repeat1)

            if (rep1 != -1):
                rev_comp_seq = rev_comp_seq[rep1 + 30:]
                rep2 = rev_comp_seq.find(repeat2)

                if rep2 != -1:
                    spacer = rev_comp_seq[rep2 - 36:rep2]

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
        f3 = mmap.mmap(forward3.fileno(), 0, access=mmap.ACCESS_READ)
        r3 = mmap.mmap(reverse3.fileno(), 0, access=mmap.ACCESS_READ)
        f4 = mmap.mmap(forward4.fileno(), 0, access=mmap.ACCESS_READ)
        r4 = mmap.mmap(reverse4.fileno(), 0, access=mmap.ACCESS_READ)
        f5 = mmap.mmap(forward5.fileno(), 0, access=mmap.ACCESS_READ)
        r5 = mmap.mmap(reverse5.fileno(), 0, access=mmap.ACCESS_READ)


        # search in the forward orientation
        location = f.find(str(spacer[0]))
        if (location != -1):
            location = location + 36
            nophage_seq.write(str(spacer[0]) + "\t" + str(spacer[1]) + "\t" + str(location) + "\t" + "top_chromosome" + "\t" + "\n")

        else:
            # search in reverse orientation
            location = r.find(str(spacer[0]))
            if (location != -1):
                pam = (len(r) - location - 35)
                nophage_seq.write(
                    str(spacer[0]) + "\t" + str(spacer[1]) + "\t" + str(pam) + "\t" + "bottom_chromosome" + "\t" + "\n")


            else:
                location = f2.find(str(spacer[0]))
                if (location != -1):
                    location = location + 36
                    nophage_seq.write(str(spacer[0]) + "\t" + str(spacer[1]) + "\t" + str(location) + "\t" + "top_pHV1" + "\t" + "\n")

                else:
                    location = r2.find(str(spacer[0]))
                    if(location != -1):
                        pam = (len(r2) - location - 35)
                        nophage_seq.write(str(spacer[0]) + "\t" + str(spacer[1]) + "\t" + str(pam) + "\t" + "bottom_pHV1" + "\t" + "\n")


                    else:
                        location = f3.find(str(spacer[0]))
                        if (location != -1):
                            location = location + 36
                            nophage_seq.write(str(spacer[0]) + "\t" + str(spacer[1]) + "\t" + str(
                                location) + "\t" + "top_pHV2" + "\t" + "\n")

                        else:
                            location = r3.find(str(spacer[0]))
                            if (location != -1):
                                pam = (len(r2) - location - 35)
                                nophage_seq.write(str(spacer[0]) + "\t" + str(spacer[1]) + "\t" + str(
                                    pam) + "\t" + "bottom_pHV2" + "\t" + "\n")


                            else:
                                location = f4.find(str(spacer[0]))
                                if (location != -1):
                                    location = location + 36
                                    nophage_seq.write(str(spacer[0]) + "\t" + str(spacer[1]) + "\t" + str(
                                        location) + "\t" + "top_pHV3" + "\t" + "\n")

                                else:
                                    location = r4.find(str(spacer[0]))
                                    if (location != -1):
                                        pam = (len(r2) - location - 35)
                                        nophage_seq.write(str(spacer[0]) + "\t" + str(spacer[1]) + "\t" + str(
                                            pam) + "\t" + "bottom_pHV3" + "\t" + "\n")


                                    else:
                                        location = f5.find(str(spacer[0]))
                                        if (location != -1):
                                            location = location + 36
                                            nophage_seq.write(str(spacer[0]) + "\t" + str(spacer[1]) + "\t" + str(
                                                location) + "\t" + "top_pHV4" + "\t" + "\n")

                                        else:
                                            location = r5.find(str(spacer[0]))
                                            if (location != -1):
                                                pam = (len(r2) - location - 35)
                                                nophage_seq.write(str(spacer[0]) + "\t" + str(spacer[1]) + "\t" + str(
                                                    pam) + "\t" + "bottom_pHV4" + "\t" + "\n")

                                            else:
                                                nophage_seq.write(str(spacer[0]) + "\t" + str(spacer[1]) + "\t" + '0' + "\t" + "not found" + "\t" + "\n")

    nophage_seq.close()

spacer_finder("HFPV1-p1_S6_L001_R1_001.fastq", "HFPV1_1_p1_acquired_test.txt", "gb_Hfvol_DS2_chromosome_for.txt", "gb_Hfvol_DS2_chromosome_rev_comp.txt", "pHV1_for.txt", "pHV1_rev_comp.txt",
              "pHV2_for.txt", "pHV2_rev_comp.txt", "pHV3_for.txt", "pHV3_rev_comp.txt", "pHV4_for.txt", "pHV4_rev_comp.txt")

