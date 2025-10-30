from Bio import SeqIO

import operator
import mmap


def spacer_finder(deep_seq_file, undetermined_deep_seq_file, name_output_file, seq_primer , genome_for, genome_revcomp):

    handle1 = open(deep_seq_file, "r")
    recs1 = SeqIO.parse(handle1, 'fastq')


    nophage_seq = open(name_output_file, 'w')

    # open the genome txt file needed to match the spacer at the end of this code
    forward = open(genome_for)
    reverse = open(genome_revcomp)

    spacers_nophage = {}
    reads = 0
    # extract spacers from determined indexes file
    for r in recs1:
        reads = reads + 1
        if reads % 100000 == 0:
            print reads
        repeat1 = 'GTTTTAGAGCTATGCTGTTTTGAATGGTCCCAAAAC'
        repeat2 = 'GTTTTAGAGCTATGCTGTTTT'

        sequence = r.seq

        # skip reads longer than 146 bp
        if len(sequence) > 146:
            pass

        else:

            # nophage experiment
            # this search for the forward primer PM171
            posF = sequence.find(seq_primer)

            if (posF != -1):
                rep1 = sequence.find(repeat1)

                if (rep1 != -1):
                    sequence = sequence[rep1 + 36:]

                    rep2 = sequence.find(repeat2)
                    if rep2 != -1:
                        spacer = sequence[rep2 - 20:rep2]

                        if spacer in spacers_nophage:
                            spacers_nophage[spacer] += 1
                        else:
                            spacers_nophage[spacer] = 1
            else:
                rev_comp_seq = sequence.reverse_complement()

                # search for the barcoded primer in the reverse complement sequence
                posR = rev_comp_seq.find(seq_primer)
                if (posR != -1):
                    rep1 = rev_comp_seq.find(repeat1)

                    if (rep1 != -1):
                        rev_comp_seq = rev_comp_seq[rep1 + 36:]

                        rep2 = rev_comp_seq.find(repeat2)
                        if rep2 != -1:
                            spacer = rev_comp_seq[rep2 - 20:rep2]

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

        # search in the forward orientation
        location = f.find(str(spacer[0]))
        if (location != -1):
            location = location + 21
            nophage_seq.write(str(spacer[0]) + "\t" + str(spacer[1]) + "\t" + str(location) + "\t" + "top_phage" + "\t" + "\n")


        else:
            # search in reverse orientation
            location = r.find(str(spacer[0]))
            if (location != -1):
                pam = (len(r) - location - 20)
                nophage_seq.write(
                    str(spacer[0]) + "\t" + str(spacer[1]) + "\t" + str(pam) + "\t" + "bottom_phage" + "\t" + "\n")

            else:
                nophage_seq.write(str(spacer[0]) + "\t" + str(spacer[1]) + "\t" + '0' + "\t" + "not found" + "\t" + "\n")

    nophage_seq.close()

spacer_finder("pPM120-pPM118control_S1_L001_R1_001.fastq","pPM120 Undetermined_S0_L001_R1_001.fastq" , "pPM120-pPM118 no phage.txt", "CAAAGTGCGATTACAAAATTTTTTAGAC" , "NM4_right_origin_for.txt", "NM4_right_origin_revcomp.txt")
spacer_finder("pPM120-pPM118NM4g4_S2_L001_R1_001.fastq","pPM120 Undetermined_S0_L001_R1_001.fastq" , "pPM120-pPM118 NM4.txt", "GTTGAGTGCGATTACAAAATTTTTTAGAC" , "NM4_right_origin_for.txt", "NM4_right_origin_revcomp.txt")
spacer_finder("pPM120-pPM118-dBgLII_S3_L001_R1_001.fastq","pPM120 Undetermined_S0_L001_R1_001.fastq", "pPM120-pPM118 NM4-dBglII.txt", "ATTGGAGTGCGATTACAAAATTTTTTAGAC" , "NM4_dBglII_right_origin_for.txt", "NM4_dBglII_right_origin_revcomp.txt")
spacer_finder("120-118-NM4-96_S1_L001_R1_001.fastq","120 Undetermined_S0_L001_R1_001.fastq", "pPM120-pPM118 NM4-96.txt", "CAAAGTGCGATTACAAAATTTTTTAGAC" , "NM4-96_right_origin_for.txt", "NM4-96_right_origin_revcomp.txt")
spacer_finder("pPM120-pPM118-98_S5_L001_R1_001.fastq","pPM120 Undetermined_S0_L001_R1_001.fastq", "pPM120-pPM118 NM4-98.txt", "GACCAGTGCGATTACAAAATTTTTTAGAC" , "NM4-98_right_origin_for.txt", "NM4-98_right_origin_revcomp.txt")
spacer_finder("120-118-dBgLII-96_S3_L001_R1_001.fastq", "120 Undetermined_S0_L001_R1_001.fastq","pPM120-pPM118 NM4-dBglII-96.txt", "ATTGGAGTGCGATTACAAAATTTTTTAGAC" , "NM4_dBglII-96_right_origin_for.txt", "NM4_dBglII-96_right_origin_revcomp.txt")
spacer_finder("pPM120-pPM118-dBg-98_S6_L001_R1_001.fastq", "pPM120 Undetermined_S0_L001_R1_001.fastq","pPM120-pPM118 NM4-dBglII-98.txt", "TTCCGAGTGCGATTACAAAATTTTTTAGAC" , "NM4_dBglII-98_right_origin_for.txt", "NM4_dBglII-98_right_origin_revcomp.txt")