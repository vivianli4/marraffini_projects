import csv


# parse through a cvs file of exce
# creates a dictionary to associate location and reads
dict_location_reads = {}
with open('.csv', 'rU') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        # if statement to add number of reads for the location if it's already in the dictionary
        if int(row['location']) in dict_location_reads:
            print("yaaaas")
            dict_location_reads[int(row['location'])] += float(row['reads'])
        else:
            dict_location_reads[int(row['location'])] = float(row['reads'])


# create an empty dictionary with each bin 0 to 40 equal to 0
bin_reads = {}
bin_key = 0
bin_reads[bin_key] = 0
while bin_key != 439:
    bin_key = int(bin_key + 1)
    bin_reads[bin_key] = 0


print
# this go through the dictionary dict_location_reads and assign the number of reads for each bin location into the empyt bin dictionary just created
for key, value in dict_location_reads.iteritems():
    bin_location = int(key / 100)
    bin_reads[bin_location] += value


# create a txt file with the first column being the bin location 0=0to502bp, 1=503to1005bp... and the second column the number of spacer reads for bin location
txtfile = open('399-5_results100bin.txt', 'w')


for key, value in bin_reads.iteritems():
    txtfile.write(str(key) + "\t" + str(value) + "\n")


txtfile.close()


