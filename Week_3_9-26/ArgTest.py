#ArgTest

import argparse
import numpy
parser=argparse.ArgumentParser()

parser.add_argument('--input_file',action='store') #
parser.add_argument('--output_fasta_file',action='store',default='FASTQ_FASTA_Output.txt') #
parser.add_argument('--output_fasta',action='store_true',default='false') #
parser.add_argument('--output_screen',action='store_true',default=False) #
parser.add_argument('--output_filename',action='store') #
parser.add_argument('--fasta',action='store_true',default=False) 
parser.add_argument('--keep_ambiguous',action='store_true',default=False) #
parser.add_argument('--ltrim',action='store',default=0) #
parser.add_argument('--rtrim',action='store',default=0) #
parser.add_argument('--qual_trim',action='store',default=0) #
parser.add_argument('--base_trim',action='store',default=25)

args=parser.parse_args()

with open(args.output_fasta_file,"w") as fasta_output:
	with open(args.input_file,"r") as fastq:
	    BP_Num=0
		Kept_Reads=0
		Read_Num=sum(1 for line in fastq)/4
		fastq.seek(0)

		row, column=2,1000
		Quality_Matrix=[[0 for x in range(column)] for y in range(row)]
		
		Read_Length=list()
		Read_Quality=list()

		for i in range(0,Read_Num):
			trim=False
			fastq.readline()
			Read=fastq.readline()		#Find Read
			
			fastq.readline()
			Quality=fastq.readline()	#Find Quality

			if args.ltrim>0:
				Read=Read[args.ltrim:]
			if args.rtrim>0:
				Read=Read[:-args.rtrim]

			Read_Quality=list()
			Read_Quality.append(numpy.mean([ord(x)-33 for x in list(Quality)[0:(len(Quality)-1)]]))

			if not args.keep_ambiguous:	#If ambiguous/excluded, skip to next block
				if 'N' in Read:
					trim=True
			if Read_Quality<args.qual_trim:
				trim=True

			if not trim:	#Read Included
				if args.output_fasta:
					fasta_output.write("> %s" % Read)
				Read_Length.append(len(Read)-1)
				BP_Num+=len(Read)-1
				Kept_Reads+=1
				for j in range(0,len(Quality)-1):
					Quality_Matrix[0][j]=Quality_Matrix[0][j]+ord(Quality[j])-33
					Quality_Matrix[1][j]=Quality_Matrix[1][j]+1

Quality_Score=Quality_Matrix[0][0:max(Read_Length)]
Quality_Count=Quality_Matrix[1][0:max(Read_Length)]

print Quality_Score
print Quality_Count
Quality_Base=[x/y for x,y in zip(Quality_Score,Quality_Count)]
print Quality_Base

with open(args.output_filename,"w") as export:
	export.write("Summary of %s\n\n" % args.input_file)
	export.write("Total Reads: %d\n" % Read_Num)
	export.write("Kept Reads: %d\n" % Kept_Reads)
	export.write("Total Base Pairs: %dbp\n" % (sum(Read_Length)))
	export.write("Mean (Min, Max) Read Length: %d (%d,%d)\n" % (numpy.mean(Read_Length),min(Read_Length),max(Read_Length)))
	export.write("Read Lengths:\n") #Distribution)
	for i in range(min(Read_Length)-1,max(Read_Length)):
		export.write("     %d: %d\n" % (i+1,Read_Length.count(i+1)))
	export.write("Mean (Min, Max) Read Quality: %d (%d,%d) \n" % (numpy.mean(Read_Quality),min(Read_Quality),max(Read_Quality)))
	export.write("Mean Quality Per Base:\n")
	for i in range(0,len(Quality_Base)):
		export.write("     %d: %d\n" % (i+1,Quality_Base[i]))
	export.close()

if args.output_screen:
	with open(args.output_filename,"r") as stat_display:
		for line in stat_display:
			print line

