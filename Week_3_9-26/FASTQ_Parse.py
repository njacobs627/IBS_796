import sys
import numpy
import argparse

source_filename=sys.argv[1]
output_filename=sys.argv[2]

#Requirements
	#Read Number							X
	#BP Number								X
	#Min, Max, Mean/Median Read Length		X
	#Read Length Distribution				X
	#Min, Max, Mean/Median Quality/Read 	X
	#Mean Quality/Base						X

with open(source_filename,"r") as fastq:
	Read_Num=sum(1 for line in fastq)/4
	fastq.seek(0)
	BP_Num=0
	
	Read_Length=list()
	Read_Quality=list()


	for i in range(0,Read_Num):
		fastq.readline() #Seq Header
		Read=fastq.readline() #Read
		fastq.readline() #Phred Symbol
		fastq.readline()
		Read_Length.append(len(Read)-1)
		
	row, column=2,max(Read_Length)
	Quality_Matrix=[[0 for x in range(column)] for y in range(row)]

	fastq.seek(0)
	for i in range(0,Read_Num):
		fastq.readline() #Seq Header
		fastq.readline() #Read
		fastq.readline() #Phred Symbol
		Quality=fastq.readline() #Quality Score
		Read_Quality.append(numpy.mean([ord(x)-33 for x in list(Quality)[0:(len(Quality)-1)]]))

		for j in range(0,len(Quality)-1):
			Quality_Matrix[0][j]=Quality_Matrix[0][j]+ord(Quality[j])-33
			Quality_Matrix[1][j]=Quality_Matrix[1][j]+1
	
	Quality_Score=Quality_Matrix[0]
	Quality_Count=Quality_Matrix[1]
	Quality_Base=[x/y for x,y in zip(Quality_Score,Quality_Count)]
	fastq.close()

with open(output_filename,"w") as export:
	export.write("Summary of %s\n\n" % source_filename)
	export.write("Total Reads: %d\n" % Read_Num)
	export.write("Total Base Pairs: %dbp\n" % (sum(Read_Length)))
	export.write("Mean (Min, Max) Read Length: %d (%d,%d)\n" % (numpy.mean(Read_Length),min(Read_Length),max(Read_Length)))
	export.write("Read Lengths:\n") #Distribution)
	for i in range(min(Read_Length)-1,max(Read_Length)):
		export.write("     %d: %d\n" % (i+1,Read_Length.count(i+1)))
	export.write("Mean (Min, Max) Read Quality: %d (%d,%d) \n" % (numpy.mean(Read_Quality),min(Read_Quality),max(Read_Quality)))
	export.write("Mean Quality Per Base:\n")
	for i in range(0,len(Quality_Base)):
		export.write("     %d: %d\n" % (i+1,Quality_Base[i]))