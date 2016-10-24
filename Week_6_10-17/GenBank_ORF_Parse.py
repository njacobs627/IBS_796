#!/usr/bin/env python

#GenBank Parser.py
def Complement(Sequence,Molecule):
	if Molecule=="DNA":
		return Sequence.translate(string.maketrans("ATCG","TAGC"))
	if Molecule=="RNA":
		return Sequence.translate(string.maketrans("AUCG","UAGC"))	

def ReverseComplement(Sequence,Molecule):
	return Complement(Sequence,Molecule)[::-1]

def ORF_Find(Input,Output):
	with open(Input,"r") as GenBank:
		with open(Output,"w") as ORF_Output:	
			
			#Extract DNA Sequence from FASTA file
			DNA=""
			DNA_Line=GenBank.readline() #Ignore Header Line
			while len(DNA_Line)>0:
				DNA_Line=GenBank.readline()
				DNA=DNA+DNA_Line.replace("\n","")

			pattern=r'(?=(ATG(?:...)*?)(?=TAG|TGA|TAA))'	#RE for START, codons, STOP
			
			#Find ORFs on + Strand
			b=list(re.finditer(pattern,DNA)) 	#Find all potential ORFs (START/STOP, bp/3 = int)

			StartPos=[a.start(0) for a in b]
			EndPos=[a.end(1) for a in b]
			
			i=0
			ORF_Count=1
			ORF_Output.write(">ORF_%05d Strand: (-) Start: %d  End: %d\n" % (ORF_Count,StartPos[i],EndPos[i]+3))
			ORF_Output.write("%s\n" % DNA[StartPos[i]:EndPos[i]+3])
			ORFStart=StartPos[0]
			ORFEnd=EndPos[0]
			
			for i in range(1,len(EndPos)):
				if StartPos[i]>ORFEnd:
					#print StartPos[i],EndPos[i], DNA[StartPos[i]:EndPos[i]+3]#,"ORFMarks",ORFStart,ORFEnd
					ORF_Count+=1
					ORF_Output.write(">ORF_%05d Strand: (+) Start: %d  End: %d\n" % (ORF_Count,StartPos[i],EndPos[i]+3))
					ORF_Output.write("%s\n" % DNA[StartPos[i]:EndPos[i]+3])
					ORFStart=StartPos[i]
					ORFEnd=EndPos[i]

			#Find ORFs on - Strand
			cDNA=ReverseComplement(DNA,"DNA")
			b=list(re.finditer(pattern,cDNA))

			StartPos=[a.start(0) for a in b]
			EndPos=[a.end(1) for a in b]
			
			i=0
			ORF_Count=ORF_Count+1
			ORF_Output.write(">ORF_%05d Strand: (-) Start: %d  End: %d\n" % (ORF_Count,len(DNA)-StartPos[i],len(DNA)-EndPos[i]+3))
			ORF_Output.write("%s\n" % cDNA[StartPos[i]:EndPos[i]+3])
			ORFStart=StartPos[0]
			ORFEnd=EndPos[0]

			for i in range(1,len(EndPos)):
				if StartPos[i]>ORFEnd:
					#print StartPos[i],EndPos[i], DNA[StartPos[i]:EndPos[i]+3]#,"ORFMarks",ORFStart,ORFEnd
					ORF_Count+=1
					ORF_Output.write(">ORF_%05d Strand: (-) Start: %d  End: %d\n" % (ORF_Count,len(DNA)-StartPos[i],len(DNA)-EndPos[i]+3))
					ORF_Output.write("%s\n" % cDNA[StartPos[i]:EndPos[i]+3])
					ORFStart=StartPos[i]
					ORFEnd=EndPos[i]


import re
import string
import argparse
import timeit
import numpy

#Get GenBank File, Output File
parser=argparse.ArgumentParser()
parser.add_argument('--input_filename',action='store') #
parser.add_argument('--output_filename',action='store') #
args=parser.parse_args()

source_file=args.input_filename
output_file=args.output_filename

n = 10
times = timeit.Timer(
    lambda: ORF_Find(args.input_filename,args.output_filename)
).repeat(repeat = n, number= 1)
for i in range(n):
	print "Run %03d: %.4f s" % (i+1,times[i])

print "Average Time: %.4f" % (sum(times)/float(n))
print "Median Time:  %.4f" % numpy.median(times)
print "Minimum Time: %.4f" % numpy.min(times)