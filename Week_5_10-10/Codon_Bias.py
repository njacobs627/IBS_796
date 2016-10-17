### Codon Bias of multi-FASTA

def BaseToNumber(Base):
	Base_Code={'A':0,
	'C':1,
	'G':2,
	'T':3}
	return Base_Code[Base]

def StringToNumber(Sequence):
	if len(Sequence)==0:
		return 0
	else:
		prefix=Sequence[0:(len(Sequence)-1)]
		suffix=Sequence[(len(Sequence)-1):len(Sequence)]
		return 4*StringToNumber(prefix)+BaseToNumber(suffix)

def NumberToString(Index,k):
	if k==1:
		return NumberToBase(Index)
	else:
		prefixIndex=Index/4
		r=Index % 4
		return NumberToString(prefixIndex,k-1) + NumberToBase(r)

def NumberToBase(Number):
	Base=['A','C','G','T']
	return Base[Number]

def Codon_Bias(Input,Output):
	Gene=' ' #Dummy Starting Value, overwritten first
	Codon_Count=[0 for x in range(64)]
	Codon_Number=0

	with open(Input,"r") as fasta:
		while Gene!='':
			#Parse Gene
			fasta.readline()	#Ignore header
			Gene=fasta.readline()
			
			#Increment Number of Codons
			Gene_Codons=len(Gene.replace('\n',''))/3
			Codon_Number=Codon_Number+Gene_Codons

			for i in range(Gene_Codons):
				Codon=Gene[(i*3):(i*3+3)]
				Codon_Count[StringToNumber(Codon)]+=1

	with open(Output,"w") as output:
		for i in range(64):
			output.write("%s %f" % (NumberToString(i,3), Codon_Count[i]*1.0/Codon_Number))
			#print "%s %f" % (NumberToString(i,3), Codon_Count[i]*1.0/Codon_Number)   #Option to Print

import argparse
import numpy
import time

parser=argparse.ArgumentParser()

parser.add_argument('--input_filename',action='store') #
parser.add_argument('--output_filename',action='store') #

args=parser.parse_args()
start=time.time()
Codon_Bias(args.input_filename,args.output_filename)
end=time.time()

print 'Runtime: %4.4f s' % float(end-start)