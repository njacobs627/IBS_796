#GenBank Parser.py
import sys
import re
import string

def Complement(Sequence,Molecule):
	if Molecule=="DNA":
		return Sequence.translate(string.maketrans("ATCG","TAGC"))
	if Molecule=="RNA":
		return Sequence.translate(string.maketrans("AUCG","UAGC"))	

def ReverseComplement(Sequence,Molecule):
	return Complement(Sequence,Molecule)[::-1]

def Transcribe(Sequence):
	return Sequence.translate(string.maketrans("ATCG","UAGC"))[::-1]

def Translate(Sequence,Start):
	#Genetic Code Dictionary
	Genetic_Code={'AAA': 'K',
	'AAC': 'N',
	'AAG': 'K',
	'AAU': 'N',
	'ACA': 'T',
	'ACC': 'T',
	'ACG': 'T',
	'ACU': 'T',
	'AGA': 'R',
	'AGC': 'S',
	'AGG': 'R',
	'AGU': 'S',
	'AUA': 'I',
	'AUC': 'I',
	'AUG': 'M',
	'AUU': 'I',
	'CAA': 'Q',
	'CAC': 'H',
	'CAG': 'Q',
	'CAU': 'H',
	'CCA': 'P',
	'CCC': 'P',
	'CCG': 'P',
	'CCU': 'P',
	'CGA': 'R',
	'CGC': 'R',
	'CGG': 'R',
	'CGU': 'R',
	'CUA': 'L',
	'CUC': 'L',
	'CUG': 'L',
	'CUU': 'L',
	'GAA': 'E',
	'GAC': 'D',
	'GAG': 'E',
	'GAU': 'D',
	'GCA': 'A',
	'GCC': 'A',
	'GCG': 'A',
	'GCU': 'A',
	'GGA': 'G',
	'GGC': 'G',
	'GGG': 'G',
	'GGU': 'G',
	'GUA': 'V',
	'GUC': 'V',
	'GUG': 'V',
	'GUU': 'V',
	'UAA': 'stop',
	'UAC': 'Y',
	'UAG': 'stop',
	'UAU': 'Y',
	'UCA': 'S',
	'UCC': 'S',
	'UCG': 'S',
	'UCU': 'S',
	'UGA': 'stop',
	'UGC': 'C',
	'UGG': 'W',
	'UGU': 'C',
	'UUA': 'L',
	'UUC': 'F',
	'UUG': 'L',
	'UUU': 'F'}
	
	Codon_Pos=Start-1+3
	Codon=Sequence[Codon_Pos:Codon_Pos+3]
	Polypeptide=list('M')

	for i in range(0,len(Sequence)/3-1):
		Polypeptide.append(Genetic_Code[Codon])
		Codon_Pos=Codon_Pos+3
		Codon=Sequence[Codon_Pos:Codon_Pos+3]
	return (''.join(Polypeptide)).replace('stop','')

#Get GenBank File, Output File
source_file=sys.argv[1]
output_file=sys.argv[2]




with open(source_file,"r") as GenBank:
	Index_Line=""
	
	#Find Origin of DNA Sequence
	while "ORIGIN" not in Index_Line:
		Index_Line=GenBank.readline()
	
	#Feed sequence into a single string, DNA
	DNA=""
	DNA_Line=" "
	while len(DNA_Line)>0:
		DNA_Line=GenBank.readline()
		DNA=DNA+DNA_Line[10:len(DNA_Line)-1].replace(" ","")

	#Start finding annotations and translating
	GenBank.seek(0)
	Index_Line=""
	with open(output_file,"w") as Annotation:
		while "ORIGIN" not in Index_Line: #ORIGIN found at end of gene products
			
			while "/locus_tag" not in Index_Line: #Find Locus Tag
				if "ORIGIN" in Index_Line:
					break
				Index_Line=GenBank.readline()
			if "ORIGIN" in Index_Line:
				break		

			Locus_Tag=re.search('"([^"]*)"',Index_Line).group().replace('"','')
			
			while "CDS" not in Index_Line: #Find CDS Index
				Index_Line=GenBank.readline()
					
			Prime5=int(re.search('\d+\.',Index_Line).group().replace('.',''))
			Prime3=int((re.search('\.>*\d+',Index_Line).group().replace('.','')).replace('>',''))
			
			if "complement" in Index_Line:
				CDS=ReverseComplement(DNA[Prime5-1:Prime3].upper(),"DNA")
			else:
				CDS=DNA[Prime5-1:Prime3].upper()
			


			#Build Product Name string
			while "/product=" not in Index_Line: #Find Product Name Header
				Index_Line=GenBank.readline()
			
			Product_Name=(re.search('"([^"]*)',Index_Line).group().replace('"','')).replace('\n',' ') #First Line

			if Index_Line.count('"')<2: 		#Multiple lines in gene name
				Index_Line=GenBank.readline()

				while '"' not in Index_Line:	#Not done yet
					Product_Name=Product_Name+Index_Line.lstrip()
					Index_Line=GenBank.readline()

				#Final Line, containing "
				Product_Name=Product_Name+re.search(r'\S*([^"]*")',Index_Line).group().replace('"','').lstrip()
				Product_Name=(Product_Name.replace('\n',' ')).replace('- ','')


			Protein=Translate(Transcribe(ReverseComplement(CDS,"DNA")),1)
			print Locus_Tag,"   ",Product_Name

			Annotation.write(">%s product=%s\n" % (Locus_Tag,Product_Name)) #Header
			Annotation.write("%s\n\n" % Protein)
			