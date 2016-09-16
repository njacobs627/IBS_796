import random
import sys
import numpy
	
nuc=["A","C","G","T"]

def random_seq(gc,base_count):
	seq=[]
	GC_Count=0

	for i in range(1,base_count):
		nuc=["A","C","G","T"]
		newchar=numpy.random.choice(["A","T","C","G"],p=[(1-gc)/2,(1-gc)/2,gc/2,gc/2])
		seq.append(newchar)
		if newchar=="G" or newchar=="C":
			GC_Count=GC_Count+1
	
	return ''.join(seq)

GC=float(sys.argv[1])
BASE_NUM=int(sys.argv[2])
SEQ_NUM=int(sys.argv[3])
output_filename=sys.argv[4]

with open(output_filename,"write") as fh:
	for i in range(0,SEQ_NUM):
		newseq=random_seq(GC,BASE_NUM)
		fh.write(">sequence_%03d gc=%f length=%d\n" % (i+1,float((list(newseq).count("C")*1.0+list(newseq).count("G")*1.0)/BASE_NUM),BASE_NUM))
		fh.write("%s\n" % newseq)