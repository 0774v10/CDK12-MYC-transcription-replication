#! /usr/bin/python

# BED generator of the genome in bins

import sys,os

def generate(genome_dictionary,bin,document):
	for chromosome in genome_dictionary.keys():
		maxlen=genome_dictionary[chromosome]
		remain= maxlen%bin 
		
		i=0
		while i<maxlen-remain:
			start=str(i)
			end=str(i+bin)
			f.write(chromosome+"\t"+start+"\t"+end+"\n")
			i+=bin
		start= str(i)
		end=str(i+remain)
		f.write(chromosome+"\t"+start+"\t"+end+"\n")


Choice=sys.argv[1]
bin=int(sys.argv[2])

mm9={"chr1":197195432,"chr2":181748087,"chr3":159599783,"chr4":155630120,"chr5":152537259,\
	 "chr6":149517037,"chr7":152524553,"chr8":131738871,"chr9":124076172,"chr10":129993255,\
	 "chr11":121843856,"chr12":121257530,"chr13":120284312,"chr14":125194864,"chr15":103494974,\
	 "chr16":98319150,"chr17":95272651,"chr18":90772031,"chr19":61342430,"chrX":166650296,"chrY":15902555,"chrM":16300}

hg19={"chr1":249250621,"chr2":243199373,"chr3":198022430,"chr4":191154276,"chr5":180915260,\
	 "chr6":171115067,"chr7":159138663,"chr8":146364022,"chr9":141213431,"chr10":135534747,\
	 "chr11":135006516,"chr12":133851895,"chr13":115169878,"chr14":107349540,"chr15":102531392,\
	 "chr16":90354753,"chr17":81195210,"chr18":78077248,"chr19":59128983, "chr20":63025520,\
	 "chr21":48129895,"chr22":51304566,"chrX":155270560,"chrY":59373566,"chrM":16571}

f=open("unsorted","w")

if Choice=="mm9":
	generate(mm9,bin,document=f)
elif Choice=="hg19":
	generate(hg19,bin,document=f)
else:
	print "Choice not valid.... exit..."
	sys.exit()

f.close()

os.system("sort -k1,1 -k2,2n unsorted > "+Choice+".bed")
os.remove("unsorted")

