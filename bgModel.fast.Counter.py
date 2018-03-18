import argparse
import os 
from collections import Counter

def revcomplE(word):
	revComp=''
	transdict={'A':'T','C':'G','G':'C','T':'A'}
	index=0
	while index <len(word):
		if word[index] == 'E' :#if is methylated C
			if index <len(word)-1 and word[index+1]=='G': #if G is next to it
				revComp += 'GE'
				index+=2   
			else:
				revComp += 'G'
				index+=1
		else:
			try: 
				revComp += transdict[word[index]]
			except KeyError:
				revComp += 'N'
			index+=1		 
	return revComp[::-1]

def revcomplEF(word):
	revComp=''
	transdict={'A':'T','C':'G','G':'C','T':'A','E':'F','F':'E'}
	index=0
	while index <len(word):
		try:
			revComp += transdict[word[index]]
		except KeyError:
			revComp += 'N'
		index += 1		 
	return revComp[::-1]

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("-gd", "--genomedir", help="a directory that contains genome fasta sequence files, ideally one chromosome per file")
	parser.add_argument("-k", "--klength",type =int, help="integer, the length of k-mers")
	parser.add_argument('--typeEF', dest='typeEF',action='store_true', help="using typeEF, default is typeE")
	parser.set_defaults(typeEF=False)

	args = parser.parse_args()
	

	genomedir = args.genomedir
	if genomedir == None:
		print "GENOMEDIR is missing, please specify"
		return
	try:
		k = int(args.klength)
	except:
		print "KLENGTH must be specified"
		return

	typeEF = args.typeEF
	print "GENOMEDIR is",genomedir
	print "klength is", k
	print "Using typeEF",typeEF

	kmercounts = Counter() # 
	total = 0
	for file in os.listdir(genomedir):
		print "Reading ...",file
		for line in open(genomedir+'/'+file):
			if line[0] == ">": # this is the name, skip
				continue
			seq = line.strip().upper()
			print len(seq)
			
			#count the kmers
			kmers = [seq[i:i+k] for i in range(len(seq)-k)]
			kmercounts.update(kmers)
			
			if typeEF: #check if type E or EF
				seq = revcomplEF(seq)
			else:
				print "revcompl E"
				seq = revcomplE(seq)

			kmers = [seq[i:i+k] for i in range(len(seq)-k)]
			kmercounts.update(kmers)

			

	#output
	total = 0

	outfile = ""
	if typeEF:
		outfile = "background_typeEF-"+str(k) +".tsv"
	else:
		outfile = "background_typeE-"+str(k) +".tsv"
	target = open(outfile,'w')
	target.write("TOTAL\t"+str(total)+'\n')
	for kmer in kmercounts:
		target.write(kmer+'\t'+str(kmercounts[kmer])+'\n')
	target.close()
		


	return 


if __name__ == "__main__":
	main()