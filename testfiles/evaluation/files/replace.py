from sys import argv

#usage : replace.py infile outfile
#replace char1 with char2

infile = argv[1]
outfile = argv[2]

seq= open(infile).read()

char1="E"
char2="m"

out = open(outfile,'w')
out.write(seq.replace(char1,char2))

out.close()

