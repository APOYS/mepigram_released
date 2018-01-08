#!/usr/bin/env python
"""
TODO: 

Non-CPG function

add basecomposition for scanning

Workflow:

    - Takes a set of sequences
    - Di-nuc shuffles it .
    - Takes the background model
    - Takes the number of of maximum motifs

"""
from sys import argv
import random as rd
import os
import sys
def load_motifs_typeE(filename):
    file=open(filename)
    seq=file.read().split("MOTIF")
    seq=seq[1:]
    motifs={}
    infos={}
    for s in range(len(seq)):
        t=seq[s].strip().split("\n")
        name=int(t[0].split('_')[0])
        motifs[name]=t[2:]
        infos[name]=t[0:2]
    for m in motifs:
        tdict={'A':[],'C':[],'G':[],'T':[],'E':[],}
        for pos in range(len(motifs[m])):
            tmp=motifs[m][pos].strip().split("\t")
            tdict['A']+=[float(tmp[0])]
            tdict['C']+=[float(tmp[1])]
            tdict['G']+=[float(tmp[2])]
            tdict['T']+=[float(tmp[3])]
            tdict['E']+=[float(tmp[4])]
        motifs[m]=tdict
    
    return motifs,infos

def load_motifs_typeEF(filename):
    file=open(filename)
    seq=file.read().split("MOTIF")
    seq=seq[1:]
    motifs={}
    infos={}
    for s in range(len(seq)):
        t=seq[s].strip().split("\n")
        name=int(t[0].split('_')[0])
        motifs[name]=t[2:]
        infos[name]=t[0:2]
    for m in motifs:
        tdict={'A':[],'C':[],'G':[],'T':[],'E':[],'F':[]}
        for pos in range(len(motifs[m])):
            tmp=motifs[m][pos].strip().split("\t")
            tdict['A']+=[float(tmp[0])]
            tdict['C']+=[float(tmp[1])]
            tdict['G']+=[float(tmp[2])]
            tdict['T']+=[float(tmp[3])]
            tdict['E']+=[float(tmp[4])]
            tdict['F']+=[float(tmp[5])]
        motifs[m]=tdict
    
    return motifs,infos

def taggingmotifs(filename,outfile):
    print "tagging m-motifs in",filename
    '''This function adds a tag to motifs with P(E) >= 0.5'''
    #filename="./test_mepigram_pipeline_complete.meme"
    #outfile=filename.replace(".meme",".tagged.meme")
    motifs,infos=load_motifs_typeE(filename)
    header='''MEME version 4.5 - modififed
ALPHABET= ACGTE 
strands: + 
Background letter frequencies
A 0.295 C 0.205 G 0.205 T 0.295 E 0.0076

'''
    #print infos
    #print the motifs with the E or non-E or EF motif tag
    target=open(outfile,'w')
    target.write(header)
    threshold=0.5
    #print header
    for m in motifs:
        modified=False
        for pos in motifs[m]['E']:
            if pos > threshold:
                modified=True
                break
        firstline=''
        if modified==True:
            firstline="MOTIF"+'\t'+infos[m][0].strip()+"_m-motif"
        else:
            #print m,"motif"
            firstline="MOTIF"+'\t'+infos[m][0].strip()
        target.write(firstline+'\n')
        target.write(infos[m][1]+'\n')
        for i in range(len(motifs[m]['A'])):
            line=[]
            for j in ['A','C','G','T','E']:
                line+=[motifs[m][j][i]]
            line='\t'.join([str(x) for x in line])
            target.write(line+'\n')
    target.close()
    print "Finished tagging"
    return

def main():
    #parsing argument
    faafile = None
    #memefile = None
    outfile = None #this contains the meme file, the scan file if chosen.
    mode = None #right now use the default mode, cpg
    backgroundfile = None
    graphdir = None #this is hardcoded for now (using only the cpg 8-mer graph)
    filter_boolean = False
    maxmotifnum = None
    enrichmentmode = "none"  #choose whether to calculate the enrichment of all motifs found or just the m-motifs, or just non-m motifs, or none

    ####Testing using hardcodes

    #faafile="../epigram-0.003/Sox2-P.faa"
    #faafile="../CrossScan_H1_CTCF/ENCFF002CDS.narrowPeak.faa"
    #graphdir="./metgraph-8mer/"
    #backgroundfile="./hg19_H1_meth_bg/background_met-8.tsv"
    maxmotifnum=200
    mode=None
    seed=rd.randint(0, 10000)

    
    usage = """USAGE: 
    %s [options] -f fastafile -m mode -b backgroundmodel 
        REQUIRED:

        -f <filename>       input file, FASTA format 
        -m typeE|typeEF     input mEpigram mode: choose 'typeE' or 'typeEF'
        -b <filename>       background model, shows the distribution of kmers in the genome
        -g <dirname>        graph directory

        OPTIONAL:
        -o <filename>       name of the output file to be created.
        -n <number> maximum number of motifs to be output, default 200.  
        -h                  print this usage message
    """ % (sys.argv[0])

    # parse command line
    i = 1
    while i < len(sys.argv):
        arg = sys.argv[i]
        if (arg == "-f"):
            i += 1
            try: faafile = sys.argv[i]
            except:
            	print " Error in -f" 
            	sys.exit(1)
        elif (arg == "-m"):
            i += 1
            try: mode = sys.argv[i]
            except:
            	print "Error in -m" 
            	sys.exit(1)
        elif (arg == "-g"):
            i += 1
            try: graphdir = sys.argv[i]
            except:
                print "Error in -g" 
                sys.exit(1)
        elif (arg == "-b"):
            i += 1
            try: backgroundfile = sys.argv[i]
            except:
            	print "Error in -b" 
            	sys.exit(1)
        elif (arg == "-o"):
            i += 1
            try: outfile = sys.argv[i]
            except:
            	print "Error in -o" 
            	sys.exit(1)
        elif (arg == "--filter"):
            i+=1
            filter_boolean = True
        elif (arg == "-n"):
            i += 1
            try: maxmotifnum = int(sys.argv[i])
            except:
                print "Error in -n" 
                sys.exit(1)
        elif (arg == "-e"):
            i += 1
            try: enrichmentmode = sys.argv[i]
            except:
                print "Error in -e" 
                sys.exit(1)
        elif (arg == "-h"):
        	print >> sys.stderr, usage; sys.exit(1)
        else:
            print >> sys.stderr, "Unknown command line argument: " + arg + " \nPlease check the usage with -h"
            sys.exit(1)
        i += 1

    # check that required arguments given
    if (faafile == None):
    	print "No input fasta file, exit."
    	print >> sys.stderr, usage; sys.exit(1)
    	sys.exit(1)
    elif (mode != 'typeE' and mode !='typeEF'):
    	print "Error in mode, please choose a valid mode."
        print >> sys.stderr, usage; sys.exit(1)
    	sys.exit(1)
    elif (enrichmentmode!="none" and enrichmentmode!="all" and enrichmentmode!="meth" and enrichmentmode!="nonmeth"):
        print "Invalid option in -e, exit."
        print >> sys.stderr, usage; sys.exit(1)
        sys.exit(1)
    if (outfile == None):
        outfile = faafile+".mepigram.meme"
    	print "No output file specified, will use the default",outfile
    if (graphdir == None):
        print "No graph directory specified, please input graph directory"
        sys.exit(1)
    if (backgroundfile == None):
    	print "No background file specified, please input your background"
        sys.exit(1)
        #background={'A':0.2,'C':0.2,'G':0.2,'T':0.2,'E':0.2}

    # Now load the data:
    # shuffle the data
    print "Reading fasta file"
    faaname=faafile
    faafile=open(faafile).read().strip().split(">")
    faafile=faafile[1:]
    faaseqs={}
    totalbasenum=0
    #remove lines with less than a certain length 
    badcount=0
    minlen=20
    for line in faafile:
        tmp=line.strip().split("\n")
        seq=''.join(tmp[1:])
        if len(seq)>=minlen:
            faaseqs[tmp[0]]=seq
            totalbasenum+=len(seq)
        else:
            badcount+=1
    if badcount>0:
        print badcount,"sequences shorter than "+str(minlen)+" base pairs. They are skipped."
    print "Number of sequences:",len(faaseqs),". Number of total bases:",totalbasenum
    
    estimatekmernum=totalbasenum-len(faaseqs)*10 #assuming kmer length is 10
    timestoshuffle=1
    minkmernumber=1 # Maybe use 20 millions here later, testing this out
    if estimatekmernum<minkmernumber:
        timestoshuffle=minkmernumber/estimatekmernum+1
        print "Less than "+str(minkmernumber)+" kmers, run Di-nuc shuffling",timestoshuffle,"times"
    else:
    	#print "More than 1000000"
    	print "Di-nuc shuffling 1 times"
    #run di-nuc shuffling
    shufflefile=faaname+'.'+str(seed)+".DS.tmp.faa"
    command=''
    if mode == 'typeE':
        command="python fasta-dinucleotide-shuffle_typeE.py -s "+str(seed)+" -f "+faaname+" -c "+str(timestoshuffle)+" > "+shufflefile
    else:
        command="python fasta-dinucleotide-shuffle_typeEF.py -s "+str(seed)+" -f "+faaname+" -c "+str(timestoshuffle)+" > "+shufflefile

    print command
    os.system(command)

    print "Running mEpigram... "
    if mode == 'typeE':
        command="python mepigram_typeE.py "+faaname+" "+shufflefile+" "+backgroundfile+" "+graphdir +" "+outfile+" "+str(maxmotifnum)
    else:
        command="python mepigram_typeEF.py "+faaname+" "+shufflefile+" "+backgroundfile+" "+graphdir +" "+outfile+" "+str(maxmotifnum)
    print command
    os.system(command)

    '''This part renames the motifs into meth and unmeth motifs'''
    memefile=outfile#+'.tagged'
    #taggingmotifs(outfile,memefile)
    #memefile="../ENCFF002CQR.2motifs.meme"

    resultdir='/'.join(memefile.split('/')[:-1])
    resultdir=memefile+".results"
    os.system("rm -r "+resultdir)
    os.system("mkdir "+resultdir)
    

    """
    Filtering step will go in here
    Combine all of the motifs together into a list, then calculate pairwise similarity distance of all of them. 
    If two motifs have similarity higher than a certain threshold, keep one, discard one. 
    """
    #print "Filtering results..."


    """
    Calculate enrichment of motifs: 2 methods
    
        * Fisher p-value: <Use this one for now>
            - Scan through to get the max score of each region
            - Determine a cutoff such that the p-value is maximized
        * Straight enrichment (Vu's method)
            - Get a score distribution of k-mers for each motif by scanning through randomly shuffled regions
            - Determine a cutoff accordingly to a p-value
            - Use the cutoff to call for matches.
    """
    print "Calculating enrichment by scanning..."
    # use method 1,

    print "Estimating background base compostion..."
    baseCompositionFile=backgroundfile+".basecomposition.tmp"
    command="python baseComposition.py "+backgroundfile+" "+baseCompositionFile
    print command
    os.system(command)
    
    if mode == 'typeE':
        print "Scanning on positive sequences..."
        command="julia quickPssmScanBestMatchLiteTypeE.jl "+memefile+" "+faaname+" "+"quickscan.positive.tmp"+" "+resultdir+" "+baseCompositionFile
        print command
        os.system(command)	
        print "Scanning on negative sequences..."
        command="julia quickPssmScanBestMatchLiteTypeE.jl "+memefile+" " +shufflefile+" "+"quickscan.negative.tmp"+" "+resultdir+" "+baseCompositionFile
        print command
        os.system(command)

    else:
        print "Scanning on positive sequences..."
        command="julia quickPssmScanBestMatchLiteTypeEF.jl "+memefile+" "+faaname+" "+"quickscan.positive.tmp"+" "+resultdir+" "+baseCompositionFile
        print command
        os.system(command)  
        print "Scanning on negative sequences..."
        command="julia quickPssmScanBestMatchLiteTypeEF.jl "+memefile+" " +shufflefile+" "+"quickscan.negative.tmp"+" "+resultdir+" "+baseCompositionFile
        print command
        os.system(command)

    #sdos.system("rm "+shufflefile)
    print "Calculating Fisher P-values..."
    #get the result files together
    scannedresults={}
    for file in os.listdir(resultdir):
        if "quickscan.positive" in file:
            motif=file.split('.quickscan')[0]
            if motif not in scannedresults:
                scannedresults[motif]={}
                scannedresults[motif]['P']=file
            else:
                if "P" in scannedresults[motif]:
                    "ERROR!!!: Already found positive scan file for motif",motif
                else:
                    scannedresults[motif]['P']=file
            
        elif "quickscan.negative" in file:
            motif=file.split('.quickscan')[0]
            if motif not in scannedresults:
                scannedresults[motif]={}
                scannedresults[motif]['N']=file
            else:
                if "N" in scannedresults[motif]:
                    print "ERROR!!!: Already found negative scan file for motif",motif
                else:
                    scannedresults[motif]['N']=file
        else:
            continue
    print "Found scanned results for",len(scannedresults),"motifs"

    fisherresults=[]
    for motif in scannedresults:
        print motif
        if len(scannedresults[motif])!=2:
            print "ERROR!!!: number of scanned files for motif",motif,"is not 2"
        command="python2.7 fisher_P-value.py "+resultdir+'/'+scannedresults[motif]['P']+" "+resultdir+'/'+scannedresults[motif]['N']+" "+resultdir+'/'+motif+".fisher.tmp" +" "+motif
        #print command
        os.system(command)
        fisherresults+=[open(resultdir+'/'+motif+".fisher.tmp").read().strip()]

    print "Calculating enrichments..."
    out=open(resultdir+'/enrichments.tsv','w')
    out.write("MOTIF\tp-value\tscoreCutoff\tPosMatches\tPosNonMatches\tNegMatches\tNegNonMatches\tEnrichment"+'\n')
    for line in fisherresults:
        out.write(line+'\n')
    out.close()
    #concat the fisher results


    os.system("mv "+memefile+" "+resultdir)
    os.system("mv "+outfile+" "+resultdir)
    print "Cleaning up temporary files..."
    os.system("rm "+shufflefile)
    #ftoremove=[]
    for file in os.listdir(resultdir):
        if "tmp" in file:
            command="rm "+resultdir+'/'+file
            print command
            os.system(command)

    
if __name__ == "__main__":
    main()



