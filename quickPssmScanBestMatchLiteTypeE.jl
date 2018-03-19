"""usage: julia quickPssmScanBestMatchLiteTypeE.jl memefile file.faa tagname resuldir"""
function score(kmer,PWM,background)
    #calculate log-odds scores for the matches
    pseudo=0.001
    s=0
    for i =1:endof(kmer)
        #s+=log(PWM[kmer[i]][i]+pseudo)
        s+=log(PWM[kmer[i]][i]+pseudo) - log(background[kmer[i]]+pseudo)
    end
    return s
end
function load_motifs(filename)
    file=open(filename)
    seq=split(strip(readstring(file)),"MOTIF")
    seq=seq[2:end];
    motifs=Dict()
    #ans=Dict{ASCIIString,Any}[]
    for s=1:endof(seq)
        t=split(strip(seq[s]),"\n")
        motifs[t[1]]=t[3:end] #skip the first 2 lines of the meme 
    end
    for m in keys(motifs)
        #println(m)
        tdict=Dict('A'=>Float64[],'C'=>Float64[],'G'=>Float64[],'T'=>Float64[],'E'=>Float64[])
        for pos=1:endof(motifs[m])
            #println(pos)
            tmp=split(strip(motifs[m][pos]),"\t")
            #print(tmp)
            push!(tdict['A'],float(tmp[1]))
            push!(tdict['C'],float(tmp[2]))
            push!(tdict['G'],float(tmp[3]))
            push!(tdict['T'],float(tmp[4]))
            push!(tdict['E'],float(tmp[5]))
        end 
        motifs[m]=tdict
        #ans[m]=tdict
    end
    return motifs
end 

function revcompl(word,transdict)
    #return word
    revComp=Char[]
    index=1
    while index <=endof(word)
        if word[index] == 'E' #if is methylated C
            #println("is E",index)
            if index < endof(word) && word[index+1]=='G' #if G is next to it
                append!(revComp,['G','E'])
                #println("found")
                index+=2   
            else
                #println("adding G, condition not met")
                append!(revComp,['G'])
                index+=1
            end 
        else
            #println("adding",transdict[word[index]])
            append!(revComp,[transdict[word[index]]])
            index+=1  
        end 
    end 
    return join(revComp[end:-1:1])
end


function main()
"""This function takes in a memefile, a set of sequences, and output the best score for each sequence by each motif in each separate file"""
	bg=Dict('A'=> 0.205, 'T' => 0.205, 'C' => 0.2874, 'G' => 0.295, 'E' => 0.0076) # this assumes about 25millions mC in the genome
    #p_cutoff = float(ARGS[1])
    motiffile = ARGS[1]
    sequencefile = ARGS[2]
    tag = ARGS[3]
    resultdir = ARGS[4] 

    bgfile =ARGS[5]
    seq=split(strip(readstring(open(bgfile))),'\n')
    bg=Dict()
    for line in seq
        tmp=split(strip(line),'\t')
        bg[collect((tmp[1]))[1]]=float(tmp[2])
    end 
    println(bg)
    #add this tag to all of the result files(at the end)
    #outfile =ARGS[3]
    #motiffile="./test.motif.meme"
    #sequencefile="./ENCFF002CIV.narrowPeak.500bps.faa"
    #negativeregions = ARGS[4]
    #output = ARGS[5] #the output are the name of the motif+ .motif
    motifs=load_motifs(motiffile)
    println("Loaded motifs")
    
    global transdict=Dict('A'=>'T','C'=>'G','G'=>'C','T'=>'A', 'N'=>'N')
    
    """
    posset = split(strip(readstring(open(sequencefile))),'>')[2:end]
    posseqs = ASCIIString[]
    posseqs_names=ASCIIString[]
    for seq in posset
        try
            tmp=split(strip(seq),"\n")
            if length(tmp)!=2
            	continue
            end
            push!(posseqs,tmp[2])
            push!(posseqs_names,tmp[1])
        catch e
            println(e)
            continue
        end 
    end
    """
    motifnames=keys(motifs)
    #println("motifs names",motifnames)
    for m in motifnames
        pwm=motifs[m]
        println("Scanning "*m)
        #outfile=m*".scanned.txt"
        target=open(resultdir*"/"*m*"."*tag*".scanned.txt","w")
        write(target,"MOTIF\t"*m*"\n")
        #scanning positive seqs
        counts=0
        k=length(pwm['A'])
        f = open(sequencefile)
        header =""
        #for i=1:length(posseqs) #for each line
        for line in eachline(f)
            if line[1]=='>'
                header=strip(line)
                
            else
                #seq = posseqs[i] #seq=strip(ln)
                seq=strip(line)
                seqname=header
                #seqname = posseqs_names[i] #seqname=header
                revcomplseq=revcompl(seq,transdict)
                sequencestoscan=[seq;revcomplseq]
                #println(seq)
                #println(revcomplseq)
                #print sequencestoscan
                bestscore=-100
                for sequence in sequencestoscan
    	            if length(sequence)<k
    	                continue
    	            else 
    	                #tmp = collect(seq)
                        #println(header)
                        #println(sequence)
    	                for i=1:length(sequence)-k
    	                    kmer = sequence[i:i+k-1]
                            #println(kmer) 
    	                    if 'N' in kmer
    	                   		continue
                                #println(kmer,"errooor")
                            else
                                s=score(kmer,pwm,bg)
                                if s > bestscore
                                    bestscore=s
                                end
                            end
    	                end
    	            end 
                #println(bestscore)
                end
            #println(seqname,'\t',bestscore)
                write(target,seqname*"\t"*string(bestscore)*"\n")
            end
        #println("TOTAL\t"*string(counts)*" matches\n")
        end
        close(target)
    end 
end

@time(main())
