import sys
import os
import time
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Seq import Seq
from Bio.Blast import NCBIXML

#Ryan P. Hillary, Stanford Medicine, rhillary@stanford.edu 2015

#TCR Sequencing Pipeline: V and J Region Identification for Multiplexed Samples.

#################THIS CODE IS DESIGNED AROUND A SPECIFIC ASSAY#######################
#######STANFORD MEDICINE CENTER FOR SLEEP SCIENCE AND MEDICINE TCR SEQUENCING########

#Previous Version: v1
#This program takes a raw merged FASTA file, searches for the 6mer sample barcode
#using a BLAST db, determines the V region sequence, and then determines the J
#region sequence in that order. If a SAMPLE ID is not found in the first 10bp, or
#if a suitable V or J region isnt found with high quality it aborts and moves to
#the next sequence.
#If a good sample barcode, J and V region have been found, BLAST supplies
#a start and end index for alignment. This codes takes the alignment
#index of the V and J and simply trims the given FASTA sequence. Specifically
#it takes the end index of the J and the beginning index of the V and trims
#all but the sequence lies between. Because the assay is multiplexed, we wish to cluster
#clones based on individual. The CDR3 region is paired with its sample ID and it
#becomes the key and entered into a dictionary of reads. After all reads have been
#processed it will output the dictionary values to an output file.

#Previous Version: v2
#v2 of this analysis incorperates completing the BLAST+ tasks on chunks of 100,000
#FASTA sequences at a time instead of individually. This greatly improved performance.
#Biopython parsing was abandoned due to its poor performance, a new parser was written.
#Attempts were made to balance speed and available resources (CPU,RAM etc)
#Specs for running 100bpx2 kits:
#v2 can process ~1,000,000 sequnces in ~20 mins.
#v1 required ~5+ days to complete a 100 mil sequence file divided into 20 chunks
#v1 required 20+ days to complete 100 mil sequences alone, no divisions of FASTA
#v2 requires ~2-3 hours to complete a 100 mil sequence file divided into 20 chunks
#v2 requires 1.5 days to complete 100 mil sequences alone, no divisions of the FASTA
#This script is designed as part of a pipeline of analyses, please note the inputs

#Previous Version: v3
#Designed to cope with both 100bpx2 and 250bpx2 Illumina Sequencing Approaches.
#v3 now requires a TCR C custom BLAST+ database input.
#Uses the TRAC/TRBC regions as an anchor to identify sample IDs.
#v2 assumed with the 100bpx2 approach that the sample would be from 4-9bp.
#Because the C region is an anchor another BLAST+ search was added.
#Various checks were added for accurate identification of the Sample ID and J region.
#v3 Is slightly slower than v2 mainly because of the added BLAST+ command.
#The cost in time allows a more versitile program to handle any approach.
#i.e: 100bpx2, 150bpx2 and 250bpx2 libraries rather than just 100bpx2.

#Previous Version: v4
#Added J_perc and V_perc to output
#Added flipping the strand on CDR3 region and outputing AA conv from C->F

#Current Version: v5
#Adjusted BLAST commands for V and J regions, slightly less stringent to allow
#for sequencing errors

#REQUIREMENTS:
#1- Biopython (pip is easiest, requries admin priv http://biopython.org/wiki/Download)
#2- BLAST+ Command Line Tools (requires admin priv, login as guest.
#2- ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ )
#2- BLAST+ excutables MUST be in the /usr/bin/ directory, if you type
#2- "which blastn" on a unix command line and it displays the path you are good to go.
#2- This code will function properly ONLY with BLAST+ 2.3.0+ or newer executables.
#3- Appropriate BLAST databases listed below.

#argv[1] input FASTA (MERGED AND QCed FASTQ -> converted to FASTA)
#argv[2] input V database (Custom BLAST+ database)
#argv[3] input J database (Custom BLAST+ database)
#argv[4] input barcode database (Custom BLAST+ database)
#argv[5] input identifiers csv (CSV dirctly from from SLEEP SERVER ["Primers" :)])
#argv[6] output file
#argv[7] chunk name string (number) (Typically the FASTA will be divided, this will be supplied by the sh script)
#argv[8] *NEWv3* input C database (Custom BLAST+ database)
#argv[9] *NEWv3* input study identifier i.e HS11 or HS92_TEST

def processBLAST(hitList,region):
    #This is for the VJ blast result
    #This assumes the highest hit must be 94% identical, 30bp long and we expect "minus" strand
    #We also gather the needed figures to calculate the V and J percentages for output.
    #This will return a python dictionary of good hits, the key is the FASTA seq name.
    
    temp_REGION = {}
    for hit in hitList:
        parsedHit = hit.split('\t')
        if hit:
            if region == "V":
                if float(parsedHit[5]) > 94 and int(parsedHit[6]) > 33 and parsedHit[13] == "plus":
                    if parsedHit[0] in temp_REGION:
                        #print "CRAP YO V" #Highly unlikely this would be satisfied
                        yo = 0
                    else:
                        makeValue = []
                        makeValue.append(parsedHit[1])
                        makeValue.append(str(parsedHit[14]))
                        makeValue.append(str(parsedHit[6]))
                        makeValue.append(str(parsedHit[10]))
                        makeValue.append(str(parsedHit[11]))
                        temp_REGION[parsedHit[0]] = str(','.join(makeValue))
            elif region == "J":
                if float(parsedHit[5]) > 94 and int(parsedHit[6]) > 28 and parsedHit[13] == "minus":
                    if parsedHit[0] in temp_REGION:
                        yo = 0
                    #print "CRAP YO J" #Highly unlikely this would be satisfied
                    else:
                        makeValue = []
                        makeValue.append(parsedHit[1])
                        makeValue.append(str(parsedHit[15]))
                        makeValue.append(str(parsedHit[6]))
                        makeValue.append(str(parsedHit[10]))
                        makeValue.append(str(parsedHit[11]))
                        temp_REGION[parsedHit[0]] = str(','.join(makeValue))
            elif region == "C":
                if float(parsedHit[5] > 70) and int(parsedHit[6]) > 15 and parsedHit[13] == "plus"):
                    if parsedHit[0] in temp_REGION:
                        yo = 0
                    #print "CRAP YO C" #Highly unlikely this would be satisfied
                    else:
                        makeValue = []
                        makeValue.append(parsedHit[1])
                        makeValue.append(parsedHit[10])
                        makeValue.append(parsedHit[11])
                        makeValue.append(parsedHit[14])
                        makeValue.append(parsedHit[15])
                        temp_REGION[parsedHit[0]] = str(','.join(makeValue))
    return temp_REGION

def processSAMPLE(hitList,C_REGIONS):
    #This is for the Sample/Barcode blast result  (Yes I could combine these two functions bla bla)
    #This ASSUMES the barcode will appear from 4 - 9 bp (1 beginning index) and will be "plus" strand
    #If the barcode passes QC it is added to a python dictionary with the FASTA name as the key
    #and the entire resulting dictionary is returned.
    
    #Get vital information from C Region
    
    temp_BARCODES = {}
    for hit in hitList:
        parsedHit = hit.split('\t')
        if hit:
            #Check if it fits to the closest match, indexwise.
            
            #print parsedHit[0]
            if parsedHit[0] in C_REGIONS:
                getCRegion = C_REGIONS.get(parsedHit[0])
                #print getCRegion
                #print hit
                parseC = getCRegion.split(',')
                #print parseC
                offset = int(parseC[1]) - 1
                adjustStart = int(parseC[3]) - offset
                sampleIDend = int(adjustStart) - 1
                #print "TCRC BEGIN: ",str(adjustStart)
                #print "SAMPLE END: ",str(sampleIDend)
                #print "SAMPLE INF: ",str(parsedHit[15])
                if int(parsedHit[15]) == int(sampleIDend) and parsedHit[13] == "plus":
                    #print "GOOD SAMPLE AND C REGION"
                    #print hit
                    if parsedHit[0] in temp_BARCODES:
                        print "CRAP YO S" #Two good barcodes that suit this? Unlikely but checking anyway.
                    else:
                        temp_BARCODES[parsedHit[0]]=parsedHit[1]
               #else:
                   #print "BAD----------------"
                   #print getCRegion
                   #print hit
    #print temp_BARCODES
    return temp_BARCODES

def makeOutFile(num,chunkName,studyName):
    #This handles making appropriate chunk FASTA file names.
    makeOut = []
    makeOut.append(str(studyName))
    makeOut.append("_temp")
    makeOut.append(str(chunkName))
    makeOut.append("_")
    makeOut.append("chunk")
    makeOut.append(str(num))
    makeOut.append(".fasta")
    return ''.join(makeOut)


def processBATCH(fasta,V,J,barcode,IDENTIFIERS,chunkNUM,C):
    #This is the main beast of the processing.
    #BLAST+ is called and the given .FASTA chunk file is handed to BLAST+, specifically to blastn.
    #Three separate BLAST runs are called for the Sample/Barcodes, the V and J regions.
    #After the BLAST runs are complete we open the FASTA chunk file and process each sequence
    #one by one until that chunk is exhuasted.
    #from Bio.Seq import Seq
    USABLE = 0
    TEMP_CDR3 = {}
    BARCODES = {}
    V_REGIONS = {}
    J_REGIONS = {}
    C_REGIONS = {}
    
    #BLAST run for the C Region, this will act as an anchor for identifying the correct barcode
    BLASTN_command_C =  NcbiblastnCommandline(query=fasta,db=C,num_threads=2,outfmt='"6 qseqid sseqid evalue bitscore length pident nident qframe sframe gaps sstart send sseq sstrand qstart qend qseq"',perc_identity=80,word_size=10,strand="plus",max_target_seqs=1)
    stdout_C, stderr_C = BLASTN_command_C()
    print "IDENTIFIED C in Chunk:        ",chunkNUM," in ",fasta
    parseBLAST_C = stdout_C.split('\n')
    #print stdout_C
    C_REGIONS = processBLAST(parseBLAST_C,"C")
    #print stdout_C
    print "PROCESSED C in Chunk:         ",chunkNUM," in ",fasta
    
    #Special Search for Short Reads, Align Sample Barcodes
    BLASTN_command_SAMPLE = NcbiblastnCommandline(query=fasta,db=barcode,max_target_seqs=5,num_threads=2,strand="plus",outfmt='"6 qseqid sseqid evalue bitscore length pident nident qframe sframe gaps sstart send sseq sstrand qstart qend qseq"',perc_identity=100,word_size=6)
    stdout_SAM, stderr_SAM = BLASTN_command_SAMPLE()
    print "IDENTIFIED BARCODE in Chunk:  ",chunkNUM," in ",fasta
    parseSAMPLE = stdout_SAM.split('\n')
    BARCODES = processSAMPLE(parseSAMPLE,C_REGIONS)#Filter and gather needed information
    #print stdout_SAM
    print "PROCESSED BARCODE in Chunk:   ",chunkNUM," in ",fasta
    
    #BLAST run for the J Region
    BLASTN_command_J = NcbiblastnCommandline(query=fasta,db=J,word_size=10,num_threads=2,max_target_seqs=1,perc_identity=90,outfmt='"6 qseqid sseqid evalue bitscore length pident nident qframe sframe gaps sstart send sseq sstrand qstart qend qseq"')
    stdout_J, stderr_J = BLASTN_command_J()
    print "IDENTIFIED J in Chunk:        ",chunkNUM," in ",fasta
    parseBLAST_J = stdout_J.split('\n')
    J_REGIONS = processBLAST(parseBLAST_J,"J")#Filter and gather needed information
    #print stdout_J
    print "PROCESSED J in Chunk:         ",chunkNUM," in ",fasta
    
    #BLAST run for the V Region
    BLASTN_command_V = NcbiblastnCommandline(query=fasta,db=V,word_size=10,perc_identity=90,num_threads=2,max_target_seqs=1, outfmt='"6 qseqid sseqid evalue bitscore length pident nident qframe sframe gaps sstart send sseq sstrand qstart qend qseq"')
    stdout_V, stderr_V = BLASTN_command_V()
    print "IDENTIFIED V in Chunk:        ",chunkNUM," in ",fasta
    parseBLAST_V = stdout_V.split('\n')
    V_REGIONS = processBLAST(parseBLAST_V,"V")#Filter and gather needed information
    #print stdout_V
    print "PROCESSED V in Chunk:         ",chunkNUM," in ",fasta
    
    currentLine = 1
    seq = ''
    header = ''
    
    #Read in and parse FASTA chunk file.
    #Yes, complexity is through the roof, this is mainly a proof of concept.
    with open(fasta) as fastaFile:
        for fasta_Seq in fastaFile:
            if fasta_Seq:
                if currentLine ==1:
                    header = fasta_Seq[1:].rstrip('\r\n')
                    currentLine+=1
                elif currentLine == 2:
                    seq = fasta_Seq.rstrip('\r\n')
                    currentLine = 1
                    
                    #TRIM THE CDR AND ADD TO/CHECK TEMP DICTIONARY
                    #The moment of truth: if its FASTA name returns from all four dictionaries
                    #we know we have a good hit so we move forward, if not we ignore.
                    if header in BARCODES and header in V_REGIONS and header in J_REGIONS and header in C_REGIONS:
                        USABLE+=1
                        #print str(C_REGIONS.get(header))
                        #print str(BARCODES.get(header))
                        #print str(J_REGIONS.get(header))
                        #print str(V_REGIONS.get(header))
                        
                        
                        
                        
                        get_BARCODE = BARCODES.get(header)
                        get_J = J_REGIONS.get(header)
                        parseJ = get_J.split(',')
                        get_V = V_REGIONS.get(header)
                        parseV = get_V.split(',')

                        jPer = float(parseJ[2])/float(parseJ[3])
                        vPer = float(parseV[2])/float(parseV[4])
                        #check beginning indexes of identifiers, adjust to convert to AA
                        
                        CDR3_Seq = seq[int(parseJ[1]):int(parseV[1])-1] #CDR3 trimmed here
                        adjust3_J = int(parseJ[1]) - (3 - int(parseJ[4]) + 1)
                        adjust3_V = int(parseV[1]) + (3 - int(parseV[3]) + 1)


                        CDR3_Seq_adj = seq[int(adjust3_J):int(adjust3_V) - 1]
                        CDR3_Seq_AA = Seq(CDR3_Seq_adj)
                        
                        getVhash = IDENTIFIERS.get(parseV[0])
                        getJhash = IDENTIFIERS.get(parseJ[0])
                        

                        #print "ADJUSTV: ",str(adjust3_V)
                        #print "AdJUSTJ: ",str(adjust3_J)
                        #print "BARCODE: ",get_BARCODE
                        #print "HEADER: ",header
                        #print "SEQ: ",seq
                        #print "J: ",parseJ
                        #print "V: ",parseV
                        #print "CDR3: ",CDR3_Seq
                        #print "CDR3_Adj: ",CDR3_Seq_AA.reverse_complement()
                        #print "TRANS_AA: ",CDR3_Seq_AA.reverse_complement().translate()
                        #print "HASHJ, ",getJhash
                        #print "HASHV: ",getVhash
                        #print "---------------------------------------------"
                        
                        #The main master dictionary of results is keyed with the CDR3 seq itself
                        #it is also paired with its sampleID number so we can identify per sample
                        #Key Structure:  <CDR3seq>_<sampleIDnum>
                        #I will be revising for more sensitivity so the key reflects:
                        #<CDR3seq>_<sampleIDnum>_<Vhash>_<Jhash>
                        #The Value in the dictionary is <sampleID>,<Vhash>,<CDR3>,<JHash>,<Count>
                        cloneID = []
                        cloneID.append(CDR3_Seq_adj)
                        cloneID.append(get_BARCODE)
                        cloneID.append(getVhash)
                        cloneID.append(getJhash)
                        cloneCheck = '_'.join(cloneID)
                        if cloneCheck in TEMP_CDR3:
                            updateFinal = TEMP_CDR3.get(cloneCheck)
                            parseFinal = updateFinal.split(',')
                            count = int(parseFinal[7]) + 1
                            parseFinal[7] = str(count)
                            TEMP_CDR3[cloneCheck] = str(','.join(parseFinal))
                        else:
                            makeFinal = []
                            makeFinal.append(get_BARCODE)
                            makeFinal.append(getVhash)
                            makeFinal.append(str(CDR3_Seq_AA.reverse_complement()))
                            makeFinal.append(getJhash)
                            makeFinal.append(str(vPer))
                            makeFinal.append(str(jPer))
                            makeFinal.append(str(CDR3_Seq_AA.reverse_complement().translate()))
                            makeFinal.append(str(1))
                            TEMP_CDR3[cloneCheck] = str(','.join(makeFinal))

#print "CHUNK PROCESSED in Chunk ",chunkNUM," in ",fasta
#Chunk has been processed, returning to be added into master dictionary.
    print "SUCCESSFULLY UTILIZED ",str(USABLE)," FASTA SEQUENCES."
    return TEMP_CDR3


#The Main Function.
def main(argv):
    
    CDR3 = {} #MASTER DICTIONARY
    IDENTIFIERS = {} #See argv[5]
    
    #This is a frequent function call, outputing a timestamp.
    print argv[1]," ",time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))
    
    
    seq_count = 0
    seq_total = 0
    #Read in and create a dictionary of Identifiers for teh VJ regions, used for vHash,jHash
    with open(argv[5]) as idFile:
        for id in idFile:
            parsedID = id.split(',')
            idHASH = str(parsedID[0])
            idNAME = str(parsedID[4])
            IDENTIFIERS[idNAME[1:-1]] = idHASH[1:-1]
    #print IDENTIFIERS
    #I decided against using Biopython to parse because of its poor performance.
    #I wrote my own parser for FASTA, will implement for FASTQ soon. This code will
    #make the code complexity nazis very angry. :)
    sequence = ''
    header=''
    chunkNum = 1
    lineNum=1
    seqCount=0
    seqBuffer=0
    chunkFile_Name = makeOutFile(chunkNum,argv[7],argv[9])
    chunkFile = open(chunkFile_Name,'wr')
    
    with open(sys.argv[1]) as fastaFile:
        for fasta in fastaFile:
            if fasta: #dumb check I know
                if seqCount == 999999: #This is the size of the chunk (199999 = 100000 FASTA seqs)
                    
                    #This Chunk file has reached its expected size, begin analysis for chunk.
                    chunkFile.write(fasta)
                    batchFile = open(chunkFile_Name,'r')
                    print "Chunk File bytes IN:  ",os.path.getsize(chunkFile_Name)," bytes in ",chunkFile_Name
                    
                    #Process Chunk
                    dictionary_return = processBATCH(chunkFile_Name,argv[2],argv[3],argv[4],IDENTIFIERS,chunkNum,argv[8])
                    
                    #Update Master dictionary with values from current batch.
                    for key, value in dictionary_return.iteritems():
                        if key in CDR3:
                            updateFinal = CDR3.get(key)
                            parseValue = value.split(',')
                            parseFinal = updateFinal.split(',')
                            count = int(parseFinal[7]) + int(parseValue[7])
                            parseFinal[7] = str(count)
                            CDR3[key] = str(','.join(parseFinal))
                        else:
                            CDR3[key] = value
                    print "Chunk File bytes OUT: ",os.path.getsize(chunkFile_Name)," bytes in ",chunkFile_Name
                    seqBuffer+=500000
                    
                    #Output Progress to User
                    print "Processed ",str(seqBuffer)," Sequences. File: ",chunkFile_Name," ",time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))
                    
                    #Make new chunk file, delete previous chunk file
                    chunkNum+=1
                    os.remove(chunkFile_Name)
                    chunkFile_Name = makeOutFile(chunkNum,argv[7],argv[9])
                    print "---------NEW CHUNKFILE: ",chunkFile_Name
                    chunkFile = open(chunkFile_Name,'wr')
                    seqCount = 0
                else:
                    #Keep writing to current chunk file.
                    chunkFile.write(fasta)
                    seqCount+=1
        else:
            #Cope with and process a partial file when the main FASTA file is exhuasted
            #print "CHUNKFILE: ",chunkFile_Name
            #print "EOF"
            print "Chunk File bytes IN:  ",os.path.getsize(chunkFile_Name)," bytes in ",chunkFile_Name
            
            #Process Chunk
            dictionary_return = processBATCH(chunkFile_Name,argv[2],argv[3],argv[4],IDENTIFIERS,chunkNum,argv[8])
            
            #Update Master dictionary with values from current batch.
            for key, value in dictionary_return.iteritems():
                #print value
                if key in CDR3:
                    
                    updateFinal = CDR3.get(key)
                    parseValue = value.split(',')
                    parseFinal = updateFinal.split(',')
                    count = int(parseFinal[7]) + int(parseValue[7])
                    parseFinal[7] = str(count)
                    CDR3[key] = str(','.join(parseFinal))
                else:
                    CDR3[key] = value
            print "Chunk File bytes OUT: ",os.path.getsize(chunkFile_Name)," bytes in ",chunkFile_Name
            os.remove(chunkFile_Name)

    #Master Dictionary Dump to Outfile. Analysis Complete.
    outFile = open(argv[6],'wr')
    for item, result in CDR3.iteritems():
        #print item
        #print result
        outFile.write(result)
        outFile.write('\n')
    print argv[1]," ",time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))
    print "**********DONE************  --->",argv[1]

#done

if __name__ == "__main__": main(sys.argv)



