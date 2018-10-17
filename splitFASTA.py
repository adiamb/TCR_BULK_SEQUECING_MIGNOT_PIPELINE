from __future__ import division
import sys

#split FASTA into chunks

def makeOutFile(num,study):
    makeOut = []
    makeOut.append(study)
    makeOut.append("_")
    makeOut.append("chunk")
    makeOut.append(str(num))
    makeOut.append(".fasta")
    return ''.join(makeOut)

def main(argv):
    
    divides = sys.argv[1]
    lines = sys.argv[2]
    study = sys.argv[3]
    
    #print lines
    #print divides
    
    sequences = float(lines) // 2
    perFile = float(sequences) // float(divides)
    remainder = float(sequences) % float(divides)
    seqChunk = perFile*2
    #print str(perFile)
    #print str(remainder)
    #print str(seqChunk)
    
    fileNum = 1
    lineNum = 1
    lineTrack = 1
    chunkTrack = 1
    chunksMade = 1
    seqCount = 1
    #makeFirst OutputFile
    outputCount = 1
    outputBuffer = 1
    outputMake = []
    
    outFile_Name = makeOutFile(fileNum,study)
    outFile = open(outFile_Name,'wr')
    #print outFile_Name
    with open(sys.argv[4]) as fastqFile:
        for fastq in fastqFile:
            if fastq:
                #print "LINE: ",str(lineNum)," ",fastq.rstrip('\r\n')
                if lineTrack == 2:
                    #print "Seq"
                    #print "SeqCount: ",str(seqCount)
                    outFile.write(fastq)
                    seqCount+=1
                    lineTrack=1
                    lineNum+=1
                    if outputCount == 1000000:
                        outputBuffer+=1000000
                        print "Processed ",str(outputBuffer-1)," FASTA seqs. ",study
                        outputCount = 1
                    else:
                        outputCount+=1

                    if chunkTrack == perFile:
                        fileNum+=1
                        newFile_Name = makeOutFile(fileNum,study)
                        outFile = open(newFile_Name,'wr')
                        #print newFile_Name
                        #print "CHUNK-------- ",str(chunksMade)
                        chunksMade+=1
                        chunkTrack=1
                    else:
                        chunkTrack+=1

                else:
                    outFile.write(fastq)
                    #print "OLD"
                    lineTrack+=1
                    #chunkTrack+=1
                    lineNum+=1
#write to file


if __name__ == "__main__": main(sys.argv)
