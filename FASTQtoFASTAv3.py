from __future__ import division
import sys

#split FASTQ into chunks

def main(argv):
    

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
    
    outFile = open(sys.argv[2],'wr')
    #print outFile_Name
    with open(sys.argv[1]) as fastqFile:
        for fastq in fastqFile:
            if fastq:
                #print "LINE: ",str(lineNum)," ",fastq.rstrip('\r\n')
                if lineTrack == 4:
                    seqCount+=1
                    lineTrack=1
                    lineNum+=1
                    if outputCount == 1000000:
                        outputBuffer+=1000000
                        print "Converted ",str(outputBuffer-1),"sequences from FASTQ to FASTA "
                        outputCount = 1
                    else:
                        outputCount+=1
                elif lineTrack == 1:
                    parseFastq = fastq.split(' ')
                    makeHeader = []
                    makeHeader.append('>')
                    makeHeader.append(parseFastq[0][1:].rstrip('\r\n'))
                    outFile.write(''.join(makeHeader))
                    outFile.write('\n')
                    lineTrack+=1
                    lineNum+=1
                elif lineTrack == 2:
                    outFile.write(fastq.rstrip('\r\n'))
                    outFile.write('\n')
                    lineTrack+=1
                    lineNum+=1
                elif lineTrack == 3:
                    lineTrack+=1
                    lineNum+=1
               



if __name__ == "__main__": main(sys.argv)
