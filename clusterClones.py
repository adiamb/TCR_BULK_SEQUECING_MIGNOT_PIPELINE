import glob
import sys
import time

study = sys.argv[1]

makePattern = []
makePattern.append(study)
makePattern.append("_chunk")
makePattern.append("*")
makePattern.append(".rcl")
pattern = ''.join(makePattern)

print "Begin Clumping"

CLONES = {}

for name in glob.glob(pattern):
    #print "Detected: ",name
    print "Clustering: ",name ," ",time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))
    with open(name) as rclFile:
        for rcl in rclFile:
            if rcl:
                #Make Dictionary Key
                #22,132,GTACTGCTCGTAAGCGAAGCTGCTGGC,8,1.0,1.0,9
                stripRcl = rcl.rstrip('\r\n')
                parsedRcl = stripRcl.split(',')
                makeKey = []
                makeKey.append(parsedRcl[2])
                makeKey.append(parsedRcl[0])
                makeKey.append(parsedRcl[1])
                makeKey.append(parsedRcl[3])
                if '_'.join(makeKey) in CLONES:
                    #print "CLONE"
                    getHit = CLONES.get('_'.join(makeKey))
                    #print getHit
                    parseHit = getHit.split(',')
                    updateCount = int(parseHit[7]) + int(parsedRcl[7])
                    parsedRcl[7] = str(updateCount)
                    CLONES['_'.join(makeKey)] = ','.join(parsedRcl)
                else:
                    CLONES['_'.join(makeKey)] = stripRcl
    print "Processed:  ",name," ",time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))


outMake = []
outMake.append(study)
outMake.append("_CLUSTERED.rcl")
outName = ''.join(outMake)
outFile = open(outName,'wr')
for item,result in CLONES.iteritems():
    outFile.write(result)
    outFile.write('\n')
                

