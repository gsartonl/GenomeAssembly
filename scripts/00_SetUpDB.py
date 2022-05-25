#!/bin/pyton3
import sys
import os
import glob

workDir=sys.argv[1]
pathRawFiles=sys.argv[2]

def getConfigGenomes(pathRawFiles,workDir) :
    os.chdir(pathRawFiles)
    fwd = glob.glob("*_1.fq.gz")
    os.chdir(workDir)
    sampleName = [f.split('_')[0] for f in fwd]
    hashm=['_'.join(f.split('_')[1:3]) for f in fwd]
    lane=[f.split('_')[3] for f in fwd]


    with open(workDir +"/setup/config2.yaml", 'w') as outFile:
        toWrite="SAMPLES:\n"
        for s in range(len(fwd)) :
            toWrite+= "  - " + sampleName[s] + "\n"

        toWrite+="LIBRARY:\n"
        for s in range(len(fwd)) :
            toWrite+= "  " + sampleName[s] + ":\n" + "    Lane: '" + lane[s] + "'\n" + "    Hash_FWD: '" + hashm[s] + "'\n" + "    Hash_REV: '" + hashm[s] + "'\n"

        print("Making config file for samples")
        outFile.write(toWrite)



getConfigGenomes(pathRawFiles,workDir)
