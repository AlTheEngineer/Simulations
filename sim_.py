import sys
import gc
from Aptamers import Aptamers
from Selection import Selection
from Amplification import Amplification
from Mutation import Mutation
from postprocess import dataAnalysis
import utils

aptamerType = str(sys.argv[1])
aptamerNum = int(sys.argv[2])
aptamerSeq = str(sys.argv[3])
seqLength = int(sys.argv[4])
selectionThreshold = long(sys.argv[5])
distanceMeasure = str(sys.argv[6])
roundNum = int(sys.argv[7])
pcrCycleNum = int(sys.argv[8])
pcrYield = float(sys.argv[9])
pcrErrorRate = float(sys.argv[10])
stringency = float(sys.argv[11])
samplingSize = int(sys.argv[12])
outputFileNames = str(sys.argv[13])
post_process = bool(sys.argv[14])

# Instantiating classes
Apt = Aptamers()
S = Selection()
Amplify = Amplification()
Mut = Mutation()

if(distanceMeasure == "hamming"):
# SELEX simulation based on random aptamer assignment, hamming-based definite selection, and
# non-ideal stochastic amplfication with no bias. 
    for r in range(roundNum):
        if(r==0):
            if(aptamerType == 'DNA'):
                alphabetSet = 'ACGT'
                aptamerSeqs, initialSeqNum = Apt.optimumAptamerGenerator(aptamerNum, alphabetSet, seqLength)
            elif(aptamerType == 'RNA'):
                alphabetSet = 'ACGU'
                aptamerSeqs, initialSeqNum = Apt.optimumAptamerGenerator(aptamerNum, alphabetSet, seqLength)
            else:
                print("Error: Simulation of %.s aptamers not supported" %aptamerType)
                break
            if(aptamerNum == 0):
                aptamerSeqs = aptamerSeq
            print("optimum sequences have been chosen")
            print("SELEX Round 1 has started")
            print("total number of sequences in initial library = "+str(initialSeqNum))
            slctdSeqs = S.stochasticHammingSelection_initial(alphabetSet, seqLength, aptamerSeqs, selectionThreshold, initialSeqNum, samplingSize, outputFileNames, r, stringency)
            print("selection carried out for R0")
            amplfdSeqs = Amplify.randomPCR_with_ErrorsAndBias_FASTv2(slctdSeqs, seqLength, pcrCycleNum, pcrYield, pcrErrorRate, aptamerSeqs, alphabetSet, distanceMeasure)
            del(slctdSeqs)
            gc.collect()
            print("Amplification carried out for R1")
            outFile = outputFileNames + "_R" + str(r+1)
            nxtRnd = open(outFile, 'w')
            print("writing R1 seqs to file")
            for seqIdx in amplfdSeqs:
                seq = Apt.pseudoAptamerGenerator(seqIdx, alphabetSet, seqLength)
                nxtRnd.write(str(seq)+'\t'+str(int(amplfdSeqs[seqIdx][0]))+'\t'+str(int(amplfdSeqs[seqIdx][1]))+'\t'+'\n') #write seqIdx, count, distance, and bias...for now
            nxtRnd.close()
        else:
            print("SELEX Round "+str(r+1)+" has started")
            totalSeqNum, uniqSeqNum = utils.seqNumberCounter(amplfdSeqs)
            print("total number of sequences in initial pool = "+str(totalSeqNum))
            print("total number of unique sequences in initial pool = "+str(int(uniqSeqNum)))
            slctdSeqs = S.stochasticHammingSelection(alphabetSet, seqLength, amplfdSeqs, selectionThreshold, uniqSeqNum, totalSeqNum, samplingSize, outputFileNames, r, stringency)
            del(amplfdSeqs)
            gc.collect()
            print("Selection carried for R"+str(r+1))
            amplfdSeqs = Amplify.randomPCR_with_ErrorsAndBias_FASTv2(slctdSeqs, seqLength, pcrCycleNum, pcrYield, pcrErrorRate, aptamerSeqs, alphabetSet, distanceMeasure)
            del(slctdSeqs)
            gc.collect()
            print("Amplification carried for R"+str(r+1))
            print("writing R"+str(r+1)+" seqs to file")
            outFile = outputFileNames + "_R" + str(r+1)
            with open(outFile, 'w') as f:
                for seqIdx in amplfdSeqs:
                    seq = Apt.pseudoAptamerGenerator(seqIdx, alphabetSet, seqLength)
                    f.write(str(seq)+'\t'+str(int(amplfdSeqs[seqIdx][0]))+'\t'+str(int(amplfdSeqs[seqIdx][1]))+'\n')
    print("SELEX completed")

elif(distanceMeasure == "basepair"):
# SELEX simulation based on random aptamer assignment, hamming-based definite selection, and
# non-ideal stochastic amplfication with no bias. 
    for r in range(roundNum):
        if(r==0):
            if(aptamerType == 'DNA'):
                alphabetSet = 'ACGT'
                aptamerSeqs, initialSeqNum = Apt.optimumAptamerGenerator(aptamerNum, alphabetSet, seqLength)
            elif(aptamerType == 'RNA'):
                alphabetSet = 'ACGU'
                aptamerSeqs, initialSeqNum = Apt.optimumAptamerGenerator(aptamerNum, alphabetSet, seqLength)
            else:
                print("Error: Simulation of %.s aptamers not supported" %aptamerType)
                break
            if(aptamerNum == 0):
                aptamerSeqs = aptamerSeq
            print("optimum sequences have been chosen")
            print("SELEX Round 1 has started")
            print("total number of sequences in initial library = "+str(initialSeqNum))
            slctdSeqs = S.stochasticBasePairSelection_initial(alphabetSet, seqLength, aptamerSeqs, selectionThreshold, initialSeqNum, samplingSize, outputFileNames, r, stringency)
            print("selection carried out for R0")
            amplfdSeqs = Amplify.randomPCR_with_ErrorsAndBias_FASTv2(slctdSeqs, seqLength, pcrCycleNum, pcrYield, pcrErrorRate, aptamerSeqs, alphabetSet, distanceMeasure)
            print("Amplification carried out for R1")
            outFile = outputFileNames + "_R" + str(r+1)
            nxtRnd = open(outFile, 'w')
            print("writing R1 seqs to file")
            for seqIdx in amplfdSeqs:
                seq = Apt.pseudoAptamerGenerator(seqIdx, alphabetSet, seqLength)
                nxtRnd.write(str(seq)+'\t'+str(int(amplfdSeqs[seqIdx][0]))+'\t'+str(int(amplfdSeqs[seqIdx][1]))+'\t'+'\n') #write seqIdx, count, distance, and bias...for now
            nxtRnd.close()
        else:
            del(slctdSeqs)
            print("SELEX Round "+str(r+1)+" has started")
            totalSeqNum, uniqSeqNum = utils.seqNumberCounter(amplfdSeqs)
            print("total number of sequences in initial pool = "+str(totalSeqNum))
            print("total number of unique sequences in initial pool = "+str(int(uniqSeqNum)))
            slctdSeqs = S.stochasticBasePairSelection(alphabetSet, seqLength, amplfdSeqs, selectionThreshold, uniqSeqNum, totalSeqNum, samplingSize, outputFileNames, r, stringency)
            print("Selection carried for R"+str(r+1))
            del(amplfdSeqs)
            amplfdSeqs = Amplify.randomPCR_with_ErrorsAndBias_FASTv2(slctdSeqs, seqLength, pcrCycleNum, pcrYield, pcrErrorRate, aptamerSeqs, alphabetSet, distanceMeasure)
            print("Amplification carried for R"+str(r+1))
            outFile = outputFileNames + "_R" + str(r+1)
            nxtRnd = open(outFile, 'w')
            print("writing R"+str(r+1)+" seqs to file")
            for seqIdx in amplfdSeqs:
                seq = Apt.pseudoAptamerGenerator(seqIdx, alphabetSet, seqLength)
                nxtRnd.write(str(seq)+'\t'+str(int(amplfdSeqs[seqIdx][0]))+'\t'+str(int(amplfdSeqs[seqIdx][1]))+'\n') #write idx and count for now
            nxtRnd.close()
    print("SELEX completed")
elif(distanceMeasure == "loop"):
# SELEX simulation based on random aptamer assignment, hamming-based definite selection, and
# non-ideal stochastic amplfication with no bias. 
    for r in range(roundNum):
        if(r==0):
            if(aptamerType == 'DNA'):
                alphabetSet = 'ACGT'
                aptamerSeqs, initialSeqNum = Apt.optimumAptamerGenerator(aptamerNum, alphabetSet, seqLength)
            elif(aptamerType == 'RNA'):
                alphabetSet = 'ACGU'
                aptamerSeqs, initialSeqNum = Apt.optimumAptamerGenerator(aptamerNum, alphabetSet, seqLength)
            else:
                print("Error: Simulation of %.s aptamers not supported" %aptamerType)
                break
            if(aptamerNum == 0):
                aptamerSeqs = aptamerSeq
            print("optimum sequences have been chosen")
            print("SELEX Round 1 has started")
            print("total number of sequences in initial library = "+str(initialSeqNum))
            slctdSeqs = S.stochasticLoopSelection_initial(alphabetSet, seqLength, aptamerSeqs, selectionThreshold, initialSeqNum, samplingSize, outputFileNames, r, stringency)
            print("selection carried out for R0")
            amplfdSeqs = Amplify.randomPCR_with_ErrorsAndBias_FASTv2(slctdSeqs, seqLength, pcrCycleNum, pcrYield, pcrErrorRate, aptamerSeqs, alphabetSet, distanceMeasure)
            print("Amplification carried out for R1")
            outFile = outputFileNames + "_R" + str(r+1)
            nxtRnd = open(outFile, 'w')
            print("writing R1 seqs to file")
            for seqIdx in amplfdSeqs:
                seq = Apt.pseudoAptamerGenerator(seqIdx, alphabetSet, seqLength)
                nxtRnd.write(str(seq)+'\t'+str(int(amplfdSeqs[seqIdx][0]))+'\t'+str(int(amplfdSeqs[seqIdx][1]))+'\t'+'\n') #write seqIdx, count, distance, and bias...for now
            nxtRnd.close()
        else:
            del(slctdSeqs)
            print("SELEX Round "+str(r+1)+" has started")
            totalSeqNum, uniqSeqNum = utils.seqNumberCounter(amplfdSeqs)
            print("total number of sequences in initial pool = "+str(totalSeqNum))
            print("total number of unique sequences in initial pool = "+str(int(uniqSeqNum)))
            slctdSeqs = S.stochasticLoopSelection(alphabetSet, seqLength, amplfdSeqs, selectionThreshold, uniqSeqNum, totalSeqNum, samplingSize, outputFileNames, r, stringency)
            print("Selection carried for R"+str(r+1))
            del(amplfdSeqs)
            amplfdSeqs = Amplify.randomPCR_with_ErrorsAndBias_FASTv2(slctdSeqs, seqLength, pcrCycleNum, pcrYield, pcrErrorRate, aptamerSeqs, alphabetSet, distanceMeasure)
            print("Amplification carried for R"+str(r+1))
            outFile = outputFileNames + "_R" + str(r+1)
            nxtRnd = open(outFile, 'w')
            print("writing R"+str(r+1)+" seqs to file")
            for seqIdx in amplfdSeqs:
                seq = Apt.pseudoAptamerGenerator(seqIdx, alphabetSet, seqLength)
                nxtRnd.write(str(seq)+'\t'+str(int(amplfdSeqs[seqIdx][0]))+'\t'+str(int(amplfdSeqs[seqIdx][1]))+'\n') #write idx and count for now
            nxtRnd.close()
    print("SELEX completed")
else:
    print("Invalid argument for distance measure")
if(post_process==True):
    print("Data post-processing has started...")
    dataAnalysis(seqLength, roundNum, outputFileNames, post_process, distanceMeasure)
    print("Data post-processing is complete")
    print("The Simulation has ended")
else:
    print("The Simulation has ended without post-processing")
