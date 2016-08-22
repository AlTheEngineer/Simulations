import sys, time
import gc
from numpy import random
import numpy as np
import linecache
from itertools import izip, imap, islice, product
import operator
from collections import OrderedDict
from scipy import stats
import Aptamers, Distance, utils
from utils import apt_loopFinder

D = Distance.Distance()
Apt = Aptamers.Aptamers()
#append path for ViennaRNA module
sys.path.append("/local/data/public/aaaa3/Simulations/ViennaRNA/lib/python2.7/site-packages/")
import RNA
from RNA import fold, bp_distance
## NEED TO CHANGE SAMPLING FOR SELECTION TO BE WEIGHTED BY COUNT OF EACH UNIQUE SEQ


class Selection: 

    def selectionProcess_loop_initial(self, slctdSeqs, aptSeq, 
                                      aptStruct, aptLoop,
                                      selectionThreshold, 
                                      alphabetSet, seqLength, 
                                      totalSeqNum, stringency):
        selectedSeqs = 0
        j=0
        while(selectedSeqs < selectionThreshold):
            randIdxs = random.randint(0, int(totalSeqNum-1), size=10000000)
            randHamms = random.randint(0, seqLength-stringency, size=10000000)
            print("sample batch drawn")
            for i, randIdx in enumerate(randIdxs):
                randSeq = Apt.pseudoAptamerGenerator(randIdx, alphabetSet, seqLength)
                randSeqDist = D.loop_func(aptSeq, aptStruct, aptLoop, randSeq, seqLength)
                if(selectedSeqs == selectionThreshold):
                    del(randIdxs)
                    del(randHamms)
                    gc.collect()
                    return slctdSeqs
                elif(randSeqDist < randHamms[i]):
                    if(randIdx in slctdSeqs):
                        slctdSeqs[randIdx][0] += 1
                    else:
                        randSeqBias = D.bias_func(randSeq, seqLength)
                        slctdSeqs[randIdx] = np.array([1, randSeqDist, randSeqBias])
                    selectedSeqs += 1

    def stochasticLoopSelection_initial(self, alphabetSet, seqLength, 
                                            aptPool, selectionThreshold, 
                                            totalSeqNum, samplingSize, 
                                            outputFileNames, rnd, stringency):
        #sampling
        print("sampling from initial library...")
        randomSamples = random.randint(0, int(totalSeqNum-1), size=samplingSize)
        sampleFileName = outputFileNames+"_samples_R"+str(rnd)
        with open(sampleFileName, 'w') as s:
            for seqIdx in randomSamples:
                seq = Apt.pseudoAptamerGenerator(seqIdx, alphabetSet, seqLength)
                s.write(seq+'\n')
        print("Sampling completed")
        #initialize seqInfo matrix
        slctdSeqs = {} 
        selectedSeqs = 0
        aptStruct = fold(aptPool)[0]
        aptLoop = apt_loopFinder(aptPool, aptStruct)
        print("Selection has started")
        #stochastic selection until threshold is met   
        slctdSeqs = self.selectionProcess_loop_initial(slctdSeqs, aptPool, 
                                                       aptStruct, aptLoop, 
                                                       selectionThreshold, 
                                                       alphabetSet, seqLength, 
                                                       totalSeqNum, stringency)
        print("sequence selection has been carried out")
        return slctdSeqs


    def stochasticLoopSelection(self, alphabetSet, seqLength, 
                                    seqPool, selectionThreshold, 
                                    uniqSeqNum, totalSeqNum, samplingSize, 
                                    outputFileNames, rnd, stringency):
        #initialize selected sequence pool
        slctdSeqs = {}
        selectedSeqs = 0
        print("seq length = "+str(seqLength))
        print("seq selection threshold = "+str(selectionThreshold))
        print("unique seq number = "+str(uniqSeqNum))
        print("sample distance = "+str(seqPool[seqPool.keys()[5]][1]))
        print("parameters for selection have been initialized")
        x = np.zeros((uniqSeqNum, 4))
        for i, seqIdx in enumerate(seqPool):
            x[i][0] = seqIdx
            x[i][1] = seqPool[seqIdx][0]
            x[i][2] = seqPool[seqIdx][1]
            x[i][3] = seqPool[seqIdx][2]
        del(seqPool)
        gc.collect()
        print("Selection sample distribution being computed...")
        #compute sampling distribution for selection
        selectionDist = utils.rvd(x, totalSeqNum, "selectionDist")
        print("Selection sample distribution computed")
        print("Sampling has started...")
        randSamples = selectionDist.rvs(size=samplingSize)
        sampleFileName = outputFileNames+"_samples_R"+str(rnd)
        with open(sampleFileName, 'w') as s:
            for seqIdx in randSamples:
                seq = Apt.pseudoAptamerGenerator(seqIdx, alphabetSet, seqLength)
                s.write(str(seq)+'\t'+str(int(x[seqIdx][1]))+'\n')
        print("Sampling has completed")
        for i, seqIdx in enumerate(x):
            x[i][1] = 0
        x = self.selectionProcess(x, selectionThreshold, 
                                  selectionDist, seqLength, 
                                  stringency)
        x = x[x[:, 1] != 0]
        for seqInfo in x:
            #change it so that seqInfo are added as on np array, without setdefault
            slctdSeqs[int(seqInfo[0])] = seqInfo[1:]
        del(x)
        gc.collect()
        print("sequence selection has been carried out")
        return slctdSeqs


    def selectionProcess_2D_initial(self, slctdSeqs, aptStruct, 
                                    selectionThreshold, 
                                    alphabetSet, seqLength, 
                                    totalSeqNum, stringency):
        selectedSeqs = 0
        j=0
        while(selectedSeqs < selectionThreshold):
            randIdxs = random.randint(0, int(totalSeqNum-1), size=10000000)
            randHamms = random.randint(0, seqLength-stringency, size=10000000)
            print("sample batch drawn")
            for i, randIdx in enumerate(randIdxs):
                randSeq = Apt.pseudoAptamerGenerator(randIdx, alphabetSet, seqLength)
                randSeqDist = D.bp_func(aptStruct, randSeq)
                if(selectedSeqs == selectionThreshold):
                    del(randIdxs)
                    del(randHamms)
                    gc.collect()
                    return slctdSeqs
                elif(randSeqDist < randHamms[i]):
                    if(randIdx in slctdSeqs):
                        slctdSeqs[randIdx][0] += 1
                    else:
                        randSeqBias = D.bias_func(randSeq, seqLength)
                        slctdSeqs[randIdx] = np.array([1, randSeqDist, randSeqBias])
                    selectedSeqs += 1


    def stochasticBasePairSelection_initial(self, alphabetSet, seqLength, 
                                            aptPool, selectionThreshold, 
                                            totalSeqNum, samplingSize, 
                                            outputFileNames, rnd, stringency):
        #sampling
        print("sampling from initial library...")
        randomSamples = random.randint(0, int(totalSeqNum-1), size=samplingSize)
        sampleFileName = outputFileNames+"_samples_R"+str(rnd)
        with open(sampleFileName, 'w') as s:
            for seqIdx in randomSamples:
                seq = Apt.pseudoAptamerGenerator(seqIdx, alphabetSet, seqLength)
                s.write(seq+'\n')
        print("Sampling completed")
        #initialize seqInfo matrix
        slctdSeqs = {} 
        selectedSeqs = 0
        aptStruct = fold(aptPool)[0]
        print("Selection has started")
        #stochastic selection until threshold is met   
        slctdSeqs = self.selectionProcess_2D_initial(slctdSeqs, aptStruct, 
                                                     selectionThreshold, 
                                                     alphabetSet, seqLength, 
                                                     totalSeqNum, stringency)
        print("sequence selection has been carried out")
        return slctdSeqs

    def stochasticBasePairSelection(self, alphabetSet, seqLength, 
                                    seqPool, selectionThreshold, 
                                    uniqSeqNum, totalSeqNum, samplingSize, 
                                    outputFileNames, rnd, stringency):
        #initialize selected sequence pool
        slctdSeqs = {}
        selectedSeqs = 0
        print("seq length = "+str(seqLength))
        print("seq selection threshold = "+str(selectionThreshold))
        print("unique seq number = "+str(uniqSeqNum))
        print("sample distance = "+str(seqPool[seqPool.keys()[5]][1]))
        print("parameters for selection have been initialized")
        x = np.zeros((uniqSeqNum, 4))
        for i, seqIdx in enumerate(seqPool):
            x[i][0] = seqIdx
            x[i][1] = seqPool[seqIdx][0]
            x[i][2] = seqPool[seqIdx][1]
            x[i][3] = seqPool[seqIdx][2]
        del(seqPool)
        gc.collect()
        print("Selection sample distribution being computed...")
        #compute sampling distribution for selection
        selectionDist = utils.rvd(x, totalSeqNum, "selectionDist")
        print("Selection sample distribution computed")
        print("Sampling has started...")
        randSamples = selectionDist.rvs(size=samplingSize)
        sampleFileName = outputFileNames+"_samples_R"+str(rnd)
        with open(sampleFileName, 'w') as s:
            for seqIdx in randSamples:
                seq = Apt.pseudoAptamerGenerator(seqIdx, alphabetSet, seqLength)
                s.write(str(seq)+'\t'+str(int(x[seqIdx][1]))+'\n')
        print("Sampling has completed")
        for i, seqIdx in enumerate(x):
            x[i][1] = 0
        x = self.selectionProcess(x, selectionThreshold, 
                                  selectionDist, seqLength, 
                                  stringency)
        x = x[x[:, 1] != 0]
        for seqInfo in x:
            #change it so that seqInfo are added as on np array, without setdefault
            slctdSeqs[int(seqInfo[0])] = seqInfo[1:]
        del(x)
        gc.collect()
        print("sequence selection has been carried out")
        return slctdSeqs

    def selectionProcess_1D_initial(self, slctdSeqs, aptPool, 
                                 selectionThreshold, 
                                 alphabetSet, seqLength, 
                                 totalSeqNum, stringency):
        selectedSeqs = 0
        j=0
        while(selectedSeqs < selectionThreshold):
            randIdxs = random.randint(0, int(totalSeqNum-1), size=10000000)
            randHamms = random.randint(0, seqLength-stringency, size=10000000)
            print("sample batch drawn")
            for i, randIdx in enumerate(randIdxs):
                randSeq = Apt.pseudoAptamerGenerator(randIdx, alphabetSet, seqLength)
                randSeqDist = D.hamming_func(randSeq, aptPool)
                if(selectedSeqs == selectionThreshold):
                    del(randIdxs)
                    del(randHamms)
                    gc.collect()
                    return slctdSeqs
                elif(randSeqDist < randHamms[i]):
                    if(randIdx in slctdSeqs):
                        slctdSeqs[randIdx][0] += 1
                    else:
                        randSeqBias = D.bias_func(randSeq, seqLength)
                        slctdSeqs[randIdx] = np.array([1, randSeqDist, randSeqBias])
                    selectedSeqs += 1



    def stochasticHammingSelection_initial(self, alphabetSet, seqLength, 
                                           aptPool, selectionThreshold, 
                                           totalSeqNum, samplingSize, 
                                           outputFileNames, rnd, stringency):
        #sampling
        print("sampling from initial library...")
        randomSamples = random.randint(0, int(totalSeqNum-1), size=samplingSize)
        sampleFileName = outputFileNames+"_samples_R"+str(rnd)
        with open(sampleFileName, 'w') as s:
            for seqIdx in randomSamples:
                seq = Apt.pseudoAptamerGenerator(seqIdx, alphabetSet, seqLength)
                s.write(seq+'\n')
        print("Sampling completed")
        #initialize seqInfo matrix
        slctdSeqs = {} 
        selectedSeqs = 0
        print("Selection has started...")
        slctdSeqs = self.selectionProcess_1D_initial(slctdSeqs, 
                                                     aptPool, selectionThreshold, 
                                                     alphabetSet, seqLength, 
                                                     totalSeqNum, stringency)
        print("sequence selection has been carried out")
        return slctdSeqs

    def stochasticHammingSelection(self, alphabetSet, seqLength, 
                                   seqPool, selectionThreshold, 
                                   uniqSeqNum, totalSeqNum, samplingSize, 
                                   outputFileNames, rnd, stringency):
        #initialize selected sequence pool
        slctdSeqs = {}
        selectedSeqs = 0
        print("seq length = "+str(seqLength))
        print("seq selection threshold = "+str(selectionThreshold))
        print("unique seq number = "+str(uniqSeqNum))
        print("eample distance = "+str(seqPool[seqPool.keys()[5]][1]))
        print("parameters for selection have been initialized")
        x = np.zeros((uniqSeqNum, 4))
        #transfer selected pool to matrix x
        for i, seqIdx in enumerate(seqPool):
            x[i][0] = seqIdx
            #seq count
            x[i][1] = seqPool[seqIdx][0]
            #seq distance
            x[i][2] = seqPool[seqIdx][1]
            #seq bias
            x[i][3] = seqPool[seqIdx][2]
        del(seqPool)
        gc.collect()
        print("Selection sample distribution being computed...")
        #distribution computed using count of each unique seq
        selectionDist = utils.rvd(x, totalSeqNum, "selectionDist")
        print("Selection sample distribution computed")
        print("Sampling has started...")
        #draw random samples from distribution
        randSamples = selectionDist.rvs(size=samplingSize)
        #write to samples file
        sampleFileName = outputFileNames+"_samples_R"+str(rnd)
        with open(sampleFileName, 'w') as s:
            for seqIdx in randSamples:
                seq = Apt.pseudoAptamerGenerator(seqIdx, alphabetSet, seqLength)
                s.write(str(seq)+'\t'+str(int(x[seqIdx][1]))+'\n')
        print("Sampling has completed")
        #reset all seq counts prior to selection
        for i, seqIdx in enumerate(x):
            x[i][1] = 0
        #draw a bunch of random seqs
        x = self.selectionProcess(x, selectionThreshold, selectionDist, 
                                  seqLength, stringency)
        #remove all seqs that haven't been selected
        x = x[x[:, 1] != 0]
        #transfer info back to selected pool
        for seqInfo in x:
            slctdSeqs[int(seqInfo[0])] = seqInfo[1:]
        del(x)
        gc.collect()
        print("sequence selection has been carried out")
        return slctdSeqs

    def selectionProcess(self, x, selectionThreshold, 
                         selectionDist, seqLength, 
                         stringency):
        selectedSeqs = 0
        while(selectedSeqs < selectionThreshold):
            randIdxs = selectionDist.rvs(size=1000000)
            randHamms = random.randint(0, seqLength-stringency, size=1000000)
            print("sample batch drawn")
            for i, randIdx in enumerate(randIdxs):
                if(selectedSeqs == selectionThreshold):
                    del(randIdxs)
                    del(randHamms)
                    gc.collect()
                    return x
                elif(int(x[randIdx][2]) < randHamms[i]):
                    x[randIdx][1] += 1
                    selectedSeqs += 1

    def definiteSelection(self, seqs, selection_rate, totalseqs):
        seqs_ordrd = OrderedDict(sorted(seqs.items(), key=lambda seqs: seqs[1][1], reverse=False)) #sort seqs by hammdist
        print("seqs have been sorted by distance value in descending order")
   #i=0
        slctd_seqs={}
        sampled_seqs=0
        print("parameters for selection have been initialized")
   # select top 20 percent in terms of hamm distance
        for seq in seqs_ordrd:
            if(sampled_seqs <= (selection_rate*totalseqs)):
       # THIS IS TAKING TOO LONG # UPDATE: NOT ANYMORE
                slctd_seqs[seq] = seqs_ordrd[seq] #add to selected pool
                sampled_seqs += slctd_seqs[seq][0] #add number of repeats
            else:
                break
       #i+=1 #increment index
   #seqs_ordrd.clear() #clear table of previous round
        print("sequence selection has been carried out")
        return slctd_seqs

    def stochasticSelection(self, seqLength, seqs, selection_rate, totalseqs):
        seqs_ordrd = OrderedDict(sorted(seqs.items(), key=lambda seqs: seqs[1][1], reverse=False)) #sort seqs by hammdist
        print("seqs have been sorted by distance value in descending order")
   #i=0
        slctd_seqs={}
        sampled_seqs=0
        print("parameters for selection have been initialized")
   # select top 20 percent in terms of hamm distance
        for seq in seqs_ordrd:
            if(sampled_seqs <= (selection_rate*totalseqs)):
                randHamm = random.randint(0, seqLength)
                if(seqs_ordrd[seq][1] <= randHamm):
       # THIS IS TAKING TOO LONG # UPDATE: NOT ANYMORE
                    slctd_seqs[seq] = seqs_ordrd[seq] #add to selected pool
                    sampled_seqs += slctd_seqs[seq][0] #add number of repeats
                else:
                    continue
            else:
                break
       #i+=1 #increment index
   #seqs_ordrd.clear() #clear table of previous round
        print("sequence selection has been carried out")
        return slctd_seqs

