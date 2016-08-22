import sys, time
import random
import linecache
from itertools import izip, imap, islice
import operator
from collections import OrderedDict
import numpy as np
from scipy import stats
import utils
from utils import apt_loopFinder
#append path for ViennaRNA module
sys.path.append("/local/data/public/aaaa3/Simulations/ViennaRNA/lib/python2.7/site-packages/")
import RNA
from RNA import fold, bp_distance

class Distance:
    
    def lavenshtein_func(self , loop1, loop2):
        if len(loop1) < len(loop2):
            return self.lavenshtein_func(loop2, loop1)
        if len(loop2) == 0:
            return len(loop1)
        loop1 = np.array(tuple(loop1))
        loop2 = np.array(tuple(loop2))
        prev_row = np.arange(loop2.size + 1)
        for nt in loop1:
            curr_row = prev_row +1
            curr_row[1:] = np.minimum(curr_row[1:], np.add(prev_row[:-1], loop2 != nt))
            curr_row[1:] = np.minimum(curr_row[1:], curr_row[0:-1] + 1)
            prev_row = curr_row
        loop2_dist = prev_row[-1]
        return loop2_dist


    def hamming_func(self, seq1, seq2):
        len(seq1) == len(seq2)
        ne = operator.ne
        seq2_dist = sum(imap(ne, seq1, seq2))
        return seq2_dist
    
    def bp_func(self, seq1_struct, seq2):
        seq2_struct = fold(seq2)[0]
        seq2_dist = bp_distance(seq1_struct, seq2_struct)
        return seq2_dist

    def loop_func(self, seq1, seq1_struct, seq1_loop, seq2, seqLength):
        seq2_struct = fold(seq2)[0]
        base = None
        baseIdx = 0
        while(base != ')' and baseIdx < seqLength-1):
            base = seq2_struct[baseIdx]
            baseIdx += 1
        if(baseIdx == seqLength-1):
            while(base != '(' and baseIdx > 0):
                base = seq2_struct[baseIdx-1]
                baseIdx -= 1
            if(baseIdx == 0):
                seq2_loop = seq2
            else:
                seq2_loop = seq2[baseIdx:]
        else:
            loop_end = baseIdx-1
            while(base != '('):
                baseIdx -= 1
                base = seq2_struct[baseIdx-1]
            seq2_loop = seq2[baseIdx:loop_end]
        seq2_loopDist = self.lavenshtein_func(seq1_loop, seq2_loop)
        seq2_bpDist = bp_distance(seq1_struct, seq2_struct)
        seq2_dist = seq2_loopDist + seq2_bpDist
        return seq2_dist
#TEST
#d = Distance()
#seq_dist = d.loop_func(he4_seq, he4_struct, ex_seq)

    def bias_func(self, seq, seqLen):
        pyrNum = 0
        for nt in seq[:-1]:
            if(nt == 'C') or (nt == 'T'):
                pyrNum += 1 #increment no. of pyrimidines
        biasScore = 0.1*(2*pyrNum - seqLen)/seqLen #compute bias
        return biasScore

    def biasedHamming_initLib(self, seqLen, optimseqs, pool_file, scale, partition):
        seqs=np.zeros((scale, 3))
        seqIdx = 0
        with open(pool_file) as p:
            while True:
                next_n_seqs = list(islice(p, partition))
                if not next_n_seqs:
                    break
                for i, seq in enumerate(next_n_seqs):
                    seqs[seqIdx][0] += 1 #increment count
                    seqs[seqIdx][1] = hamming_func(optimseqs, seq[:-1])
                    for nt in seq[:-1]:
                        if(nt == 'C') or (nt == 'T'):
                            seqs[seqIdx][2] += 1 #increment no. of pyrimidines
                    seqs[seqIdx][2] = 0.1*(2*seqs[seqIdx][2] - seqLen)/seqLen #compute bias
                    seqIdx += 1
        p.close()

        return seqs

##TEST AREA

#d = Distance()
#seqs = d.biasedHamming_initLib(20, 'AAAAAAAAAAAAAAAAAAAA', 'random_initLib_1KM', 100000000, 1000000)


# NOTE: DISTANCES ARE NOW APPENDED AS A 3RD VALUE, NOT 2ND DUE TO BIAS
# This method computes the hamming distance of each sequence to the optimum set
# and their PCR bias. It requires sequence length as 3rd argument
    def seqsHamming_and_Bias(self, optimseqs, pool_file, seqLen):
        seqs={}
        seqsfile = open(pool_file)
        seq=seqsfile.readline()
        while(seq):
            if seq not in seqs:
                seqs.setdefault(seq, []).append(1) #append seq count
                seqs.setdefault(seq, []).append(0) #append seq bias
                for nt in seq:
                    if(nt == 'C') or (nt == 'T'):
                        seqs[seq][1]+=1 #increment no. of pyrimidines
                seqs[seq][1] = 0.1*(2*seqs[seq][1] - seqLen)/seqLen #compute bias
                for optimseq in optimseqs:
                    hammdist = hamming_func(optimseq, seq) #compute distance
                    seqs.setdefault(seq, []).append(hammdist) #append distance
            else:
                seqs[seq][0]+=1 #increment count
            seq=seqsfile.readline()
        seqsfile.close()
        print("distance calculations and bias scoring completed")
        return seqs

    def seqsHamming(self, optimseqs, pool_file):
        seqs=np.zeros((scale, 3))
        seqsfile = open(pool_file)
        seq=seqsfile.readline()

        while(seq):
            if seq not in seqs:
                seqs.setdefault(seq, []).append(1)
            
                for optimseq in optimseqs:
                    hammdist = hamming_func(optimseq, seq)
                    seqs.setdefault(seq, []).append(hammdist)
            else:
                seqs[seq][0]+=1
            seq=seqsfile.readline()
        seqsfile.close()
        print("distance calculations passed")
        return seqs

