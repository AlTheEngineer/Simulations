import time
import random
import linecache
import math
from math import factorial
from itertools import izip, imap
import operator
from collections import OrderedDict
from scipy import stats
import numpy as np
from scipy.misc import comb
import matplotlib
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import Distance


def seqNumberCounter(seqPool):
    totalSeqNum = int(0)
    uniqSeqNum = int(0)
    for seqIdx in seqPool:
        totalSeqNum += seqPool[seqIdx][0]
        uniqSeqNum += 1
    return int(totalSeqNum), int(uniqSeqNum)

def binomCoeff(n, k):
    binom = factorial(n)/(factorial(k)*factorial(n-k))
    return binom


def convert_to_distribution(x, y, distName):
    xDist = stats.rv_discrete(name=distName, values=(x, y))
    return xDist

def apt_loopFinder(apt_seq, apt_struct):
    base = None
    baseIdx = 0
    while(base != ')' and baseIdx < seqLength):
        base = apt_struct[baseIdx]
        baseIdx += 1
    if(baseIdx == seqLength):
        while(base != '('and baseIdx > 1):
            baseIdx -= 1
            base = apt_struct[baseIdx-1]
        if(baseIdx == 1):
            apt_loop = apt_seq
            return apt_loop
        else:
            apt_loop = apt_seq[baseIdx:]
            return apt_loop
    else:
        loop_end = baseIdx-1
        while(base != '(' and baseIdx > 1):
            baseIdx -= 1
            base = apt_struct[baseIdx-1]
        if(baseIdx == 1):
            apt_loop = apt_seq[:loop_end]
            return apt_loop
        else:
            apt_loop = apt_seq[baseIdx:loop_end]
            return apt_loop

# Add method for computing the L1 norm 

def rvd(X, X_sum, distName):
    seqIdxs = np.zeros(X.shape[0])
    probs = np.zeros(X.shape[0])
    for i, seq in enumerate(X):
        seqIdxs[i] = i
        probs[i] = seq[1]/X_sum
    dist = stats.rv_discrete(name=distName, values=(seqIdxs, probs))
    return dist

def bias_avg(seqFile, seqLength):
    bias = 0
    w_bias = 0
    totalSeqs = 0
    uniqSeqs = 0
    with open(seqFile, 'r') as f:
        for line in f:
            row = line.split()
            seq = row[0]
            bias += d.bias_func(seq, seqLength)
            w_bias += int(row[1])*d.bias_func(seq, seqLength)
            totalSeqs += int(row[1])
            uniqSeqs += 1
    avg_bias = bias/uniqSeqs
    weighted_avg_bias = w_bias/totalSeqs
    return weighted_avg_bias, avg_bias

#bias_avg("window_R14", 20)
def seq_div_hamm(seqLength, alphabetSet):
    uniqSeqNum_per_dist = np.zeros(seqLength+1)
    for h in xrange(seqLength+1):
        uniqSeqNum_per_dist[h] = (len(alphabetSet)-1)**(h)*comb(seqLength, h)
    hammDistAxis = np.linspace(0, seqLength, seqLength+1)
    hammDistAxis_smooth = np.linspace(0, seqLength, 200)
    uniqSeqNum_smooth = spline(hammDistAxis, uniqSeqNum_per_dist, hammDistAxis_smooth)
    fig, ax = plt.subplots(1,1)
    ax.plot(hammDistAxis_smooth, uniqSeqNum_smooth)
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    fig.text(0.5, 0.95, 'Unique Sequences', ha='center')
    fig.text(0.5, 0.04, 'Hamming Distance', ha='center')
    fig.text(0.04, 0.5, 'Frequency', va='center', rotation='vertical')
    fig.savefig("SELEX_Analytics_seqDiv_20nt", format='pdf')
    return 0
#seq_div_hamm(20, 'ACGT')

def distance_range(scale, ref_seq, seqLength, alphabetSet):
    ref_struct = fold(ref_seq)[0]
    ref_loop = apt_loopFinder(ref_seq, ref_struct)
    hamm_dist_array = np.zeros(int(seqLength*1.5))
    bp_dist_array = np.zeros(int(seqLength*1.5))
    loop_dist_array = np.zeros(int(seqLength*1.5))
    randIdxs = random.randint(0, 4**(20)-1, size=scale)
    for i in xrange(scale):
        randIdx = randIdxs[i]
        randSeq = apt.pseudoAptamerGenerator(randIdx, alphabetSet, seqLength)
        randHammDist = d.hamming_func(randSeq, ref_seq)
        randbpDist = d.bp_func(ref_struct, randSeq)
        randLoopDist = d.loop_func(ref_seq, ref_struct, ref_loop, randSeq, seqLength)
        hamm_dist_array[randHammDist] += 1
        bp_dist_array[randbpDist] += 1
        loop_dist_array[randLoopDist] += 1
    for dist in xrange(int(seqLength*1.5)):
        hamm_dist_array[dist] /= scale
        bp_dist_array[dist] /= scale
        loop_dist_array[dist] /= scale
    fig, axis = plt.subplots(1,1)
    distAxis = np.linspace(0, int(seqLength+9), int(seqLength+10))
    distAxis_smooth = np.linspace(0, int(seqLength+9), 200)
    hamm_dist_smooth = spline(distAxis, hamm_dist_array, distAxis_smooth)
    bp_dist_smooth = spline(distAxis, bp_dist_array, distAxis_smooth)
    loop_dist_smooth = spline(distAxis, loop_dist_array, distAxis_smooth)
    axis.plot(distAxis_smooth, hamm_dist_smooth, label='Hamming')
    axis.plot(distAxis_smooth, bp_dist_smooth, label='Base-Pair')
    axis.plot(distAxis_smooth, loop_dist_smooth, label='Loop')
    axis.set_xlim([0, 25])
    axis.set_ylim([0, 0.4])
    axis.legend()
    fig.text(0.5, 0.04, 'Distance', ha='center')
    fig.text(0.04, 0.5, 'Fractional Frequency', va='center', rotation='vertical')
    fig.text(0.5, 0.95, 'Distance Distributions', ha='center')
    fig.savefig("SELEX_Analytics_distance_distributions", format='pdf')
    return hamm_dist_array



#distance_range(10000, "GTACGACAGTCATCCTACAC", 20, 'ACGT')

