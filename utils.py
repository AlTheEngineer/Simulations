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
    while(base != ')'):
        base = apt_struct[baseIdx]
        baseIdx += 1
    loop_end = baseIdx-1
    while(base != '('):
        baseIdx -= 1
        base = apt_struct[baseIdx-1]
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
        uniqSeqNum_per_dist[h] = (len(alphabetSet)-1)*comb(seqLength, h)
    hammDistAxis = np.linspace(0, seqLength+1, seqLength+1)
    fig, ax = plt.subplots(1,1)
    ax.plot(hammDistAxis, uniqSeqNum_per_dist)
    ax.grid()
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    fig.text(0.5, 0.95, 'Unique Sequences', ha='center')
    fig.text(0.5, 0.04, 'Hamming Distance', ha='center')
    fig.text(0.04, 0.5, 'Frequency', va='center', rotation='vertical')
    fig.savefig("SELEX_Analytics_seqDiv_20nt", format='pdf')
    return 0
#seq_div_hamm(20, 'ACTG')
