## NOTE: put all input params and output components in one numpy
#        matrix. Use gather() to collect local objects of the 
#        matrix. 
## This script is used to generate a data set of frequency dists after PCR under different 
#   initial conditions. It then approximates each dist using a Gaussian Mixture model (GMM).
#   A histogram plot of the dist and GMM fitting is generated.
#   It follows the design of trainingSetGenerator.py but normalises
#   the gmm parameters based on expectation (see BruteGMMnormed). It also sorts the 
#   parameters of each gmm according to mean magnitudes



import numpy as np
import Amplification
import sys
from mpi4py import MPI

maxInitialCount = int(sys.argv[1])
maxPCRcycles = int(sys.argv[2])
maxYield = float(sys.argv[3])
yieldInc = float(sys.argv[4]) #yield increment
maxGaussNum = int(sys.argv[5])
outFile = str(sys.argv[6])

rank = MPI.COMM_WORLD.Get_rank()
size = MPI.COMM_WORLD.Get_size()


amp = Amplification.Amplification()

inputDim = 3

sampleNum, outputDim = (maxInitialCount*maxPCRcycles*(maxYield/yieldInc)), (maxGaussNum*3)

dimension = outputDim + inputDim

#initialize input vectors
data = np.zeros((size+1, (sampleNum/size), dimension))

y = float() #initialize yield index
y = 0.0

if rank == 0:
    #initialize local sample index
    j = 0

    for i in range(maxInitialCount):
        for n in range(maxPCRcycles):
            modelParams = amp.BruteGMMnormed(i+1, n+1, 0.1, 10000, maxGaussNum)
            data[rank+1][j][2] = 0.1
            data[rank+1][j][0] = i+1
            data[rank+1][j][1] = n+1
            for k, mu in enumerate(modelParams):
                data[rank+1][j][(k+1)*3] = modelParams[k][0]
                data[rank+1][j][((k+1)*3)+1] = modelParams[k][1]
                data[rank+1][j][((k+1)*3)+2] = modelParams[k][2]
            j+=1 #increment sample index

    subOutfileName = outFile + str(rank)
    subOutfile = open(subOutfileName, 'w')
    writeString = str()


    for l, param in enumerate(data[rank+1]):
        subOutfile.write(str(param[0])+'\t'+str(param[1])+'\t'+str(param[2])+'\t'+str(param[3])+'\t'+str(param[4])+'\t'+str(param[5])+'\t'+str(param[6])+'\t'+str(param[7])+'\t'+str(param[8])+'\t'+str(param[9])+'\t'+str(param[10])+'\t'+str(param[11])+'\t'+str(param[12])+'\t'+str(param[13])+'\t'+str(param[14])+'\t'+'\n')

    subOutfile.close()

if rank == 1:
    #initialize local sample index
    j = 0

    for i in range(maxInitialCount):
        for n in range(maxPCRcycles):
            modelParams = amp.BruteGMMnormed(i+1, n+1, 0.2, 10000, maxGaussNum)
            data[rank+1][j][2] = 0.2
            data[rank+1][j][0] = i+1
            data[rank+1][j][1] = n+1
            for k, mu in enumerate(modelParams):
                data[rank+1][j][(k+1)*3] = modelParams[k][0]
                data[rank+1][j][((k+1)*3)+1] = modelParams[k][1]
                data[rank+1][j][((k+1)*3)+2] = modelParams[k][2]
            j+=1 #increment sample index

    subOutfileName = outFile + str(rank)
    subOutfile = open(subOutfileName, 'w')
    writeString = str()


    for l, param in enumerate(data[rank+1]):
        subOutfile.write(str(param[0])+'\t'+str(param[1])+'\t'+str(param[2])+'\t'+str(param[3])+'\t'+str(param[4])+'\t'+str(param[5])+'\t'+str(param[6])+'\t'+str(param[7])+'\t'+str(param[8])+'\t'+str(param[9])+'\t'+str(param[10])+'\t'+str(param[11])+'\t'+str(param[12])+'\t'+str(param[13])+'\t'+str(param[14])+'\t'+'\n')

    subOutfile.close()

if rank == 2:
    #initialize local sample index
    j = 0

    for i in range(maxInitialCount):
        for n in range(maxPCRcycles):
            modelParams = amp.BruteGMMnormed(i+1, n+1, 0.3, 10000, maxGaussNum)
            data[rank+1][j][2] = 0.3
            data[rank+1][j][0] = i+1
            data[rank+1][j][1] = n+1
            for k, mu in enumerate(modelParams):
                data[rank+1][j][(k+1)*3] = modelParams[k][0]
                data[rank+1][j][((k+1)*3)+1] = modelParams[k][1]
                data[rank+1][j][((k+1)*3)+2] = modelParams[k][2]
            j+=1 #increment sample index

    subOutfileName = outFile + str(rank)
    subOutfile = open(subOutfileName, 'w')
    writeString = str()


    for l, param in enumerate(data[rank+1]):
        subOutfile.write(str(param[0])+'\t'+str(param[1])+'\t'+str(param[2])+'\t'+str(param[3])+'\t'+str(param[4])+'\t'+str(param[5])+'\t'+str(param[6])+'\t'+str(param[7])+'\t'+str(param[8])+'\t'+str(param[9])+'\t'+str(param[10])+'\t'+str(param[11])+'\t'+str(param[12])+'\t'+str(param[13])+'\t'+str(param[14])+'\t'+'\n')

    subOutfile.close()

if rank == 3:
    #initialize local sample index
    j = 0

    for i in range(maxInitialCount):
        for n in range(maxPCRcycles):
            modelParams = amp.BruteGMMnormed(i+1, n+1, 0.4, 10000, maxGaussNum)
            data[rank+1][j][2] = 0.4
            data[rank+1][j][0] = i+1
            data[rank+1][j][1] = n+1
            for k, mu in enumerate(modelParams):
                data[rank+1][j][(k+1)*3] = modelParams[k][0]
                data[rank+1][j][((k+1)*3)+1] = modelParams[k][1]
                data[rank+1][j][((k+1)*3)+2] = modelParams[k][2]
            j+=1 #increment sample index

    subOutfileName = outFile + str(rank)
    subOutfile = open(subOutfileName, 'w')
    writeString = str()


    for l, param in enumerate(data[rank+1]):
        subOutfile.write(str(param[0])+'\t'+str(param[1])+'\t'+str(param[2])+'\t'+str(param[3])+'\t'+str(param[4])+'\t'+str(param[5])+'\t'+str(param[6])+'\t'+str(param[7])+'\t'+str(param[8])+'\t'+str(param[9])+'\t'+str(param[10])+'\t'+str(param[11])+'\t'+str(param[12])+'\t'+str(param[13])+'\t'+str(param[14])+'\t'+'\n')

    subOutfile.close()

if rank == 4:
    #initialize local sample index
    j = 0

    for i in range(maxInitialCount):
        for n in range(maxPCRcycles):
            modelParams = amp.BruteGMMnormed(i+1, n+1, 0.5, 10000, maxGaussNum)
            data[rank+1][j][2] = 0.5
            data[rank+1][j][0] = i+1
            data[rank+1][j][1] = n+1
            for k, mu in enumerate(modelParams):
                data[rank+1][j][(k+1)*3] = modelParams[k][0]
                data[rank+1][j][((k+1)*3)+1] = modelParams[k][1]
                data[rank+1][j][((k+1)*3)+2] = modelParams[k][2]
            j+=1 #increment sample index

    subOutfileName = outFile + str(rank)
    subOutfile = open(subOutfileName, 'w')
    writeString = str()


    for l, param in enumerate(data[rank+1]):
        subOutfile.write(str(param[0])+'\t'+str(param[1])+'\t'+str(param[2])+'\t'+str(param[3])+'\t'+str(param[4])+'\t'+str(param[5])+'\t'+str(param[6])+'\t'+str(param[7])+'\t'+str(param[8])+'\t'+str(param[9])+'\t'+str(param[10])+'\t'+str(param[11])+'\t'+str(param[12])+'\t'+str(param[13])+'\t'+str(param[14])+'\t'+'\n')

    subOutfile.close()
if rank == 5:
    #initialize local sample index
    j = 0

    for i in range(maxInitialCount):
        for n in range(maxPCRcycles):
            modelParams = amp.BruteGMMnormed(i+1, n+1, 0.6, 10000, maxGaussNum)
            data[rank+1][j][2] = 0.6
            data[rank+1][j][0] = i+1
            data[rank+1][j][1] = n+1
            for k, mu in enumerate(modelParams):
                data[rank+1][j][(k+1)*3] = modelParams[k][0]
                data[rank+1][j][((k+1)*3)+1] = modelParams[k][1]
                data[rank+1][j][((k+1)*3)+2] = modelParams[k][2]
            j+=1 #increment sample index

    subOutfileName = outFile + str(rank)
    subOutfile = open(subOutfileName, 'w')
    writeString = str()


    for l, param in enumerate(data[rank+1]):
        subOutfile.write(str(param[0])+'\t'+str(param[1])+'\t'+str(param[2])+'\t'+str(param[3])+'\t'+str(param[4])+'\t'+str(param[5])+'\t'+str(param[6])+'\t'+str(param[7])+'\t'+str(param[8])+'\t'+str(param[9])+'\t'+str(param[10])+'\t'+str(param[11])+'\t'+str(param[12])+'\t'+str(param[13])+'\t'+str(param[14])+'\t'+'\n')

    subOutfile.close()
if rank == 6:
    #initialize local sample index
    j = 0

    for i in range(maxInitialCount):
        for n in range(maxPCRcycles):
            modelParams = amp.BruteGMMnormed(i+1, n+1, 0.7, 10000, maxGaussNum)
            data[rank+1][j][2] = 0.7
            data[rank+1][j][0] = i+1
            data[rank+1][j][1] = n+1
            for k, mu in enumerate(modelParams):
                data[rank+1][j][(k+1)*3] = modelParams[k][0]
                data[rank+1][j][((k+1)*3)+1] = modelParams[k][1]
                data[rank+1][j][((k+1)*3)+2] = modelParams[k][2]
            j+=1 #increment sample index

    subOutfileName = outFile + str(rank)
    subOutfile = open(subOutfileName, 'w')
    writeString = str()


    for l, param in enumerate(data[rank+1]):
        subOutfile.write(str(param[0])+'\t'+str(param[1])+'\t'+str(param[2])+'\t'+str(param[3])+'\t'+str(param[4])+'\t'+str(param[5])+'\t'+str(param[6])+'\t'+str(param[7])+'\t'+str(param[8])+'\t'+str(param[9])+'\t'+str(param[10])+'\t'+str(param[11])+'\t'+str(param[12])+'\t'+str(param[13])+'\t'+str(param[14])+'\t'+'\n')

    subOutfile.close()
if rank == 7:
    #initialize local sample index
    j = 0

    for i in range(maxInitialCount):
        for n in range(maxPCRcycles):
            modelParams = amp.BruteGMMnormed(i+1, n+1, 0.8, 10000, maxGaussNum)
            data[rank+1][j][2] = 0.8
            data[rank+1][j][0] = i+1
            data[rank+1][j][1] = n+1
            for k, mu in enumerate(modelParams):
                data[rank+1][j][(k+1)*3] = modelParams[k][0]
                data[rank+1][j][((k+1)*3)+1] = modelParams[k][1]
                data[rank+1][j][((k+1)*3)+2] = modelParams[k][2]
            j+=1 #increment sample index

    subOutfileName = outFile + str(rank)
    subOutfile = open(subOutfileName, 'w')
    writeString = str()

    for l, param in enumerate(data[rank+1]):
        subOutfile.write(str(param[0])+'\t'+str(param[1])+'\t'+str(param[2])+'\t'+str(param[3])+'\t'+str(param[4])+'\t'+str(param[5])+'\t'+str(param[6])+'\t'+str(param[7])+'\t'+str(param[8])+'\t'+str(param[9])+'\t'+str(param[10])+'\t'+str(param[11])+'\t'+str(param[12])+'\t'+str(param[13])+'\t'+str(param[14])+'\t'+'\n')

    subOutfile.close()
if rank == 8:
    #initialize local sample index
    j = 0

    for i in range(maxInitialCount):
        for n in range(maxPCRcycles):
            modelParams = amp.BruteGMMnormed(i+1, n+1, 0.9, 10000, maxGaussNum)
            data[rank+1][j][2] = 0.9
            data[rank+1][j][0] = i+1
            data[rank+1][j][1] = n+1
            for k, mu in enumerate(modelParams):
                data[rank+1][j][(k+1)*3] = modelParams[k][0]
                data[rank+1][j][((k+1)*3)+1] = modelParams[k][1]
                data[rank+1][j][((k+1)*3)+2] = modelParams[k][2]
            j+=1 #increment sample index

    subOutfileName = outFile + str(rank)
    subOutfile = open(subOutfileName, 'w')
    writeString = str()

    for l, param in enumerate(data[rank+1]):
        subOutfile.write(str(param[0])+'\t'+str(param[1])+'\t'+str(param[2])+'\t'+str(param[3])+'\t'+str(param[4])+'\t'+str(param[5])+'\t'+str(param[6])+'\t'+str(param[7])+'\t'+str(param[8])+'\t'+str(param[9])+'\t'+str(param[10])+'\t'+str(param[11])+'\t'+str(param[12])+'\t'+str(param[13])+'\t'+str(param[14])+'\t'+'\n')

    subOutfile.close()
if rank == 9:
    #initialize local sample index
    j = 0

    for i in range(maxInitialCount):
        for n in range(maxPCRcycles):
            modelParams = amp.BruteGMMnormed(i+1, n+1, 0.95, 10000, maxGaussNum)
            data[rank+1][j][2] = 0.95
            data[rank+1][j][0] = i+1
            data[rank+1][j][1] = n+1
            for k, mu in enumerate(modelParams):
                data[rank+1][j][(k+1)*3] = modelParams[k][0]
                data[rank+1][j][((k+1)*3)+1] = modelParams[k][1]
                data[rank+1][j][((k+1)*3)+2] = modelParams[k][2]
            j+=1 #increment sample index

    subOutfileName = outFile + str(rank)
    subOutfile = open(subOutfileName, 'w')
    writeString = str()


    for l, param in enumerate(data[rank+1]):
        subOutfile.write(str(param[0])+'\t'+str(param[1])+'\t'+str(param[2])+'\t'+str(param[3])+'\t'+str(param[4])+'\t'+str(param[5])+'\t'+str(param[6])+'\t'+str(param[7])+'\t'+str(param[8])+'\t'+str(param[9])+'\t'+str(param[10])+'\t'+str(param[11])+'\t'+str(param[12])+'\t'+str(param[13])+'\t'+str(param[14])+'\t'+'\n')

    subOutfile.close()
