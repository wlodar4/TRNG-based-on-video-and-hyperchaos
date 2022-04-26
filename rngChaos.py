from pickletools import uint8
import wave
import numpy as np
import matplotlib.pyplot as plt
import math
import struct

def getSamples(file):
    nframes=file.getnframes()
    frame_data = file.readframes(nframes)
    if frame_data:
        sample_width = file.getsampwidth()
        nb_samples = len(frame_data) // sample_width
        format = {1:"%db", 2:"<%dh", 4:"<%dl"}[sample_width] % nb_samples
        return np.uint8(list(struct.unpack(format, frame_data)))
    else:
        return ()

def get3LSBS(data):
    data_copy=data
    for i in range (0, len(data_copy)):
        data_copy[i] = data_copy[i] & 7
    return data_copy
    
def saveArrayToFile(data, file):
    for i in range(0,len(data)):
        file.write("{0:b}".format(data[i]))
    file.close()

def entropyAndHistogram(file):
    data=file.read()
    n = 8
    chunks = [data[i:i+n] for i in range(0, len(data), n)]
    values = [int(chunks[i],2) for i in range(0, len(chunks))]
    plt.hist(values, bins=256, range=[0,256])
    plt.show()
    counter = np.zeros(256)
    total = 0
    for i in values:
        counter[i] += 1
        total += 1
    E = 0
    for i in counter:
        if(i != 0):
            E += (i/total) * math.log2(1 / (i/total))
    return E

def tentMap(x):
    alpha = 1.99999
    if (x<0):
        return 0
    else: 
        if (x<0.5):
            return (alpha*x)
        else:
            if (x<1):
                return (alpha*(1-x))
            else:
                return 0

def float2bin(f):
    [d] = struct.unpack(">Q", struct.pack(">d", f))
    return f'{d:064b}'

def swapBits(number):
    x=int(float2bin(number),2)
    if (x==0):
        return 0
    else:
        binaryNum = str(bin(x))
        print (binaryNum)
        n=(int(len(binaryNum)/2))
        p1 = n-1
        p2 = 0
        set1 =  (x >> p1) & ((1<< n) - 1)
        set2 =  (x >> p2) & ((1 << n) - 1)
        xor = (set1 ^ set2)
        xor = (xor << p1) | (xor << p2)
        result = x ^ xor
        return result


def postProcessing(length, plik):
    data = plik
    gamma = 4
    l = 8
    s = 64
    epsilon = 0.05
    n = (2*length)/s
    omega = 0.5
    b = 3
    arrX = [[0.141592, 0.653589, 0.793238, 0.462643, 0.383279, 0.502884, 0.197169, 0.399375], [0,0,0,0,0,0,0,0], [0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0]]
    counter = 0
    z=[0,0,0,0,0,0,0,0,0]
    processedData = ""
    while (len(processedData)< length):
        for i in range (0,(l-1)):
            arrX[0][i] = (((omega*data[counter])/((2**b)-1))*(1/(1+omega)))
            counter=counter+1
        for t in range (0,(gamma - 1)):
            for i in range (0,(l-1)):
                arrX [t+1][i] = (((1-epsilon)*tentMap(arrX[t][i]))+((epsilon/2)*tentMap(arrX[t][(i+1)%l]))+(tentMap(arrX[t][(i-1)%l])))
        for i in range (0,(l-1)):
            z[i]=arrX[gamma-1][i]
            arrX[0][i]=arrX[gamma-1][i]
        for i in range (0,(int((l/2)-1))):
            z[i] =  (z[i] + swapBits(z[int(i+(l/2))]) % 64)
        for i in range (0,(int((l/2)-1))):
            processedData += str(float2bin(z[i]))

    return processedData






plik = wave.open("q.wav")
probki=getSamples(plik) 
outfile = open("data.txt", 'w')
outfile.write(str(probki))
outfile.close()
trzy_bity=get3LSBS(probki)
outfile = open("dane.txt","w")
saveArrayToFile(trzy_bity, outfile)
outfile = open("dane.txt", "r")
entropia = entropyAndHistogram(outfile)
print(entropia)
processed=postProcessing(1000000,probki)
outfile2 = open("dane2.txt", "w")
outfile2.write(str(processed))
outfile2 = open("dane2.txt", "r")
entropia2 = entropyAndHistogram(outfile2)
print(entropia2)
