import numpy as np
import matplotlib.pyplot as plt
import csv
import pandas as pd

inputfile = "results/attack_success.csv"
df = pd.read_csv(inputfile,header=None)

QBERrange = np.arange(0.01,0.12,0.01).round(2)
bitRange = np.arange(50,1000,20)


for QBER in QBERrange:
    QBERdf = df[df[1] == QBER]
    iter1 = []
    iter2 = []
    iter3 = []
    for N in bitRange:
        Ndf = QBERdf[QBERdf[0] == N]
        for rep in range(3):
            repdf = Ndf[Ndf[2] == rep]
            success_bits = repdf[4].to_numpy()
            success_rate = np.average(success_bits)/N
            if rep == 0:
                iter1.append(success_rate)
            if rep == 1:
                iter2.append(success_rate)
            if rep == 2:
                iter3.append(success_rate) 

    plt.figure()
   # plt.ylim([0,1])
    plt.plot(bitRange,iter1,label="Iteration 1")
    plt.plot(bitRange,iter2,label="Iteration 2")
    plt.plot(bitRange,iter3,label="Iteration 3")
    plt.xlabel("Key Length")
    plt.ylabel("Key recovery rate")
    plt.title("Key recovery for QBER: {0}".format(QBER))
    plt.legend()
    
    figname = "figures/QBER_{0}_iterations.png".format(QBER)
    plt.savefig(figname)


for N in bitRange:
    Ndf = df[df[0] == N]
    iter1 = []
    iter2 = []
    iter3 = []
    for QBER in QBERrange:
        QBERdf = Ndf[Ndf[1] == QBER]
        for rep in range(3):
            repdf = QBERdf[QBERdf[2] == rep]
            success_bits = repdf[4].to_numpy()
            success_rate = np.average(success_bits)/N
            if rep == 0:
                iter1.append(success_rate)
            if rep == 1:
                iter2.append(success_rate)
            if rep == 2:
                iter3.append(success_rate) 

    plt.figure()
   # plt.ylim([0,1])
    plt.plot(QBERrange,iter1,label="Iteration 1")
    plt.plot(QBERrange,iter2,label="Iteration 2")
    plt.plot(QBERrange,iter3,label="Iteration 3")
    plt.xlabel("QBER")
    plt.ylabel("Key recovery rate")
    plt.title("Key recovery for Key Length of {0} bits".format(N))
    plt.legend()
    figname = "figures/N_{0}_iterations.png".format(N)
    plt.savefig(figname)



            
            
            

