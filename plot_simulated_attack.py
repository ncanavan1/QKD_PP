import numpy as np
import matplotlib.pyplot as plt
import csv
import pandas as pd

inputfile = "attack_success.csv"
df = pd.read_csv(inputfile,header=None)

QBERrange = np.arange(0.01,0.12,0.01)
bitRange = np.arange(50,201,5)


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
    plt.plot(bitRange,iter1)
    plt.plot(bitRange,iter2)
    plt.plot(bitRange,iter3)
    figname = "figures/QBER_{0}_iterations.png".format(QBER)
    plt.savefig(figname)

            
            
            

