import numpy as np
import random


def gen_sifted_keys(N,QBER):
    xA = np.zeros([1,N])
    xB = np.zeros([1,N])
    N_errs = int(np.floor(N/(100*QBER)))
    err_pos = random.sample(range(0,N),N_errs)
    for i in range(0,N):
        xA[0,i] = np.random.choice([0,1])
        xB[0,i] = xA[0,i]
        if i in err_pos:
            xB[0,i] = (xB[0,i] + 1) % 2
    return xA, xB

##calculates parity of block k
##outputs parity and also the estimates power trace of a register containing the current sum of the parity check
def parity_check(k):
    hw = []
    parity = 0
    for i in range(k.shape[1] - 1):
        if i == 0:
            parity = int(k[0,0]) ^ int(k[0,1])
        else:
            parity = parity ^ int(k[0,i+1])

        if parity == 0:
            hw.append(0)
        else:
            hw.append(1)
    return parity, hw



QBER = 0.1
N = 10
X,Y = gen_sifted_keys(N,QBER)
print(X)

parity, hw = parity_check(X)
print(parity)
print(hw)
