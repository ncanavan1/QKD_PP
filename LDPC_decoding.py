import numpy as np
import random
from pyldpc import make_ldpc, decode


###reference of which checksum(new checksum per row) connects to which message bits
def compact_row(H, d_c):
    rows, cols = H.shape
    R = np.zeros([rows,d_c])
    for k in range(0, rows):
        pos = 0
        for j in range(0, cols):
            if H[k,j] == 1:
                R[k,pos] = j
                pos = pos + 1
    return R

##Position of R_ref is simply an increasing by one refernce from 0
def set_R_ref(H, d_c):
    rows = H.shape[0]
    R_ref = np.zeros([rows, d_c])
    c = 0
    for k in range(0, rows):
        for j in range(0, d_c):
            R_ref[k,j] = R_ref[k,j] + c
            c = c + 1
    return R_ref    

def set_r_ref(H, d_v):
    rows, cols = H.shape
    r_ref = np.zeros([cols, d_v])
    c = 0
    col_pos = np.zeros([cols])
    for i in range(0, rows):
        for j in range(0, cols):
            if H[i,j] == 1:
                r_ref[j,int(col_pos[j])] = c
                c = c + 1
                col_pos[j] = col_pos[j] + 1
    return r_ref



###reference of which message bits(new message bit per row) connects to which check sum bits
def compact_col(H, d_v):
    rows, cols = H.shape
    r = np.zeros([cols,d_v])
    for i in range(0, cols):
        pos = 0
        for j in range(0, rows):
            if H[j,i] == 1:
                r[i,pos] = j
                pos = pos + 1
    return r


def get_CheckSum(H,x):
    return np.matmul(H,x.transpose())%2

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

def get_sign_vector(C,D):
    m = C.shape
    v = np.zeros(m)
    for k in range(m[0]):
        if C[k] == D[k]:
            v[k] = 1
        else:
            v[k] = -1
    return v

def get_t(H):
    rows, cols = H.shape
    t = np.zeros([rows])
    for k in range(0,rows):
        for j in range(0, cols):
            if H[k,j] == 1:
                t[k] = t[k] + 1
    return t

def get_m(H):
    rows, cols = H.shape
    m = np.zeros([cols])
    for i in range(0,cols):
        for j in range(0, rows):
            if H[j,i] == 1:
                m[i] = m[i] + 1
    return m




def initialise_BP(QBER,C,D,L):
    L = L*(1-2*QBER)
 #   v = get_sign_vector(C,D) ###Done in verify step
    lam = (1-QBER)/QBER
    return L, lam


def step1(L,R_ref,t,v):
    M = R_ref.shape[0]
    L_update = np.zeros(L.shape)
    for k in range(M):
        for j in range(int(t[k])):
            product = v[k,0]
            for i in range(int(t[k])): ###unsure on this line
                if i != j:
                    product = product*L[int(R_ref[k,i])]
            L_update[int(R_ref[k,j])] = product
    return L_update


def step2(L,r_ref,m,lam,Y):
    ##convert from a to f value
    N = r_ref.shape[0]
    for i in range(L.shape[0]):
        L[i] = (1+L[i])/(1-L[i])
    
    Y_dash = Y.copy()
    for i in range(N):
        F_i = lam
        for j in range(int(m[i])):
            F_i = F_i * L[int(r_ref[i,j])]
        if F_i < 1:
            Y_dash[0,i] = (Y[0,i] - 1) % 2
    
    L_update = np.zeros(L.shape)
    for i in range(N):
        for j in range(int(m[i])):
            product = lam
            for k in range(int(m[i])):
                if k != j:
                    product = product * L[int(r_ref[i,k])]
            L_update[int(r_ref[i,j])] = product  
    L = L_update
    for i in range(L.shape[0]):
        L[i] = (L[i] - 1)/(L[i] + 1)

    return Y_dash, L

def verifyY(Y,C,H):
    D = get_CheckSum(H,Y)
    v = get_sign_vector(C,D)
    if (D == C).all():
        return True, v
    else:
        return False, v


def gen_belief_reg(H):
    ones = 0
    for i in range(0,H.shape[0]):
        for j in range(0,H.shape[1]):
            if H[i,j] == 1:
                ones = ones + 1
    return np.ones([ones])


def decode_my(Y,H,C,D,max_iter,QBER,X):

    L = gen_belief_reg(H)
    R_ref = set_R_ref(H, d_c)
    r_ref = set_r_ref(H, d_v)
    t = get_t(H)
    m = get_m(H)
    L, lam = initialise_BP(QBER, C, D, L)
    matching, v = verifyY(Y,C,H)
    iter = 0
    Y_dash = -1
    while matching == False and iter < max_iter:
        L = step1(L,R_ref,t,v)
        Y_dash, L = step2(L,r_ref,m,lam,Y)
        diff = get_difference(X,Y_dash)
        print("{0} Difference in bits: {1}".format(Y_dash, diff))
        matching, v = verifyY(Y_dash,C,H)
        iter = iter + 1
    return Y_dash

def get_difference(X,Y):
    diff = 0
    for i in range(X.shape[1]):
        if X[0,i] != Y[0,i]:
            diff = diff + 1
    return diff

QBER = 0.1
n = 10
d_v = 4
d_c = 5
H, G = make_ldpc(n, d_v, d_c, systematic=True, sparse=True, seed=0)
X,Y = gen_sifted_keys(n,QBER)
#X = np.array([[1,1,0,1,0,0,1,0,0,1,0,0,0,0,0]])
#Y = np.array([[1,1,0,1,0,0,1,0,0,1,0,0,1,0,0]])
print(H)
C = get_CheckSum(H,X)
D = get_CheckSum(H,Y)
print(C)

max_iter = 10

diff = get_difference(X,Y)
print(X)
print("{0} Difference in bits: {1}".format(Y, diff))
print("###################################################")

Y_decode = decode_my(Y,H,C,D,max_iter,QBER,X)



#print("###################################################")
#Y_decode_alt = decode(H,Y,0.01)
#diff = get_difference(X,Y_decode_alt)
#print("{0} Difference in bits: {1}".format(Y_decode_alt, diff))



