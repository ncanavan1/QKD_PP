import numpy as np
from pyldpc import utils

np.random.seed(0)


def gen_ldpc_matrix(m,n,d):
    h = np.zeros([m,n])
    for i in range(0,m):
        used_index = []
        while len(used_index) < d:
            next_set = np.random.randint(0,n)
            if used_index.__contains__(next_set) == False:
                h[i,next_set] = 1
                used_index.append(next_set)
    return h

def parity_check_matrix(n_code, d_v, d_c, seed=None):
    """
    Build a regular Parity-Check Matrix H following Callager's algorithm.

    Parameters
    ----------
    n_code: int, Length of the codewords.
    d_v: int, Number of parity-check equations including a certain bit.
        Must be greater or equal to 2.
    d_c: int, Number of bits in the same parity-check equation. d_c Must be
        greater or equal to d_v and must divide n.
    seed: int, seed of the random generator.

    Returns
    -------
    H: array (n_equations, n_code). LDPC regular matrix H.
        Where n_equations = d_v * n / d_c, the total number of parity-check
        equations.

    """
    rng = utils.check_random_state(seed)

    if d_v <= 1:
        raise ValueError("""d_v must be at least 2.""")

    if d_c <= d_v:
        raise ValueError("""d_c must be greater than d_v.""")

 #   if n_code % d_c:
 #       raise ValueError("""d_c must divide n for a regular LDPC matrix H.""")

    n_equations = (n_code * d_v) // d_c

    block = np.zeros((n_equations // d_v, n_code), dtype=int)
    H = np.empty((n_equations, n_code))
    block_size = n_equations // d_v

    # Filling the first block with consecutive ones in each row of the block

    for i in range(block_size):
        for j in range(i * d_c, (i+1) * d_c):
            block[i, j] = 1
    H[:block_size] = block

    # reate remaining blocks by permutations of the first block's columns:
    for i in range(1, d_v):
        H[i * block_size: (i + 1) * block_size] = rng.permutation(block.T).T
    H = H.astype(int)
    return H

def gen_x(n):
    x = np.zeros([1,n])
    for i in range(0,n):
        x[0,i] = np.random.randint(0,2)
    return x

def reduce_h(h,d):
    h_dash = np.zeros([h.shape[0],d])
    for i in range(0,h.shape[0 ]):
        count = 0
        for j in range(0,h.shape[1]):
            if h[i,j] == 1:
                h_dash[i,count] = j
                count = count + 1
    return h_dash


############Hamming weight model goes here
def calc_syndrome_case1(h_dash,m,d,x):

    xor_recordings = np.zeros([m,d]) ##stores 0 if jth xor on ith value of s is 0, 1 otherwise

    s = np.zeros([m,1])
    for i in range(0,m):
        for j in range(0,d):
            x_index = int(h_dash[i,j])
            s[i,0] = int(s[i,0]) ^ int(x[0,x_index])
            if s[i,0] == 1:
                xor_recordings[i,j] = 1
    return s, xor_recordings




n = 16
m = 8
d = 6
#h, g = code.make_ldpc(n,d,m)
#h = parity_check_matrix(16,d,8)
h_hard = [[0,0,1,0,0,1,1,1,0,0,0,0],
          [1,1,0,0,1,0,0,0,0,0,0,1],
          [0,0,0,1,0,0,0,0,1,1,1,0],
          [0,1,0,0,0,1,1,0,0,1,0,0],
          [1,0,1,0,0,0,0,1,0,0,1,0],
          [0,0,0,1,1,0,0,0,1,0,0,1],
          [1,0,0,1,1,0,1,0,0,0,0,0],
          [0,0,0,0,0,1,0,1,0,0,1,1],
          [0,1,1,0,0,0,0,0,1,1,0,0]]

h_hard = [[0,1,0,0,0,1,1,0,0,0,0,0,0,1,1,1],
           [0,1,1,1,0,0,0,0,0,1,0,0,0,1,1,0],
           [0,1,0,0,0,1,1,0,0,1,0,1,0,0,0,1],
           [0,0,1,0,1,0,1,1,1,0,0,1,0,0,0,0],
           [0,0,1,1,1,0,0,0,0,0,1,0,1,0,0,0],
           [1,0,1,1,1,0,0,0,0,0,1,0,1,0,0,0],
           [1,0,0,0,0,1,0,0,1,0,1,1,1,0,0,0],
           [1,0,0,0,0,0,0,1,0,1,1,0,1,1,0,0],
           [0,0,0,1,1,0,0,1,1,0,0,0,0,0,1,1]]
h = np.asarray(h_hard)






print(h.transpose())
#h = gen_ldpc_matrix(m,n,d)
#print(h)
x = gen_x(n)
print("Code word projection: (0 vector for valid code words w.r.t H)\n", np.matmul(h,x.transpose())%2)

h_dash = reduce_h(h,d)

s1, xor_record = calc_syndrome_case1(h_dash,m,d,x)
#print(s1)



###Checks s1 calculation against full matrix multiplication
s_real = np.matmul(h,x.transpose())%2
if(s_real == s1).all():
    print("Valid syndrome calculation")
else:
    print("Invalid syndrome calculation")




###Assume now that adversary has access to xor_record matrix (equivalent to sub trace clustering values) and h_bar.
###They wish to recover s and wish to recover x


def recover_x_and_s_s1(xor_record, h_dash, m, n, d):
    x_recovered = np.zeros([n,1])
    s_recovered = np.zeros([m,1])
    for i in range(0,m):
        prev_s = 0
        prev_x = 0
        for j in range(0,d):
            if j == 0:
                prev_x = xor_record[i,j]
                prev_s = prev_x
                x_index = int(h_dash[i,j])
                x_recovered[x_index] = prev_x

            if j > 0 and j < d-1:
                prev_x = int(prev_s) ^ int(xor_record[i,j])
                prev_s = prev_x
                x_index = int(h_dash[i,j])
                x_recovered[x_index] = prev_x

            if j == d-1:
                prev_x = int(prev_s) ^ int(xor_record[i,j])
                s_recovered[i,0] = xor_record[i,j]
                x_index = int(h_dash[i,j])
                x_recovered[x_index] = int(x_recovered[x_index]) ^ int(prev_x)
    return x_recovered, s_recovered

x_recovered, s_recovered = recover_x_and_s_s1(xor_record, h_dash, m, n,d)

if(x == x_recovered).all():
    print("Successful x recovery")
else:
    print("Failed x recovery")
    print(x.tolist())
    print(x_recovered.tolist())


if(s1 == s_recovered).all():
    print("Successful s recovery")
else:
    print("Failed s recovery")

print(h_dash)
print(h)

