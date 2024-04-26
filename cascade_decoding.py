import numpy as np
import random
from colorama import Fore

def gen_sifted_keys(N,QBER):
    xA = np.zeros([N])
    xB = np.zeros([N])
    N_errs = int(np.floor(N*QBER))
    err_pos = random.sample(range(N),N_errs)
    for i in range(N):
        xA[i] = np.random.choice([0,1])
        xB[i] = xA[i]
        if i in err_pos:
            xB[i] = (xB[i] + 1) % 2
    return xA, xB

##calculates parity of block k
##outputs parity and also the estimates power trace of a register containing the current sum of the parity check
def parity_check(k):
    hw = []
    parity = int(k[0])
    for i in range(k.shape[0] - 1):
        parity = parity ^ int(k[i+1])
        if parity == 0:
            hw.append(0)
        else:
            hw.append(1)
    return parity, hw

def shuffle_key(perm):
  #  length = perm.shape[0]
    np.random.shuffle(perm)
  #  temp = np.zeros(key.shape) ##do not save key in shuffled form, save as a copy to retain original positions/
  #  for i in range(length):
  #      temp[i] = key[int(perm[i])]
  #  return temp, perm
    return perm

def split_blocks(N, iter , blocksize_prev, QBER, perm):
    
    ##set current blocksize (using floor func)
    blocksize = int(2*blocksize_prev)
    if iter == 0:
        blocksize = int(0.73/QBER)
    
    blocks = []
    bit_positions = []
    num_blocks = int(np.ceil(N/blocksize)) ##must round up to account for smaller end blocks
    for i in range(num_blocks):
        if i != num_blocks-1:
           # block_i = key[i*blocksize:(i+1)*blocksize]
            bit_positions_i = perm[i*blocksize:(i+1)*blocksize]
            block_i = np.zeros(bit_positions_i.shape)
            for j in range(bit_positions_i.shape[0]):
                block_i[j] = Y[int(bit_positions_i[j])]

        else:
            #block_i = key[i*blocksize:]
            bit_positions_i = perm[i*blocksize:]
            block_i = np.zeros(bit_positions_i.shape)
            for j in range(bit_positions_i.shape[0]):
                block_i[j] = Y[int(bit_positions_i[j])]

        blocks.append(block_i)
        bit_positions.append(bit_positions_i)
    return blocks, blocksize, bit_positions


##Bob only calls this function
def top_level_parity_check(blocks):
    parities = []
    for block in blocks:
        parity, hw = parity_check(block) #hw not needed yet. Usefull for attack on bobs side
        parities.append(parity)
    return parities

def query_correct_parity(X,bit_positions):

    ##build bits to be checked into numpy array of correct form
    check_message = []
    for i in bit_positions:
        check_message.append(X[i])
    check_message = np.asarray(check_message)
  #  print(check_message)
    parity, hw = parity_check(check_message)
    #################
    ## Side channel attack based on hw fits here





    #################
    return parity

def binary_bit_flip(block, bit_positions, X, Y):
    ##split block in two
    new_blocksize = int(np.ceil(block.shape[0]/2))
    block_L = block[0:new_blocksize]
    bit_positions_L = bit_positions[0:new_blocksize]
    block_R = block[new_blocksize:]
    bit_positions_R = bit_positions[new_blocksize:]

    correct_L_parity = query_correct_parity(X,bit_positions_L)
    current_L_parity, hw = parity_check(block_L)

    print("{0} : Left sub block.  ---> Map to bits of Y in position: {1}".format(block_L, bit_positions_L))
    print("{0} : Right sub block. ---> Map to bits of Y in position: {1}".format(block_R, bit_positions_R))
    print("Current left parity: {0}, Correct left parity {1}". format(current_L_parity, correct_L_parity))

    if correct_L_parity != current_L_parity:
        if new_blocksize == 1:
            print("---> Flip bit in left block\n")
            Y[int(bit_positions_L[0])] = (Y[int(bit_positions_L[0])] - 1) % 2 ##flip incorrect bit in Y_orig
            return Y
        else:
            print("---> binary algorithm on left sub block\n")
            binary_bit_flip(block_L,bit_positions_L,X,Y)
    
    else:
        if new_blocksize == 1:
            print("---> Flip bit in right block\n")
            Y[int(bit_positions_R[0])] = (Y[int(bit_positions_R[0])] - 1) % 2##flip incorrect bit
            return Y
        else:
            print("---> binary algorithm on right sub block\n")
            binary_bit_flip(block_R,bit_positions_R,X,Y)
    return Y




def cascade_EC(X,Y,QBER, max_iter):
    N = X.shape[0]
    orig_perm = np.arange(0,N,1)
    iter = 0
    blocksize = 0
    while iter < max_iter:
        if iter == 0:
            perm = orig_perm
        else:
            perm = shuffle_key(perm) ##returns shuffled indexing of Y

        print("\n#### Iteration {0} ####\n{1} : Permutation of Bob's key".format(iter+1, perm))

        blocks, blocksize, bit_positions = split_blocks(N,iter,blocksize,QBER, perm)
        for i in range(len(blocks)):
            print("Block {0}: {1}  ---> Map to bits of Y in position: {2}".format(i, blocks[i], bit_positions[i]))

        current_parities = top_level_parity_check(blocks)
        correct_parities = []
        for bit_pos in bit_positions:
            correct_parity = query_correct_parity(X,bit_pos)
            correct_parities.append(correct_parity)

        print("Current parities of Bob's blocks:   {0}".format(current_parities))
        print("Correct parities of Alice's blocks: {0}".format(correct_parities))

        for i in range(len(current_parities)):
            if current_parities[i] != correct_parities[i]:
                print("\n#### Preforming binary algorithm on block {0} #####".format(i))
                Y = binary_bit_flip(blocks[i], bit_positions[i], X, Y)
        
        print("{0} : Alice's original key".format(X))
        print("{0} : Bob's corrected key".format(Y))
        iter = iter + 1

    return Y



########################### Eves Tools ########################


def reveal_bits(Z, x_e, x_e_pos, bit_positions, trace, parity_bit):

    ## set x_e in Z
    Z[x_e_pos] = x_e

    for i in range(bit_positions):
        start = False ##variable to allow inference past error bit
        ##start inference from revealed position
        if bit_positions[i] == x_e_pos or start == True:

            if i == bit_positions - 1:
                x_e_p1 = parity_bit
            else:
                x_e_p1 = trace[i]
            start = True
            if Z[bit_positions[i]] == 0 and x_e_p1 == 0:
                Z[bit_positions[i]] = 0
            elif Z[bit_positions[i]] == 0 and x_e_p1 == 1:
                Z[bit_positions[i]] = 1
            elif Z[bit_positions[i]] == 1 and x_e_p1 == 0:
                Z[bit_positions[i]] = 1
            elif Z[bit_positions[i]] == 1 and x_e_p1 == 1:
                Z[bit_positions[i]] = 0
    return Z



################################################################







##### ToDo: Cascade effect #####

QBER = 0.3
N = 16
X,Y = gen_sifted_keys(N,QBER)
#X is Alices key
#Y is Bobs key
print("{0} : Alice's original key".format(X))
print("{0} : Bob's original key".format(Y))
print("########### Beginnig Error Correction ###############\n\n")

Y_ec = cascade_EC(X,Y,QBER,3)

print("\n\n########### Finished Error Correction ###############\n")

print("{0} : Bobs final key".format(Y_ec))
if (X == Y_ec).all():
    print("Successfully Corrected Errors!")
else:
    print("Unsuccessful Error Correction")

