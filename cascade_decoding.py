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



def binary_bit_flip(block, bit_positions, X, Y, Z, Eves_traces):
    ##split block in two
    new_blocksize = int(np.ceil(block.shape[0]/2))
    block_L = block[0:new_blocksize]
    bit_positions_L = bit_positions[0:new_blocksize]
    block_R = block[new_blocksize:]
    bit_positions_R = bit_positions[new_blocksize:]

    correct_L_parity, Z, Eves_traces = query_correct_parity(X,bit_positions_L, Eves_traces, Z)
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
            binary_bit_flip(block_L,bit_positions_L,X,Y,Z,Eves_traces)
    
    else:
        if new_blocksize == 1:
            print("---> Flip bit in right block\n")
            Y[int(bit_positions_R[0])] = (Y[int(bit_positions_R[0])] - 1) % 2##flip incorrect bit
            return Y
        else:
            print("---> binary algorithm on right sub block\n")
            binary_bit_flip(block_R,bit_positions_R,X,Y,Z,Eves_traces)
    return Y




def cascade_EC(X,Y,QBER, max_iter):

    Z = np.ones(X.shape)*-1 ###Eve's learned key. -1 is an unknown value
    Eves_traces = []

    N = X.shape[0]
    orig_perm = np.arange(0,N,1)
    iter = 0
    blocksize = 0
    while iter < max_iter:
        if iter == 0:
            perm = orig_perm
        else:
            #perm2 = np.asarray([9,2,0,7,15,11,14,10,1,12,4,6,17,13,18,19,5,8,16,3])
            #perm3 = np.asarray([11,16,18,1,19,10,4,17,0,7,6,3,15,8,12,9,5,13,2,14])
            perm = shuffle_key(perm) ##returns shuffled indexing of Y
            #if iter == 1:
             #   perm = perm2
            #if iter == 2:
              #  perm = perm3
        print("\n#### Iteration {0} ####\n{1} : Permutation of Bob's key".format(iter+1, perm))

        blocks, blocksize, bit_positions = split_blocks(N,iter,blocksize,QBER, perm)
        for i in range(len(blocks)):
            print("Block {0}: {1}  ---> Map to bits of Y in position: {2}".format(i, blocks[i], bit_positions[i]))

        current_parities = top_level_parity_check(blocks)
        correct_parities = []
        for bit_pos in bit_positions:
            ##we can pass Z and Eves_traces to learn key here as this attack is used on Alices hardware
            correct_parity, Z, Eves_traces = query_correct_parity(X,bit_pos,Eves_traces,Z)
            correct_parities.append(correct_parity)

        print("Current parities of Bob's blocks:   {0}".format(current_parities))
        print("Correct parities of Alice's blocks: {0}".format(correct_parities))

        for i in range(len(current_parities)):
            if current_parities[i] != correct_parities[i]:
                print("\n#### Preforming binary algorithm on block {0} #####".format(i))
                Y = binary_bit_flip(blocks[i], bit_positions[i], X, Y, Z, Eves_traces)
        
        print("{0} : Alice's original key".format(X))
        print("{0} : Bob's corrected key".format(Y))
        print("{0} : Eve's infered key". format(Z))
        iter = iter + 1

    return Y, Z, Eves_traces

def query_correct_parity(X,bit_positions,Eves_traces,Z):

    ##build bits to be checked into numpy array of correct form
    check_message = []
    for i in bit_positions:
        check_message.append(X[i])
    check_message = np.asarray(check_message)
  #  print(check_message)
    parity, hw = parity_check(check_message)

    print("\nparity check on bits at: {0}".format(bit_positions))
    print(check_message)
    print("   {0}".format(hw))
    print("\n")

    #################
    ## Side channel attack based on hw fits here

    ##let an entry in Eves traces correspond to both the hamming weights and corresponding bit positions and the corresponding parity bit
    if bit_positions.shape[0] > 1:
        trace_entry = [hw, bit_positions, parity] ##maybe let this run for ==1 condition too
        Eves_traces.append(trace_entry.copy()) 

    if bit_positions.shape[0] == 1:
        for trace_info in Eves_traces:
            Z = reveal_bits(Z,parity,bit_positions[0],trace_info)


    #################
    return parity, Z, Eves_traces

########################### Eves Tools ########################


def reveal_bits(Z, x_revealed, x_revealed_pos, trace_info):

    ## set x_e in Z
    Z[x_revealed_pos] = x_revealed

    trace = trace_info[0]
    bit_positions = trace_info[1]

    start = False ##variable to allow inference past error bit
    for i in range(bit_positions.shape[0]-1):
        ##start inference from revealed position
        if bit_positions[i] == x_revealed_pos or start == True:

            ##condition satisifed on first inference
          #  if start == False:
            if i == 0:
                if x_revealed == 0 and trace[i] == 0:
                    x_revealed_next = 0
                if x_revealed == 0 and trace[i] == 1:
                    x_revealed_next = 1
                if x_revealed == 1 and trace[i] == 0:
                    x_revealed_next = 1
                if x_revealed == 1 and trace[i] == 1:
                    x_revealed_next = 0

            ##condition satisfied on subsequent inferences
            else:
                if trace[i-1] == 0 and trace[i] == 0:
                    x_revealed_next = 0
                if trace[i-1] == 0 and trace[i] == 1:
                    x_revealed_next = 1
                if trace[i-1] == 1 and trace[i] == 0:
                    x_revealed_next = 1     
                if trace[i-1] == 1 and trace[i] == 1:
                    x_revealed_next = 0

  
            start = True
            Z[bit_positions[i+1]] = x_revealed_next

            



            ##set last parity bit in trace calculation to the output parity
            ##set the parity of the known bit and the next unknown bit to parity infered from the trace
         #   x_revealed = int(x_revealed) ^ int(x_r_p1_parity)
          #  Z[bit_positions[i+1]] = x_revealed
            """"
            if x_revealed == 0 and x_r_p1_parity == 0:
                Z[bit_positions[i+1]] = 0
            elif x_revealed == 0 and x_r_p1_parity == 1:
                Z[bit_positions[i+1]] = 1
            elif x_revealed == 1 and x_r_p1_parity == 0:
                Z[bit_positions[i+1]] = 1
            elif x_revealed == 1 and x_r_p1_parity == 1:
                Z[bit_positions[i+1]] = 0
            """
    return Z


def infer_remaining_bits(Z,Eves_traces,X):

    ##search for remaining unknown Z values
    for i in range(Z.shape[0]):
        if Z[i] == -1:
            for trace in Eves_traces:
                bit_positions = trace[1]
                if np.any(bit_positions == i) and bit_positions[0] != i:
                    parity_trace = trace[0]

                    j = np.where(bit_positions == i)[0][0]

                    ##for condition where only one parity trace bit is known
                    if len(parity_trace) == 1:
                        if j == 0:
                            x_known = Z[bit_positions[1]]
                        else:
                            x_known = Z[bit_positions[0]]

                        ##can only infer unkown bit if the other bit is known
                        if x_known != -1:
                            x_infer = int(x_known) ^ int(parity_trace[0])
                            Z[bit_positions[j]] = x_infer

                    ###must include condition that allows inference only if previous bit in parity check is known

                    if len(parity_trace) > 1 and Z[int(bit_positions[j-1])] != -1:                        
                      #  if parity_trace[j-2] == 0 and parity_trace[j-1] == 0:
                      #      infered_z = 0
                      #  if parity_trace[j-2] == 0 and parity_trace[j-1] == 1:
                      #      infered_z = 1
                      #  if parity_trace[j-2] == 1 and parity_trace[j-1] == 0:
                      #      infered_z = 1
                      #  if parity_trace[j-2] == 1 and parity_trace[j-1] == 1:
                      #      infered_z = 0
                        
                     #   Z[i] = infered_z
                        if j == 1:
                            Z_prev = Z[int(bit_positions[j-1])]
                            ##condition satisifed on first inference
                            
                            if Z_prev == 0 and parity_trace[j-1] == 0:
                                x_revealed_next = 0
                            if Z_prev == 0 and parity_trace[j-1] == 1:
                                x_revealed_next = 1
                            if Z_prev == 1 and parity_trace[j-1] == 0:
                                x_revealed_next = 1
                            if Z_prev == 1 and parity_trace[j-1] == 1:
                                x_revealed_next = 0


                       
                        
                        else:
                            if parity_trace[j-2] == 0 and parity_trace[j-1] == 0:
                                x_revealed_next = 0
                            if parity_trace[j-2] == 0 and parity_trace[j-1] == 1:
                                x_revealed_next = 1
                            if parity_trace[j-2] == 1 and parity_trace[j-1] == 0:
                                x_revealed_next = 1
                            if parity_trace[j-2] == 1 and parity_trace[j-1] == 1:
                                x_revealed_next = 0

                        Z[i] = x_revealed_next
                        if Z[i] != X[i]:
                            h = 9 ##debugg

                        break
    return Z

                    



################################################################



##### ToDo: Cascade effect #####

QBER = 0.1
N = 20
X,Y = gen_sifted_keys(N,QBER)
#X = np.asarray([0,0,0,1,1,1,1,0,0,0,1,0,0,1,1,1,1,0,0,0])
#Y = np.asarray([0,0,0,1,1,0,1,1,0,0,1,0,0,1,1,1,1,0,0,0])


#X is Alices key
#Y is Bobs key
print("{0} : Alice's original key".format(X))
print("{0} : Bob's original key".format(Y))
print("########### Beginnig Error Correction ###############\n\n")

Y_ec, Z, Eves_traces = cascade_EC(X,Y,QBER,3)

print("\n\n########### Finished Error Correction ###############\n")

print("{0} : Bobs final key".format(Y_ec))
if (X == Y_ec).all():
    print("Successfully Corrected Errors!")
else:
    print("Unsuccessful Error Correction")


print("\n\n########### Eves Infered Key ###############\n")
print("{0}".format(Z))

infered_total = 0
infered_correct = 0
for i in range(Z.shape[0]):
    if Z[i] != -1:
        infered_total = infered_total + 1
    if Z[i] == X[i]:
        infered_correct = infered_correct + 1

print("Total bits infered: {0}, {1}%% of which are correct".format(infered_total, (100*infered_correct/infered_total)))

Z = infer_remaining_bits(Z,Eves_traces,X)



infered_total = 0
infered_correct = 0
for i in range(Z.shape[0]):
    if Z[i] != -1:
        infered_total = infered_total + 1
    if Z[i] == X[i]:
        infered_correct = infered_correct + 1

print("\n{0}".format(Z))
print("Total bits infered after cleanup algorithm: {0}, {1}%% of which are correct".format(infered_total, (100*infered_correct/infered_total)))


Z = infer_remaining_bits(Z,Eves_traces,X)

infered_total = 0
infered_correct = 0
for i in range(Z.shape[0]):
    if Z[i] != -1:
        infered_total = infered_total + 1
    if Z[i] == X[i]:
        infered_correct = infered_correct + 1

print("\n{0}".format(Z))
print("Total bits infered after cleanup algorithm, 2nd iteration: {0}, {1}%% of which are correct".format(infered_total, (100*infered_correct/infered_total)))

