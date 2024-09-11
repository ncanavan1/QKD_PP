import numpy as np
import random
import linalgtools
import noise_modelling

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
def parity_check(k,sigma):
    parity = 0
    pwr_trace = []
    pwr_thresh = -4.55 ##for stm32
    #pwr_thresh = 0.5
    for i in range(k.shape[0]):
        
        pwr = noise_modelling.simulate_xor_trace_stm32_kim(parity,int(k[i]),sigma)
        parity = parity ^ int(k[i])
        
        if pwr <= pwr_thresh + sigma*0.8:
            pwr_trace.append(0)
        else:
            pwr_trace.append(1)
    return parity, pwr_trace

def shuffle_key(perm):
  #  length = perm.shape[0]
    np.random.shuffle(perm)
  #  temp = np.zeros(key.shape) ##do not save key in shuffled form, save as a copy to retain original positions/
  #  for i in range(length):
  #      temp[i] = key[int(perm[i])]
  #  return temp, perm
    return perm

def split_blocks(N, iter , blocksize_prev, QBER, perm, Y):
    
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
def top_level_parity_check(blocks, sigma):
    parities = []
    for block in blocks:
        parity, pwr_trace = parity_check(block, sigma) #hw not needed yet. Usefull for attack on bobs side
        parities.append(parity)
    return parities



def binary_bit_flip(block, bit_positions, X, Y, Eves_traces,sigma,record):
    ##split block in two
    new_blocksize = int(np.ceil(block.shape[0]/2))
    block_L = block[0:new_blocksize]
    bit_positions_L = bit_positions[0:new_blocksize]
    block_R = block[new_blocksize:]
    bit_positions_R = bit_positions[new_blocksize:]

    correct_L_parity, Eves_traces = query_correct_parity(X,bit_positions_L, Eves_traces,sigma,record)
    current_L_parity, pwr_trace = parity_check(block_L, sigma)

    #print("{0} : Left sub block.  ---> Map to bits of Y in position: {1}".format(block_L, bit_positions_L))
    #print("{0} : Right sub block. ---> Map to bits of Y in position: {1}".format(block_R, bit_positions_R))
    #print("Current left parity: {0}, Correct left parity {1}". format(current_L_parity, correct_L_parity))

    if correct_L_parity != current_L_parity:
        if new_blocksize == 1:
            #print("---> Flip bit in left block\n")
            Y[int(bit_positions_L[0])] = (Y[int(bit_positions_L[0])] - 1) % 2 ##flip incorrect bit in Y_orig
            return Y

        else:
           # print("---> binary algorithm on left sub block\n")
            binary_bit_flip(block_L,bit_positions_L,X,Y,Eves_traces,sigma,record)
    
    else:
        if new_blocksize == 1:
          #  print("---> Flip bit in right block\n")
            Y[int(bit_positions_R[0])] = (Y[int(bit_positions_R[0])] - 1) % 2##flip incorrect bit
            return Y
        else:
         #   print("---> binary algorithm on right sub block\n")
            binary_bit_flip(block_R,bit_positions_R,X,Y,Eves_traces,sigma,record)
    return Y


def cascade_effect(X, key_index_list, key_value_list, corrected_pos, block_list, bit_positions, Eves_traces, Eve_mode, sigma, record):
    iterations = len(key_index_list)
    key_curr = key_index_list[-1]
    #key_prev = shuffled_key_list[-2]

    for key in range(len(key_index_list)-1):
            ##reconstruct blocks
            for i in range(len(bit_positions[key])):
                for j in range(len(bit_positions[key][i])):
                    block_list[key][i][j] = key_value_list[key][bit_positions[key][i][j]]

    Y_old = key_value_list[-2].copy()
    key_value_list[-2], Eves_traces = error_detect_correct_blocks(X,key_value_list[-2],block_list[-2],bit_positions[-2],Eves_traces,Eve_mode,sigma,record)

    ##record what bits have been flipped
    corrected_pos_new = []
    for i in range(len(Y_old)):
        if key_value_list[-2][i] != Y_old[i]:
            corrected_pos_new.append(i)

    if len(key_index_list) == 2:
        return key_value_list[-2], Eves_traces
    
    else:
        cascade_effect(X, key_index_list[:-1], key_value_list[:-1], corrected_pos[:-1], block_list[:-1], bit_positions[:-1], Eves_traces, Eve_mode, sigma, record)
        return key_value_list[-2], Eves_traces

def error_detect_correct_blocks(X,Y,blocks,bit_positions,Eves_traces,Eve_mode,sigma,record):
        current_parities = top_level_parity_check(blocks, sigma)
        correct_parities = []
        for block_pos in bit_positions:

            if Eve_mode == 0:
                if iter == 0:
                    record = True
                else:
                    record = False
            
            correct_parity, Eves_traces = query_correct_parity(X,block_pos,Eves_traces,sigma,record)
            correct_parities.append(correct_parity)

        for i in range(len(current_parities)):
            if current_parities[i] != correct_parities[i]:
        
                if Eve_mode == 0: ##top block trace recording only
                    record = False
                Y = binary_bit_flip(blocks[i], bit_positions[i], X, Y, Eves_traces,sigma,record)
        return Y, Eves_traces

def cascade_EC(X, Y, QBER, max_iter, sigma, Eve_mode):

    Eves_traces = []
    N = X.shape[0]
    orig_perm = np.arange(0,N,1)
    iter = 0
    blocksize = 0
    permutations_list = []
    iteration_blocks_pos = []
    iteration_correction = []
    key_value_list = []
    block_list = []


    record = True
    if Eve_mode == 0:
        record = False
    elif Eve_mode == 1:
        record == True

    while iter < max_iter:
        if iter == 0:
            perm = orig_perm
            permutations_list.append(perm)

        else:
            perm = shuffle_key(perm.copy()) ##returns shuffled indexing of Y
            permutations_list.append(perm)

        blocks, blocksize, bit_positions = split_blocks(N,iter,blocksize,QBER, perm,Y)
        block_list.append(blocks)
        iteration_blocks_pos.append(bit_positions)

        Y_old = Y.copy()

     
        Y, Eves_traces = error_detect_correct_blocks(X,Y,blocks,bit_positions,Eves_traces,Eve_mode,sigma,record)

        key_value_list.append(Y)
        ##record what bits have been flipped
        corrected = []
        for i in range(len(Y)):
            if Y[i] != Y_old[i]:
                corrected.append(i)
        iteration_correction.append(corrected)

        if iter > 0:
            Y, Eves_traces = cascade_effect(X,permutations_list,key_value_list,iteration_correction, block_list ,iteration_blocks_pos, Eves_traces,Eve_mode, sigma, record)

        iter = iter + 1

    return Y, Eves_traces

def query_correct_parity(X,bit_positions,Eves_traces,sigma,record):

    ##build bits to be checked into numpy array of correct form
    check_message = []
    for i in bit_positions:
        check_message.append(X[i])
    check_message = np.asarray(check_message)
    parity, pwr_trace = parity_check(check_message, sigma)


    trace_entry = [pwr_trace, bit_positions, parity]
    if record == True:
        Eves_traces.append(trace_entry.copy()) 

    #################
    return parity, Eves_traces



X = np.asarray([0,1,1,0,0,0,1,1,0,1,1])
Y = np.asarray([0,1,1,1,0,0,1,0,0,1,1])
Y_Ec, Eves_Traces = cascade_EC(X,Y,0.1,2,0.1,1)

EveKey = np.asarray([0,-1,1,-1,0,-1,1,-1,0,-1,1])