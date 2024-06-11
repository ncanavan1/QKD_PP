import numpy as np
import matplotlib.pyplot as plt
import random
from collections import Counter

#########################################################################
#########################################################################
#########################################################################

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
    parity = 0
    for i in range(k.shape[0]):
        parity = parity ^ int(k[i])
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
def top_level_parity_check(blocks):
    parities = []
    for block in blocks:
        parity, hw = parity_check(block) #hw not needed yet. Usefull for attack on bobs side
        parities.append(parity)
    return parities



def binary_bit_flip(block, bit_positions, X, Y, Eves_traces,XOR_noise,record):
    ##split block in two
    new_blocksize = int(np.ceil(block.shape[0]/2))
    block_L = block[0:new_blocksize]
    bit_positions_L = bit_positions[0:new_blocksize]
    block_R = block[new_blocksize:]
    bit_positions_R = bit_positions[new_blocksize:]

    correct_L_parity, Eves_traces = query_correct_parity(X,bit_positions_L, Eves_traces,XOR_noise,record)
    current_L_parity, hw = parity_check(block_L)

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
            binary_bit_flip(block_L,bit_positions_L,X,Y,Eves_traces,XOR_noise,record)
    
    else:
        if new_blocksize == 1:
          #  print("---> Flip bit in right block\n")
            Y[int(bit_positions_R[0])] = (Y[int(bit_positions_R[0])] - 1) % 2##flip incorrect bit
            return Y
        else:
         #   print("---> binary algorithm on right sub block\n")
            binary_bit_flip(block_R,bit_positions_R,X,Y,Eves_traces,XOR_noise,record)
    return Y




def cascade_EC(X,Y,QBER, max_iter,XOR_noise,Eve_mode):

    Eves_traces = []

    record = True
    if Eve_mode == 0:
        record = False
    elif Eve_mode == 1:
        record == True

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
            perm = shuffle_key(perm.copy()) ##returns shuffled indexing of Y
            #if iter == 1:
             #   perm = perm2
            #if iter == 2:
              #  perm = perm3
        #print("\n#### Iteration {0} ####\n{1} : Permutation of Bob's key".format(iter+1, perm))

        blocks, blocksize, bit_positions = split_blocks(N,iter,blocksize,QBER, perm,Y)
        #for i in range(len(blocks)):
         #   print("Block {0}: {1}  ---> Map to bits of Y in position: {2}".format(i, blocks[i], bit_positions[i]))

        current_parities = top_level_parity_check(blocks)
        correct_parities = []
        for bit_pos in bit_positions:
            ##we can pass Z and Eves_traces to learn key here as this attack is used on Alices hardware
            
            if Eve_mode == 0:
                if iter == 0:
                    record = True
                else:
                    record = False
            
            correct_parity, Eves_traces = query_correct_parity(X,bit_pos,Eves_traces,XOR_noise,record)
            correct_parities.append(correct_parity)

        #print("Current parities of Bob's blocks:   {0}".format(current_parities))
        #print("Correct parities of Alice's blocks: {0}".format(correct_parities))

        for i in range(len(current_parities)):
            if current_parities[i] != correct_parities[i]:
        
        ##      print("\n#### Preforming binary algorithm on block {0} #####".format(i))
                if Eve_mode == 0: ##top block trace recording only
                    record = False
                Y = binary_bit_flip(blocks[i], bit_positions[i], X, Y, Eves_traces,XOR_noise,record)
        
        #print("{0} : Alice's original key".format(X))
        #print("{0} : Bob's corrected key".format(Y))
        iter = iter + 1

    return Y, Eves_traces

def query_correct_parity(X,bit_positions,Eves_traces,XOR_noise,record):

    ##build bits to be checked into numpy array of correct form
    check_message = []
    for i in bit_positions:
        check_message.append(X[i])
    check_message = np.asarray(check_message)
    parity, hw = parity_check(check_message)
    for i in range(len(hw)):
        ###add noise effect
        noise_p = np.random.uniform(0,1)
        if noise_p < XOR_noise:
            hw[i] = int(hw[i]) ^ 1

    trace_entry = [hw, bit_positions, parity]
    if record == True:
        Eves_traces.append(trace_entry.copy()) 

    #################
    return parity, Eves_traces

########################### Eves Tools ########################

               



################################################################




#########################################################################
#########################################################################
#########################################################################


def gen_secret(N):
    xA = np.zeros([N])
    for i in range(N):
        xA[i] = np.random.choice([0,1])
    return xA
 
def gen_Trace(secret,N,err):
    trace = np.zeros([N])
    trace[0] = int(secret[0]) ^ 0
    for i in range(1,N):
        trace[i] = int(trace[i-1]) ^ int(secret[i])
 
    N_errs = int(np.floor(N*err))
    err_pos = random.sample(range(N),N_errs)
    for i in range(N):
        if i in err_pos:
            trace[i] = (trace[i] + 1) % 2
    return trace
 
def recover_secret_majoirty(traces,N,sigma):
    eve_key_all = np.empty(N,dtype=object)
    for i in range(0,N):
        eve_key_all[i] = []

    for trace in traces:
        parity = trace[0]
        pos = trace[1]
        eve_key_all[pos[0]].append(int(parity[0]) ^ 0)
        for i in range(1,pos.shape[0]):
            eve_key_all[pos[i]].append(int(parity[i]) ^ int(parity[i-1]))
 
    final_key = np.zeros(N)
    for i in range(N):
        ones = 0
        zeors = 0
        for vote in eve_key_all[i]:
            if vote == 0:
                zeors = zeors + 1
            else:
                ones = ones + 1
        if ones > zeors:
            final_key[i] = 1
        else:
            final_key[i] = 0
 
    return final_key, eve_key_all


def Eve_parity_checks(traces,eve_key):
    ##loop through all traces and compare Eves parity with Alices parity
    failed = []
    failed_pos = []
    for trace in traces:
        pos = trace[1]
        parityAlice = trace[2]
        block = []
        for p in pos:
            block.append(eve_key[p])
        parityEve, hw = parity_check(np.asarray(block))
        if parityAlice != parityEve:
            failed.append(trace)
            failed_pos.extend(pos)


    return failed, failed_pos, len(failed)

        


def confidence_flip(eve_key_all, eve_key_maj, traces, N, thresh):
    conf_rating = np.ones(N) * 0.5
    err_c = 0
    for trace in traces:
        if len(trace[0]) == 1:
            
            if eve_key_maj[int(trace[1][0])] != trace[0][0]:
                err_c = err_c + 1

            eve_key_maj[int(trace[1][0])] = trace[0][0]

            if trace[0][0] == 0:
                conf_rating[int(trace[1][0])] = 0
            else:
                conf_rating[int(trace[1][0])] = 1

    
    for i in range(N):
        bit_avg = np.average(np.asarray(eve_key_all[i]))
        conf_rating[i] = bit_avg

    ###Insert parity sum error checks
    ###Set threshold error to iterate through
    failed_parity_traces, failed_pos, N_Failed_min = Eve_parity_checks(traces,eve_key_maj)
    sorted_failed = [item for items, c in Counter(failed_pos).most_common() for item in [items] * c]
    failed_ranked = []
    for pos in sorted_failed:
        if not failed_ranked.__contains__(pos):
            failed_ranked.append(pos)


    ###for the top most ranked bits in failure checks, generate all possible variations

    eve_key_new = eve_key_maj.copy()
    for p in failed_ranked[0:5]:
        eve_key_new[p] = int(eve_key_new[p]) ^ 1

    failed_parity_traces, failed_pos, N_Failed_new = Eve_parity_checks(traces,eve_key_new)

    
    return eve_key_new
 
def recover_secret_basic(traces,N,sigma):
    eve_key = np.zeros(N)
    
    for trace in traces:
        parity = trace[0]
        pos = trace[1]
        eve_key[pos[0]] = int(parity[0]) ^ 0
        for i in range(1,pos.shape[0]):
            eve_key[pos[i]] = int(parity[i]) ^ int(parity[i-1])
    return eve_key
 

"""
N = 10000
secret = gen_secret(N)
err_range = np.arange(0,0.2001,0.01)
sr_maj = []
sr_niev = []
for err in err_range:
    traces = []
    for i in range(5):
        traces.append(gen_Trace(secret,N,err))
    eve_key_maj = recover_secret_majoirty(traces,N,0)
    eve_key_niev = recover_secret_basic(traces[0],N,0)
 
    success_count_maj = 0
    success_count_niev = 0
 
    for i in range(N):
        if secret[i] == eve_key_maj[i]:
            success_count_maj = success_count_maj + 1
    success_rate_maj = success_count_maj/N
    sr_maj.append(success_rate_maj)
 
    for i in range(N):
        if secret[i] == eve_key_niev[i]:
            success_count_niev = success_count_niev + 1
    success_rate_niev = success_count_niev/N
    sr_niev.append(success_rate_niev)
 
 
plt.plot(err_range,sr_maj)
plt.plot(err_range,sr_niev)
 
 
plt.show()
#print("Success rate: {0}".format(success_rate))
"""

def run_EC_majority(N,QBER,XOR_noise,max_iter):
    X,Y = gen_sifted_keys(N,QBER)
    Y_ec, Eves_traces = cascade_EC(X,Y,QBER,max_iter,XOR_noise,1)
    if (X == Y_ec).all():
        print("Successfully Corrected Errors!")
        ###Run Attack Here

        final_key, eve_key_all = recover_secret_majoirty(Eves_traces,N,0)
        success_count = 0
        for i in range(N):
            if X[i] == final_key[i]:
                success_count = success_count + 1
        success_rate = success_count/N
        return success_rate
       # print("Eve Success rate: {0}".format(success_rate))
       # print(final_key)
        
    else:
        #print("Unsuccessful Error Correction")
        return -1

def run_EC_topblock_only(N,QBER,XOR_noise,max_iter):
    X,Y = gen_sifted_keys(N,QBER)
    Y_ec, Eves_traces = cascade_EC(X,Y,QBER,max_iter,XOR_noise,0)
    if (X == Y_ec).all():
        #print("Successfully Corrected Errors!")
        ###Run Attack Here

        final_key = recover_secret_basic(Eves_traces,N,0)
        success_count = 0
        for i in range(N):
            if X[i] == final_key[i]:
                success_count = success_count + 1
        success_rate = success_count/N
        #print("Eve Success rate: {0}".format(success_rate))
       # print(final_key)
        return success_rate
    else:
#        print("Unsuccessful Error Correction")
        return -1
    
def run_EC_confidence_flip(N,QBER,XOR_noise,max_iter):
    X,Y = gen_sifted_keys(N,QBER)
    Y_ec, Eves_traces = cascade_EC(X,Y,QBER,max_iter,XOR_noise,1)
    if (X == Y_ec).all():
        #print("Successfully Corrected Errors!")
        ###Run Attack Here

        final_key,eve_key_all = recover_secret_majoirty(Eves_traces,N,0)
        conf_key = confidence_flip(eve_key_all,final_key,Eves_traces,N,0.1)
        success_count = 0
        for i in range(N):
            if X[i] == conf_key[i]:
                success_count = success_count + 1
        success_rate = success_count/N
        #print("Eve Success rate: {0}".format(success_rate))
       # print(final_key)
        return success_rate
    else:
#        print("Unsuccessful Error Correction")
        return -1
    

def plot_params(N,QBER,max_iter, err_range):
    tb_only = []
    majority = []
    conf = []

    for err in err_range:
        sr_tb = -1
        sr_maj = -1
        sr_conf = -1
        while sr_tb == -1:
            sr_tb = run_EC_topblock_only(N,QBER,err,max_iter)
        tb_only.append(sr_tb)

        while sr_maj == -1:
            sr_maj = run_EC_majority(N,QBER,err,max_iter)
        majority.append(sr_maj)

        while sr_conf == -1:
            sr_conf = run_EC_confidence_flip(N,QBER,err,max_iter)
        conf.append(sr_conf)

    return tb_only, majority, conf




err_range = np.arange(0,0.2001,0.01)
N = 512
max_iter = 5
QBER = 0.01
avg_iter = 10
tb_avg = np.zeros(shape=(21,avg_iter))
maj_avg = np.zeros(shape=(21,avg_iter))
conf_avg = np.zeros(shape=(21,avg_iter))

for i in range(avg_iter):
    tb_only, majority, conf = plot_params(N,QBER,max_iter,err_range)
    tb_avg[:,i] = np.asarray(tb_only)
    maj_avg[:,i] = np.asarray(majority)
    conf_avg[:,i] = np.asarray(conf)
majority = np.average(maj_avg,axis=1)
tb_only = np.average(tb_avg,axis=1)
conf = np.average(conf_avg,axis=1)

plt.plot(err_range,majority,label="Majority Rule, QBER={0}".format(QBER))
plt.plot(err_range,tb_only, label="Top Block Only, QBER={0}".format(QBER))
plt.plot(err_range,conf, label="Confidence, QBER={0}".format(QBER))

QBER = 0.1
tb_avg = np.zeros(shape=(21,avg_iter))
maj_avg = np.zeros(shape=(21,avg_iter))
conf_avg = np.zeros(shape=(21,avg_iter))

for i in range(avg_iter):
    tb_only, majority,conf = plot_params(N,QBER,max_iter,err_range)
    tb_avg[:,i] = np.asarray(tb_only)
    maj_avg[:,i] = np.asarray(majority)
    conf_avg[:,i] = np.asarray(conf)
majority = np.average(maj_avg,axis=1)
tb_only = np.average(tb_avg,axis=1)
conf = np.average(conf_avg,axis=1)

plt.plot(err_range,majority,label="Majority Rule, QBER={0}".format(QBER))
plt.plot(err_range,tb_only, label="Top Block Only, QBER={0}".format(QBER))
plt.plot(err_range,conf, label="Confidence, QBER={0}".format(QBER))

plt.xlabel("Trace Error Rate")
plt.ylabel("Key Recovery Success Rate")
plt.title("N={0}, Max Iterations={1}".format(N,max_iter))
plt.legend()
plt.savefig("results/Eve_strat_N_{0}_maxiter_{1}_conf.png".format(N,max_iter))
