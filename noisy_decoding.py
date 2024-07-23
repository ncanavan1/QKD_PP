import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpathces
import random
from decimal import Decimal
import math
from collections import Counter
import linalgtools as my_solver
from gf2 import solve_gf2

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


def cascade_effect(X, key_index_list, key_value_list, corrected_pos, block_list, bit_positions, Eves_traces, Eve_mode, XOR_noise, record):
    iterations = len(key_index_list)
    key_curr = key_index_list[-1]
    #key_prev = shuffled_key_list[-2]

    for key in range(len(key_index_list)-1):
            ##reconstruct blocks
            for i in range(len(bit_positions[key])):
                for j in range(len(bit_positions[key][i])):
                    block_list[key][i][j] = key_value_list[key][bit_positions[key][i][j]]

    Y_old = key_value_list[-2].copy()
    key_value_list[-2], Eves_traces = error_detect_correct_blocks(X,key_value_list[-2],block_list[-2],bit_positions[-2],Eves_traces,Eve_mode,XOR_noise,record)

    ##record what bits have been flipped
    corrected_pos_new = []
    for i in range(len(Y_old)):
        if key_value_list[-2][i] != Y_old[i]:
            corrected_pos_new.append(i)

    if len(key_index_list) == 2:
        return key_value_list[-2], Eves_traces
    
    else:
        cascade_effect(X, key_index_list[:-1], key_value_list[:-1], corrected_pos[:-1], block_list[:-1], bit_positions[:-1], Eves_traces, Eve_mode, XOR_noise, record)
        return key_value_list[-2], Eves_traces

def error_detect_correct_blocks(X,Y,blocks,bit_positions,Eves_traces,Eve_mode,XOR_noise,record):
        current_parities = top_level_parity_check(blocks)
        correct_parities = []
        for block_pos in bit_positions:

            if Eve_mode == 0:
                if iter == 0:
                    record = True
                else:
                    record = False
            
            correct_parity, Eves_traces = query_correct_parity(X,block_pos,Eves_traces,XOR_noise,record)
            correct_parities.append(correct_parity)

        for i in range(len(current_parities)):
            if current_parities[i] != correct_parities[i]:
        
                if Eve_mode == 0: ##top block trace recording only
                    record = False
                Y = binary_bit_flip(blocks[i], bit_positions[i], X, Y, Eves_traces,XOR_noise,record)
        return Y, Eves_traces

def cascade_EC(X, Y, QBER, max_iter, XOR_noise, Eve_mode):

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

     
        Y, Eves_traces = error_detect_correct_blocks(X,Y,blocks,bit_positions,Eves_traces,Eve_mode,XOR_noise,record)

        key_value_list.append(Y)
        ##record what bits have been flipped
        corrected = []
        for i in range(len(Y)):
            if Y[i] != Y_old[i]:
                corrected.append(i)
        iteration_correction.append(corrected)

        if iter > 0:
            Y, Eves_traces = cascade_effect(X,permutations_list,key_value_list,iteration_correction, block_list ,iteration_blocks_pos, Eves_traces,Eve_mode, XOR_noise, record)

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

def top_offenders(N,traces,P_A,P_E):
    num_checks = len(traces)
    offence_count = np.zeros(N)
    pos_list = np.arange(0,N,1)
    for i in range(num_checks):
        if P_A[i] != P_E[i]:
            offence_trace = traces[i][1]
            for pos in offence_trace:
                offence_count[pos] = offence_count[pos] + 1
    combined =  np.vstack((pos_list,offence_count)).T
    combined2 = np.argsort(-combined[:,1])
    combined = combined[combined2]
    return combined
    

def recover_secret_majoirty(traces,N,sigma):
    M = len(traces)
    H_E = np.zeros([M,N])
    P_A = np.zeros(M)
    Included_matrix = np.zeros([M,N])
    Conf_vector = np.zeros(N)
    final_key = np.zeros(N)
    alpha = 0.4 ##emperic
    beta = 2.5


    check = 0
    for trace in traces:
        inter_parity = trace[0]
        pos = trace[1]
        out_parity = trace[2]
        P_A[check] = out_parity
        if pos.size == 1:
            Conf_vector[pos[0]] = 1
            final_key[pos[0]] = out_parity
        for i in reversed(range(0,pos.size)): ##this was 1,pos.size WHY??
            H_E[check,pos[i]] = inter_parity[i] ^ inter_parity[i-1]
            Included_matrix[check,pos[i]] = 1
        H_E[check,pos[0]] = inter_parity[0]
        check = check + 1

   # test = Included_matrix[:,0:check]
    #incl_inv = np.linalg.inv(test)
    
    random_choice_pos = []

    for i in range(N):
        if Conf_vector[i] != 1:
            ones = 0
            zeros = 0
            for check in range(M):
                if Included_matrix[check,i] == 1:
                    if H_E[check,i] == 0:
                        zeros = zeros + 1
                    else:
                        ones = ones + 1
            tot = ones+zeros
            if ones > zeros:
                final_key[i] = 1
              #  Conf_vector[i] = ((ones/tot)**(1-sigma) * (zeros/tot)**(sigma)) - ((zeros/tot)**(1-sigma) * (ones/tot)**(sigma))
               # Conf_vector[i] = ((1-sigma)**(ones) * sigma**(zeros)) / ((1-sigma)**(zeros) * sigma**(ones))
               # Conf_vector[i] = math.log10(Conf_vector[i])
                #Conf_vector[i] = np.tanh(alpha*(ones-zeros))**(1/(1-sigma)**2)
                Conf_vector[i] = (1-sigma)**(ones/(ones+zeros)) * (sigma)**(zeros/(ones+zeros))

            elif zeros > ones:
                final_key[i] = 0
               # Conf_vector[i] = ((zeros/tot)**(1-sigma) * (ones/tot)**(sigma)) - ((ones/tot)**(1-sigma) * (zeros/tot)**(sigma))
             #   Conf_vector[i] = ((1-sigma)**(zeros) * sigma**(ones))/((1-sigma)**(ones) * sigma**(zeros))
               # Conf_vector[i] = math.log10(Conf_vector[i])
                #Conf_vector[i] = np.tanh(alpha*(zeros-ones))**(1/(1-sigma)**2)
                Conf_vector[i] = (1-sigma)**(zeros/(ones+zeros)) * (sigma)**(ones/(ones+zeros))



            else: #ones==zeros
                #final_key[i] = random.choice([0,1])
                #Conf_vector[i] = 0
                final_key[i] = -1
                random_choice_pos.append(i)

                if final_key[i] == 0:
                    Conf_vector[i] = (1-sigma)**(zeros/(ones+zeros)) * (sigma)**(ones/(ones+zeros))
                else:
                    Conf_vector[i] = (1-sigma)**(ones/(ones+zeros)) * (sigma)**(zeros/(ones+zeros))
                Conf_vector[i] = 0.1

    if len(random_choice_pos) > 0:
        print("Key contains {0} random choices".format(len(random_choice_pos)))
        
        ####ensure that we don't try too large spaces
        if len(random_choice_pos) <= 12:
            diff, min_list = test_random_choices(final_key,random_choice_pos,traces)
        else:
            for pos in random_choice_pos:
                final_key[pos] = random.choice([0,1])


    return final_key, H_E, Included_matrix, Conf_vector, P_A

def test_random_choices(key,random_pos,traces):
    rand_len = len(random_pos)
    rand_combs = gen_binary_comb(rand_len)
    min_list = []

    for pos in random_pos:
        key[pos] = 0
    failed, failed_pos, diff = Eve_parity_checks(traces,key)
    min_list.append(rand_combs[0])

    for comb in rand_combs[1:]:
        for i in range(rand_len):
            key[random_pos[i]] = comb[i]
        failed, failed_pos, new_diff = Eve_parity_checks(traces,key)
      
        if new_diff == diff:
            min_list.append(comb)
        elif new_diff < diff:
            min_list = []
            min_list.append(comb)
            diff =  new_diff

    chosen_comb = min_list[0]
    if len(min_list) > 0:
        chosen_comb = random.choice(min_list)
    for i in range(rand_len):
        key[random_pos[i]] = chosen_comb[i]
        

    return diff, min_list  

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


def calc_PE(H_E,N,M):
   # ones_vec = np.ones([N,1])
    #PE = np.matmul(H_E,ones_vec)%2
    PE=np.zeros([1,M])
    for m in range(M):
        p, hw = parity_check(H_E[m,:])
        PE[0,m] = p
    return np.reshape(PE,M)

def calc_PE_from_key(key,traces,M):
    PE = np.zeros([1,M])
    i = 0
    for trace in traces:
        pos = trace[1]
        block = []
        for p in pos:
            block.append(key[p])
        parityEve, hw = parity_check(np.asarray(block))
        PE[0,i] = parityEve
        i = i+1
    return np.reshape(PE,M)





def gen_binary_comb(n):
    bin_strings = []
    def get_bin(n,bs=''):
        if len(bs) == n:
            bin_strings.append(bs)
        else:
            get_bin(n, bs + '0')
            get_bin(n, bs + '1')
    get_bin(n)
    return bin_strings

def get_trace_res1(traces,bit):
    res = []
    for trace in traces:
        parity = trace[0]
        bits_pos = trace[1]
        for i in range(len(bits_pos)):
            if bit == bits_pos[i]:
                res.append(parity[i])
    return res

def get_trace_res2(H_E,Included_matrix,bit):
    res = []
    Inc_vect = Included_matrix[:,bit]
    for check in range(len(Inc_vect)):
        if Inc_vect[check] == 1:
            res.append(H_E[check,bit])
    return res


def plot_confidence(Conf_vector,unconf_pos,misses,sigma,mode):

    invalid_count = 0
    plt.figure()
    bars = plt.bar(np.arange(N),Conf_vector,width=1,edgecolor='black',linewidth=0)
    for i in range(N):
        if misses.__contains__(i):
            bars[i].set_color('black')
            invalid_count = invalid_count + 1

        if unconf_pos.__contains__(i):
            bars[i].set_color('red')
            if misses.__contains__(i):
                bars[i].set_color('green')
                invalid_count = invalid_count - 1

    plt.xlabel("Bit")
    plt.ylabel("Confidence")



    blk_p = mpathces.Patch(color='black',label='Error - Confident')
    blu_p = mpathces.Patch(color='blue', label="Correct - Confident")
    red_p = mpathces.Patch(color='red',label = "Correct - Unconfident")
    grn_p = mpathces.Patch(color='green',label="Error - Unconfident")

    plt.legend(handles=[blk_p,grn_p,blu_p,red_p])
        
    if mode==0:
        plt.title("Bit Confidence for N={0}, $\sigma$={1}".format(N,sigma))
        plt.savefig("results/bit_conf_N_{0}_sig_{1}.png".format(N,sigma))

    if mode==1:
        plt.title("Bit Confidence for N={0}, $\sigma$={1}, with parity offenders".format(N,sigma))
        plt.savefig("results/bit_conf_N_{0}_sig_{1}_offenders.png".format(N,sigma))


    if invalid_count == 0:
        print("The key can be found by flipping {0} bits -> {1} keys".format(len(unconf_pos),'{:.2e}'.format(2**len(unconf_pos))))
    else:
        print("This key will not be found, {0} incorrect bits were assumed confident".format(invalid_count))
    return invalid_count


def partition (list_in, n):
    random.shuffle(list_in)
    return [list_in[i::n] for i in range(n)]

def gather_votes(N,votes,valid_bs,occourances,invalid_count,unconf_pos,majority_key,traces,diff_len,unconf_len,start_split,low_score):

    groups = partition(unconf_pos,int(unconf_len/start_split)) ##random partition
    for group in groups:
        grouplen = len(group)
        bs_list = gen_binary_comb(grouplen)
        for bs in bs_list:
            new_key_temp = majority_key.copy()
            for i in range(len(bs)):
                new_key_temp[group[i]] = int(new_key_temp[group[i]]) ^ int(bs[i])
            failed, failed_pos, new_diff = Eve_parity_checks(traces,new_key_temp)

            if new_diff == 0: ##identifies if comes across a valid solution
                valid_bs.append(bs)
            elif new_diff < diff_len:
                for g in group:
                    if new_key_temp[g] != majority_key[g]:
                        votes[g] = votes[g] + (np.abs(diff_len - new_diff)/occourances[int(g)])
                if low_score == [[],[],[]] or new_diff < low_score[0] :
                    low_score[0] = new_diff
                    low_score[1] = bs
                    low_score[2] = group

            #elif new_diff > diff_len:
             #   for g in group:
              #      if new_key_temp[g] != majority_key[g]:
               #         votes[g] = votes[g] - (np.abs(diff_len - new_diff)/occourances[int(g)])

    return votes, valid_bs, groups, low_score


def recount_votes(N,groups_list,valid_bs,occourances,majority_key,traces,diff_len):
    
    votes = np.zeros(N)
    for groups in groups_list:
        for group in groups:
            grouplen = len(group)
            bs_list = gen_binary_comb(grouplen)

            for bs in bs_list:
                new_key_temp = majority_key.copy()
                for i in range(len(bs)):
                    new_key_temp[group[i]] = int(new_key_temp[group[i]]) ^ int(bs[i])
                failed, failed_pos, new_diff = Eve_parity_checks(traces,new_key_temp)

                if new_diff == 0: ##identifies if comes across a valid solution
                    valid_bs.append(bs)
                elif new_diff < diff_len:
                    for g in group:
                        if new_key_temp[g] != majority_key[g]:
                            votes[g] = votes[g] + (np.abs(diff_len - new_diff)/occourances[int(g)])

                elif new_diff > diff_len:
                    for g in group:
                        if new_key_temp[g] != majority_key[g]:
                            votes[g] = votes[g] - (np.abs(diff_len - new_diff)/occourances[int(g)])
    return votes, valid_bs



def confidence_flip(N,M,majority_key,H_E,Included_matrix,Conf_vector,P_A,start_thresh, misses,traces,sigma):
    #P_E = calc_PE(H_E,N,M)
    P_E = calc_PE_from_key(majority_key,traces,M)
    diff = (P_E + P_A)%2
    diff_len = 0
    for d in diff:
        if d==1:
            diff_len = diff_len + 1

    lowConf = 1
    secondLow = 1

    for c in Conf_vector:
        if c < lowConf:
            lowConf = c
        if c > lowConf and c < secondLow:
            secondLow = c

    unconf_pos = []
    start_thresh = 0.5
   # if (P_E != P_A).any(): ##P_E == P_A is true for some cases of faulty keys
    iterations = 3
    for i in range(N): ##Gather unconfident values
  #      if Conf_vector[i] <= secondLow + 0.01 or Conf_vector[i] <= start_thresh: ##need better system for this
        if Conf_vector[i] < start_thresh:
            unconf_pos.append(i)

    invalid_count = plot_confidence(Conf_vector,unconf_pos,misses,sigma,0)


    # Solve for x using linear algebra techniques
    #################

    Included_matrix_new = Included_matrix.copy()
    P_A_new = P_A.copy()
    dels = 0
    for i in range(N):
        if not unconf_pos.__contains__(i):
            conf_x_vec = majority_key[i]*Included_matrix[:,i]
            delpos = i - dels
            Included_matrix_new = np.delete(Included_matrix_new,delpos,1)
            P_A_new = (P_A_new - conf_x_vec) % 2
            dels = dels + 1

    ###remove any all zero rows (occours when all paity values are known)
    temp = Included_matrix_new.copy()
    temp_PA = P_A_new.copy()
    dels = 0
    for i in range(Included_matrix_new.shape[0]):
        if (Included_matrix_new[i,:] == np.zeros(len(unconf_pos))).all():
            temp = np.delete(temp,i-dels,0)
            temp_PA = np.delete(temp_PA,i-dels,0)
            dels = dels + 1
    Included_matrix_new = temp
    P_A_new = temp_PA


    row_ech,b,row_order = my_solver.row_echelon_form(Included_matrix_new,P_A_new)
    for i in range(row_ech.shape[0]):
        if (row_ech[i,:] == np.zeros(len(unconf_pos))).all():
            row_ech = row_ech[:i,:]
            b = b[:i]
            break
    
    if row_ech.shape[0] == row_ech.shape[1]:
        solution_x = np.linalg.solve(row_ech,b)
    

    #################


        ##This code is still required, but having it out makes developing the flip vote easier
    
    offenders = top_offenders(N,traces,P_A,P_E)
    updated_unconf_pos = []
    offender_thresh = 1
    for o in offenders:
        if o[1] > offender_thresh and unconf_pos.__contains__(o[0]):
            updated_unconf_pos.append(int(o[0]))
        elif Conf_vector[int(o[0])] == 0.1:
            updated_unconf_pos.append(int(o[0]))
    
    invalid_count = plot_confidence(Conf_vector,updated_unconf_pos,misses,sigma,1) 
    unconf_pos = updated_unconf_pos.copy()


    valid_bs = []
    if invalid_count == 0:
        max_unconf = 10
        start_split = 4
        unconf_len = len(unconf_pos)
        if unconf_len >= max_unconf:

            ##count number of occourances in traces
            ##use this to normalise vote values

            occourances = np.zeros(N)
            for trace in traces:
                pos = trace[1]
                for p in pos:
                    occourances[p] = occourances[p] + 1

            majority_key_old = majority_key.copy()
            for it in range(1):
                groupslist = []
                low_score_list = []
                
                votes = np.zeros(N)
                for v in range(10):
                    low_score = [[],[],[]]
                    votes, valid_bs, groups, low_score = gather_votes(N,votes,valid_bs,occourances,invalid_count,unconf_pos,majority_key,traces,diff_len,unconf_len,4,low_score)
                    groupslist.append(groups)
                    low_score_list.append(low_score)
            #   for v in range(5):
            #      votes, valid_bs = gather_votes(N,votes,valid_bs,occourances,invalid_count,unconf_pos,majority_key,traces,diff_len,unconf_len,3)
        #     for v in range(5):
            #        votes, valid_bs = gather_votes(N,votes,valid_bs,occourances,invalid_count,unconf_pos,majority_key,traces,diff_len,unconf_len,4)
        #       for v in range(5):
        #          votes, valid_bs = gather_votes(N,votes,valid_bs,occourances,invalid_count,unconf_pos,majority_key,traces,diff_len,unconf_len,5)

                plt.figure()
                bars = plt.bar(np.arange(N),votes,width=1,edgecolor='black',linewidth=0)
                for i in range(N):
                    if misses.__contains__(i):
                        bars[i].set_color('black')
                        invalid_count = invalid_count + 1

                    if unconf_pos.__contains__(i):
                        bars[i].set_color('red')
                        if misses.__contains__(i):
                            bars[i].set_color('green')
                            invalid_count = invalid_count - 1
                plt.title("Votes for changing bits for N={0}, $\sigma$={1}, iter {2}".format(N,sigma,it))
                plt.savefig("results/flip_votes_N_{0}_sig_{1}_iter{2}.png".format(N,sigma,it))
                

                ##flip max value from votes
                """
                maxpos = np.argmax(votes)
                majority_key[maxpos] = int(majority_key[maxpos]) ^ 1
                unconf_pos.remove(maxpos)
                for groups in groupslist:
                    for group in groups:
                        if group.__contains__(maxpos):
                            group.remove(maxpos) """


                votes, valid_bs = recount_votes(N,groupslist,valid_bs,occourances,majority_key,traces,diff_len)


                plt.figure()
                bars = plt.bar(np.arange(N),votes,width=1,edgecolor='black',linewidth=0)
                for i in range(N):
                    if misses.__contains__(i):
                        bars[i].set_color('black')
                        invalid_count = invalid_count + 1

                    if unconf_pos.__contains__(i):
                        bars[i].set_color('red')
                        if misses.__contains__(i):
                            bars[i].set_color('green')
                            invalid_count = invalid_count - 1
                plt.title("Recount for changing bits for N={0}, $\sigma$={1}, iter {2}".format(N,sigma,it))
                plt.savefig("results/recount_votes_N_{0}_sig_{1}_iter{2}.png".format(N,sigma,it))
                


    ###This code computes the valid keys, commented out for now
    valid_bs = []

    """
    bs_list = gen_binary_comb(len(unconf_pos))
    for bs in bs_list:
        new_key_temp = majority_key.copy()
        for i in range(len(bs)):
            new_key_temp[unconf_pos[i]] = int(new_key_temp[unconf_pos[i]]) ^ int(bs[i])
        failed, failed_pos, no_failed = Eve_parity_checks(traces,new_key_temp)
        if no_failed == 0:
            valid_bs.append(bs)
    """

    return valid_bs, unconf_pos






 
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

        final_key, H_E, Included_matrix, Conf_vector,P_A = recover_secret_majoirty(Eves_traces,N,XOR_noise)
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
        print("\nNoise: {0}".format(XOR_noise))

        #print("Successfully Corrected Errors!")
        ###Run Attack Here

        majority_key, H_E, Included_matrix, Conf_vector,P_A = recover_secret_majoirty(Eves_traces,N,XOR_noise)
        misses = []
        conf_key = majority_key.copy()
        if (majority_key == X).all():
            print("Majority, Successful")
        else:
            for bit in range(N):
                if X[bit] != majority_key[bit]:
                    misses.append(bit)
            valid_bs, unconf_pos = confidence_flip(N,len(Eves_traces),majority_key,H_E,Included_matrix,Conf_vector,P_A,0.25,misses,Eves_traces,XOR_noise) ##misses for debugging
            included_misses = 0
            for m in misses:
                if unconf_pos.__contains__(m):
                    included_misses = included_misses + 1
            if included_misses == len(misses): #and len(valid_bs) > 0:
                print("Valid List, length: {0}".format(len(valid_bs)))
                return len(valid_bs)
            else:
                print("Invalid List, length: {0}".format(len(valid_bs)))
                return 0
            """
            if conf_key[0] == -1:
                print("Full Info, failed")
                return -1

            else:
                misses_new = []
                if (conf_key == X).all():
                    print("Full Info, Success")
                else:
                    for bit in range(N):
                        if X[bit] != conf_key[bit]:
                            misses_new.append(bit)
        
                    success_count = 0
                    for i in range(N):
                        if X[i] == conf_key[i]:
                            success_count = success_count + 1
                    success_rate = success_count/N
                    #print("Eve Success rate: {0}".format(success_rate))
                # print(final_key)
                    return success_rate"""

    

def plot_params(N,QBER,max_iter, err_range):
    tb_only = []
    majority = []
    conf = []

    for err in err_range:
        sr_tb = -1
        sr_maj = -1
        sr_conf = -1
        """
        while sr_tb == -1:
            sr_tb = run_EC_topblock_only(N,QBER,err,max_iter)
        tb_only.append(sr_tb)

        while sr_maj == -1:
            sr_maj = run_EC_majority(N,QBER,err,max_iter)
        majority.append(sr_maj)"""

        #while sr_conf == -1:
        sr_conf = run_EC_confidence_flip(N,QBER,err,max_iter)
       # conf.append(sr_conf)

    return tb_only, majority, conf




err_range = np.arange(0,0.2001,0.01)
N = 128
max_iter = 3
QBER = 0.1
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
