import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpathces
import random
from decimal import Decimal
import math
from collections import Counter
import linalgtools as my_solver
import cascade_EC as cascade


########################### Eves Tools ########################

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
        parityEve, hw = cascade.parity_check(np.asarray(block))
        if parityAlice != parityEve:
            failed.append(trace)
            failed_pos.extend(pos)


    return failed, failed_pos, len(failed)


def calc_PE(H_E,N,M):
   # ones_vec = np.ones([N,1])
    #PE = np.matmul(H_E,ones_vec)%2
    PE=np.zeros([1,M])
    for m in range(M):
        p, hw = cascade.parity_check(H_E[m,:])
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
        parityEve, hw = cascade.parity_check(np.asarray(block))
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



def confidence_flip(N,M,majority_key,H_E,Included_matrix,Conf_vector,P_A,start_thresh, misses,traces,sigma, Y):
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

    print(Included_matrix_new)
    row_ech,b,row_order = my_solver.row_echelon_form(Included_matrix_new,P_A_new)
    for i in range(row_ech.shape[0]):
        if (row_ech[i,:] == np.zeros(len(unconf_pos))).all():
            row_ech = row_ech[:i,:]
            b = b[:i]
            break

    correct_soln = []
    for i in unconf_pos:
        correct_soln.append(int(Y[i]))

    calc_soln = my_solver.find_solution(row_ech,b)

    for soln in calc_soln:
        solnl = soln.tolist()
        if (solnl == correct_soln):
            print("Correct solution contained in solution space of length: {0}".format(calc_soln.shape[0]))
            break
    

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
 


def run_EC_majority(N,QBER,XOR_noise,max_iter):
    X,Y = cascade.gen_sifted_keys(N,QBER)
    Y_ec, Eves_traces = cascade.cascade_EC(X,Y,QBER,max_iter,XOR_noise,1)
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
    X,Y = cascade.gen_sifted_keys(N,QBER)
    Y_ec, Eves_traces = cascade.cascade_EC(X,Y,QBER,max_iter,XOR_noise,0)
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
    X,Y = cascade.gen_sifted_keys(N,QBER)
    Y_ec, Eves_traces = cascade.cascade_EC(X,Y,QBER,max_iter,XOR_noise,1)
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
            valid_bs, unconf_pos = confidence_flip(N,len(Eves_traces),majority_key,H_E,Included_matrix,Conf_vector,P_A,0.25,misses,Eves_traces,XOR_noise,Y) ##misses for debugging
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
N = 300
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
