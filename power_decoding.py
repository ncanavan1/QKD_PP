import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpathces
import random
from decimal import Decimal
import math
from collections import Counter
import linalgtools as my_solver
from cascade.run_EC import run_reconciliation
import noise_modelling
import csv
#import LinAlg


########################### Eves Tools ########################


def gen_system_from_trace(traces,N):

    M = len(traces)
    P_A = np.zeros(M)
    Included_matrix = np.zeros([M,N]) ##linear system matrix
    partial_key = np.ones(N)*-1

    check = 0
    for trace in traces:
        #inter_parity = trace[0]
        pos = trace[1]
        out_parity = trace[2]
        P_A[check] = out_parity
        if pos.size == 1:
            partial_key[pos[0]] = out_parity
        for i in reversed(range(0,pos.size)): 
            Included_matrix[check,pos[i]] = 1
        check = check + 1
    return Included_matrix, P_A, partial_key


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


def power_recovery_only(N, Eves_traces, a_base, b_base, majority=False):
    recovered_key = [[] for _ in range(N)]

    for trace in Eves_traces:
        powers = trace[0]
        bit_pos = trace[1]

        cur_parity = 0
        for i in range(len(bit_pos)):
            if np.abs(powers[i] - a_base) < np.abs(powers[i] - b_base):
                bit = int(cur_parity) ^ 0
                recovered_key[bit_pos[i]].append(bit)
                cur_parity = 0
            else:
                bit = int(cur_parity) ^ 1
                recovered_key[bit_pos[i]].append(bit)
                cur_parity = 1
    
    recovered_key_return = np.zeros(N)
    if majority == False:
        for i in range(N):
            recovered_key_return[i] = recovered_key[i][0]

    else:
        for i in range(N):
            avg = sum(recovered_key[i])/len(recovered_key[i])
            if avg < 0.5:
                recovered_key_return[i] = 0
            elif avg > 0.5:
                recovered_key_return[i] = 1
            else:
                recovered_key_return[i] = np.random.choice([0,1])
    
    return recovered_key_return


    





def calculate_intermediate_parities(Traces, partial_solution):

    parity_list = []

    ##calculate from behind using p_i-1 = 0
    for trace in Traces:
        bit_pos = trace[1]
        parities = []
        parities.append(0)
        unkown_found = False
        count = 0
        for i in bit_pos:

            if partial_solution[i] == -1:
                unkown_found = True

            if unkown_found == True:
                parities.append(-1)

            else:
                parities.append(int(parities[count] ^ int(partial_solution[i])))
            count = count + 1
        parities.pop()
        parities.append(trace[2])
        parity_list.append(parities)
    
    ##calculate from end using known output parity
    
    current = 0
    for trace in Traces:
        bit_pos = trace[1]
        count = len(parity_list[current])-1

        for i in reversed(bit_pos):

            if partial_solution[i] != -1:
                parity_list[current][count-1] = parity_list[current][count] ^ int(partial_solution[i])
                count = count - 1

            else:
                break

        current = current + 1
    
    return parity_list
    

##makes a list of power samples for parities involving each bit and a known auxilary parity value
def sort_power_samples(Traces,parity_list,partial_solution,Y):

    power_samples = []
    empty = True
    eqn_count = 0
    for trace in Traces:
        power = trace[0]
        bit_pos = trace[1]
        parities = parity_list[eqn_count]

        parity_count  = 0
        for i in bit_pos:
            ##x_i is unknown and preceding parity is known
            if(partial_solution[i] == -1 and parities[parity_count] != -1):
            #if(parities[parity_count] != -1):
                power_samples.append([int(i),int(parities[parity_count]),power[parity_count]])
            parity_count = parity_count + 1

        
        parity_count = len(parities)-1
        for i in reversed(bit_pos):
            if(partial_solution[i] == -1 and parities[parity_count] != -1 and parities[parity_count-1] == -1):
                power_samples.append([int(i),int(parities[parity_count]),power[int(parity_count-2)]])
            parity_count = parity_count-1
        
        eqn_count = eqn_count + 1

        
    
        

    power_samples = np.asarray(power_samples)
   # for p in power_samples:
    #    print("{0}, {1}".format(p,Y[int(p[0])]))
    try:
        power_samples = power_samples[power_samples[:,0].argsort()]
        empty = False
    except:
        power_samples = np.asarray([])
        empty = True

    return power_samples, empty


def mirror_power(y,mid):
    mid_diff = np.abs(y-mid)
    if y < mid:
        return mid + mid_diff
    else:
        return mid - mid_diff     
        
def guess_bit_with_samples(power_samples,mid,sigma,power_resolution,start_power,end_power,a_base,b_base,pos,Y):
    
    S = []
    for s in power_samples:
        if s[1] == 1:
            S.append(mirror_power(s[2],mid))
        else:
            S.append(s[2])


    x = np.arange(start_power,end_power,10**(-power_resolution))
    x = np.round(x,decimals=power_resolution)
    cdf_a, cdf_b = noise_modelling.calc_cdfs(x,a_base,b_base,sigma)

    prob_0, prob_1 = noise_modelling.prob_given_S(S,x,cdf_a,cdf_b,power_resolution)

    if power_samples[0,1] == 1:
        tmp = prob_0
        prob_0 = prob_1
        prob_1 = tmp

    ##xor is likey XOR->0
    if prob_0 > prob_1:
        guess = 0 ^ int(power_samples[0,1])

    if prob_0 < prob_1:
        guess = 1 ^ int(power_samples[0,1])
    else:
        guess = -1
        confidence = 0
    
    try:
        confidence = np.abs(prob_0-prob_1)
    except:
        k=7
    return guess, confidence

def check_solution(partial_solution,Y):
    cor = 0
    inc = 0
    for i in range(len(partial_solution)):
        if partial_solution[i] != -1:
            # print("{0} {1}".format(partial_solution[i],Y[i]))  
            if partial_solution[i] == Y[i]: 
                cor = cor + 1
            else:
                inc = inc + 1
                print(i)

    print("{0}% of key correctly discovered, {1}% incorrect".format((cor*100)/len(Y), (inc*100)/len(Y))) 
    return cor, inc 



def linear_system_resolve(N,Included_matrix,P_A,partial_solution,Y):

    debugging = True
     # Solve for x using linear algebra techniques
    #################

    Included_matrix_new = Included_matrix.copy()
    P_A_new = P_A.copy()
    dels = 0

    bit_pos = []

    ###build a new matrix that includes only the variables (bits) still unconfident
    for i in range(N):
        if partial_solution[i] != -1:
            conf_x_vec = partial_solution[i]*Included_matrix[:,i]
            delpos = i - dels
            Included_matrix_new = np.delete(Included_matrix_new,delpos,1)
            P_A_new = (P_A_new - conf_x_vec) % 2
            dels = dels + 1
        else:
            bit_pos.append(i)

    ###remove any all zero rows (occours when all bits are confident in the row)
    temp = Included_matrix_new.copy()
    temp_PA = P_A_new.copy()
    dels = 0
    for i in range(Included_matrix_new.shape[0]):
        if (Included_matrix_new[i,:].max() == 0):
            temp = np.delete(temp,i-dels,0)
            temp_PA = np.delete(temp_PA,i-dels,0)
            dels = dels + 1
    Included_matrix_new = temp
    P_A_new = temp_PA

    ##convert linear system to reduced row echelon form
    print("Converting to RREF...")
    row_ech,b,row_order = my_solver.row_echelon_form(Included_matrix_new,P_A_new)

  
    ##Get solution space for linear system
    print("Solving solution of size {0}x{1}...".format(row_ech.shape[0],row_ech.shape[1]))

    newfound = 0
    for row in range(row_ech.shape[0]):
        if np.sum(row_ech[row]) == 1:
            bit = np.where(row_ech[row] == 1)[0][0]
            partial_solution[bit_pos[bit]] = P_A_new[row]
            newfound = newfound+1
           # if(Y[bit_pos[bit]] == partial_solution[bit_pos[bit]]):
            #    print("x{0} correct".format(row))
            #else:
             #   print("x{0} incorrect".format(row))
    print("Found {0} new bit values".format(newfound))
    
    return partial_solution

def compare_output_parities(partial_key, Eves_traces):
    Err_list = []
    trace_ctr = 0
    for trace in Eves_traces:
        bit_pos = trace[1]
        true_parity = trace[2]
        completed = True
        calc_parity = 0
        for bit in bit_pos:
            if (partial_key[bit] == -1):
                completed = False
                break
            else:
                calc_parity = calc_parity ^ int(partial_key[bit])
        if completed == True:
            if calc_parity != true_parity:
                Err_list.append(trace_ctr)
        trace_ctr = trace_ctr+1
    return Err_list
        
            




##runs error correction algorithm and returns how Eve has perfromed in terms of key recovery
def run_EC_power_flip(N,QBER,sigma,max_iter,cascade_params, delta_p = 1, initial_only=False):
    #X,Y = cascade.gen_sifted_keys(N,QBER)
    #Y_ec, Eves_traces = cascade.cascade_EC(X,Y,QBER,max_iter,sigma,1)

    algo = cascade_params[0]
    error_type = cascade_params[1]
    

    a_base = 0
    b_base = delta_p
    mid = (a_base+b_base)/2
    start_power = np.floor(a_base - 5*sigma)
    end_power = np.ceil(b_base + 5*sigma)
    power_resolution = 3
    power_params = [sigma, a_base, b_base]


    Y, X, Y_ec, Eves_traces, remaining_errors = run_reconciliation(algo,N,error_type,QBER,power_params)
   # P_A = np.asarray([trace[2] for trace in Eves_traces],dtype=int)
    X = np.asarray([*X._bits.values()])
    Y = np.asarray([*Y._bits.values()])
    Y_ec = np.asarray([*Y_ec._bits.values()])


 

    naive_key = power_recovery_only(N,Eves_traces,a_base,b_base,majority=True)
    power_only_cor, power_only_inc = check_solution(naive_key,Y)

    if remaining_errors == 0:
        print("\nN:  {0}, QBER: {1}, Noise: {2}".format(N,QBER,sigma))

        #print("Successfully Corrected Errors!")
        ###Run Attack Here

        Included_matrix, P_A, partial_solution = gen_system_from_trace(Eves_traces,N)

        initial_cor, initial_inc = check_solution(partial_solution,Y)

        ####solve system in completely unknown 
        #row_ech, b, row_order = my_solver.row_echelon_form(Included_matrix.copy(),P_A.copy())
        #partial_solution = my_solver.read_off_basic_variables(row_ech,b) 
        
        partial_solution = linear_system_resolve(N,Included_matrix,P_A,partial_solution,Y)

        cor, inc = check_solution(partial_solution,Y)

        remaining_start = (partial_solution == -1).sum()
        finished = False
        confidence_thresh_start = 0.99
        confidence_thresh = 0.99
        sample_min = 2

        if initial_only == False:
            while finished == False:
                remaining_start = (partial_solution == -1).sum()
                parity_list = calculate_intermediate_parities(Eves_traces,partial_solution)
                power_samples, empty = sort_power_samples(Eves_traces,parity_list,partial_solution,Y)
                if empty == True:
                    break



               # start_power = -10.675
                #end_power = -2.675
                #a_base = -6.9
                #b_base = -6.45
                #mid = -6.675

                confidence_ratings = np.zeros(N)

                s_curr = power_samples[0,0]
                s_curr_start = 0
                s_curr_end = 0
                pos = 0
                for s in power_samples:
                    if s[0] != s_curr:
                        s_curr_end = pos

                        if  s_curr_end-s_curr_start >= sample_min:
                            guess, confidence = guess_bit_with_samples(power_samples[s_curr_start:s_curr_end,:],mid,sigma,power_resolution,start_power,end_power,a_base,b_base,s_curr,Y)
                            if confidence > confidence_thresh:
                                partial_solution[int(s_curr)] = guess
                                if guess != Y[int(s_curr)]:
                                    err_list = compare_output_parities(partial_solution,Eves_traces)
                                    guess_bit_with_samples(power_samples[s_curr_start:s_curr_end,:],mid,sigma,power_resolution,start_power,end_power,a_base,b_base,s_curr,Y)

                                confidence_ratings[int(s_curr)] = confidence
                    
                        s_curr_start = pos
                        s_curr = s[0]
                    pos = pos+1

                check_solution(partial_solution,Y)

                partial_solution = linear_system_resolve(N,Included_matrix,P_A,partial_solution,Y)
                check_solution(partial_solution,Y)

                remaining_end = (partial_solution == -1).sum()
                if (remaining_end == remaining_start or remaining_end == 0):
                    if confidence_thresh < 0.95 or remaining_end == 0:
                        finished = True
                    else:
                        confidence_thresh = confidence_thresh*(0.99)
                else:
                    confidence_thresh = confidence_thresh_start

            cor, inc = check_solution(partial_solution,Y)

        ##if unkowns remain, find by exhaustively matching parity vectors
        #if (inc == 0 and N-cor < 10):

      #  arr = top_offenders(partial_solution,Eves_traces,P_A)
        return cor, inc, power_only_cor, power_only_inc
        """
        linalg2  = True
        if linalg2 == False:
            if len(random_choice_pos) > 0:
                print("Key contains {0} random choices".format(len(random_choice_pos)))
                
                ####ensure that we don't try too large spaces
                if len(random_choice_pos) <= 8:
                    diff, min_list = test_random_choices(partial_key,random_choice_pos,Eves_trac es)
                else:
                    for pos in random_choice_pos:
                        partial_solution[pos] = random.choice([0,1])
                        """
        



    
            
    else:
        return 0, 0, 0, 0


def calc_PE_from_key(key,traces,M):
    PE = np.zeros(M,dtype=int)
    for i in range(M):
        bit_pos = traces[i][1]
        cur_parity = 0
        for j in range(len(bit_pos)):
            if key[bit_pos[j]] == -1:
                cur_parity = -1
                break
            else:
                cur_parity = cur_parity ^ int(key[bit_pos[j]])
        PE[i] = cur_parity
    return PE


def top_offenders(recovered_key,traces,P_A):
    num_checks = len(traces)
    P_E = calc_PE_from_key(recovered_key,traces,num_checks)
    offence_count = np.zeros(N)
    pos_list = np.arange(0,N,1)
    for i in range(num_checks):
        if P_E[i] != -1:
            if P_A[i] != P_E[i]:
                offence_trace = traces[i][1]
                for pos in offence_trace:
                    offence_count[pos] = offence_count[pos] + 1
    combined =  np.vstack((pos_list,offence_count)).T
    combined2 = np.argsort(-combined[:,1])
    combined = combined[combined2]
    return combined

##this experiment aims to show the percentage of key successfully recovered by the attack for various N and noise
def exp2():
    Nrange = [100,200,500]#,500,1000]
    noise_range = np.arange(0.01,0.81,0.05)
    repeats = 5
    max_exhaust = 10
    QBER = 0.1


    plt.figure()

    sr = []
    for N in Nrange:
        srN = []
        for sigma in noise_range:
            SR_AVG = 0
            for i in range(repeats):
                cor, inc = run_EC_power_flip(N,QBER,sigma,3)
                if inc == 0 and (N-cor < max_exhaust):
                    SR_AVG = SR_AVG + 1
            SR_AVG = SR_AVG/repeats       
            srN.append(SR_AVG)     
        sr.append(srN)

    plt.plot(noise_range,sr[0],label="N={0}".format(Nrange[0]))
    plt.plot(noise_range,sr[1],label="N={0}".format(Nrange[1]))
    plt.plot(noise_range,sr[2],label="N={0}".format(Nrange[2]))
    #plt.plot(noise_range,sr[3],label="N={0}".format(Nrange[3]))


    plt.xlabel("Trace noise")
    plt.ylabel("Full recovery rate")

    plt.legend()
    plt.title("Effect of noise on attack on general hardware")#with \n STM32-F405RGT6 (ARM Cortex-M4)")
    plt.savefig("results/exp2_large_general.png",bbox_inches = "tight")






def exp_delp_to_sig_ratio():
    N_range = [250]
    QBER_range = np.array([0.15])#,0.15])
    repeats = 3
    delta_p_range = np.array([0.5,1,2])
    del_p_noise_ratios = np.arange(0.1,5.1,0.5)
    algorithms = ["original"]#, "yanetal", "biconf", "option3", "option4", "option7", "option8"]
    QBER_type = "bernoulli"

    with open("results/vary_dpn_ratio_exp.csv","w") as csvfile:
        writer = csv.writer(csvfile,delimiter=",")
        header = ["N", "QBER", "delta_p", "dpn_ratio",  "Algorithm", "Error Type" , "Result Type", "Result"]
        writer.writerow(header)
        for N in N_range:
            for alg in algorithms:
                cascade_params = [alg,QBER_type]
                for QBER in QBER_range:
                    for delta_p in delta_p_range:
                        for dpn in del_p_noise_ratios:
                            sigma = delta_p/dpn
                            for _ in range(repeats):
                                print("N={0}, QBER = {1}, sigma={2}, del_p = {3}".format(N,QBER,sigma,delta_p))
                                cor, inc, power_only_cor, power_only_inc = run_EC_power_flip(N,QBER,sigma,2,cascade_params, delta_p = delta_p, initial_only=False)
                                writer.writerow([N, QBER, delta_p, dpn] + cascade_params + ["Correct", cor])
                                writer.writerow([N, QBER, delta_p, dpn] + cascade_params + ["Incorrect", inc])
                                writer.writerow([N, QBER, delta_p, dpn] + cascade_params + ["Correct Power Only", power_only_cor])
                                writer.writerow([N, QBER, delta_p, dpn] + cascade_params + ["Incorrect Power Only", power_only_inc])




def plot_del_p_to_sig_ratio():
    plt.style.use("seaborn-v0_8-darkgrid")
    df = pd.read_csv("results/vary_dpn_ratio_exp.csv")
    print(df)

    ##averaged results
    results = df.groupby(['N', 'Result Type', 'delta_p', 'dpn_ratio', 'QBER'])["Result"].mean()
    results = results.reset_index()
    N_range = results['N'].unique()
    Result_types = results['Result Type'].unique()
    dpn_ratio_range = results['dpn_ratio'].unique()
    QBER_range = results['QBER'].unique()
    delta_p_range = results['delta_p'].unique()


    for N in N_range:
        for result_type in Result_types:
            for QBER in QBER_range:
                for delta_p in delta_p_range:
                    filtered = results[results["Result Type"] == result_type]
                    filtered = filtered[filtered["QBER"] == QBER]
                    filtered = filtered[filtered["N"] == N]
                    filtered = filtered[filtered["delta_p"] == delta_p]

                    if result_type == "Correct" or result_type == "Incorrect":
                        if QBER == 0.15:
                            plt.plot(filtered["dpn_ratio"],filtered["Result"]/N,label="{0}: QBER={1:.2f},$\Delta_p$={2:.2f}".format(result_type, QBER,delta_p))
    plt.legend(bbox_to_anchor=(1.05,1.2))
  #  plt.show()
    plt.xlabel("$\Delta_p$/$\sigma$",fontsize=16)
    plt.ylabel("Revealed Rate",fontsize=16)
    plt.savefig("results/vary_dpn_ratio.png",bbox_inches="tight")
    k=7



def exp_vary_del_p():
    N_range = [250]
    QBER_range = np.array([0.05,0.15])
    repeats = 5
    sigma_range = np.arange(0.2,2.1,0.2)
    algorithms = ["original"]#, "yanetal", "biconf", "option3", "option4", "option7", "option8"]
    QBER_type = "bernoulli"
    delta_p_range = np.arange(0.2,2.1,0.2)

    with open("results/vary_delta_p_exp.csv","w") as csvfile:
        writer = csv.writer(csvfile,delimiter=",")
        header = ["N", "QBER", "delta_p", "Sigma",  "Algorithm", "Error Type" , "Result Type", "Result"]
        writer.writerow(header)

        for N in N_range:
            for alg in algorithms:
                cascade_params = [alg,QBER_type]
                for QBER in QBER_range:
                    for delta_p in delta_p_range:
                        for sigma in sigma_range:
                            for _ in range(repeats):
                                print("N={0}, QBER = {1}, sigma={2}, del_p = {3}".format(N,QBER,sigma,delta_p))
                                cor, inc, power_only_cor, power_only_inc = run_EC_power_flip(N,QBER,sigma,2,cascade_params, delta_p = delta_p, initial_only=False)
                                writer.writerow([N, QBER, delta_p, sigma] + cascade_params + ["Correct", cor])
                                writer.writerow([N, QBER, delta_p, sigma] + cascade_params + ["Incorrect", inc])
                                writer.writerow([N, QBER, delta_p, sigma] + cascade_params + ["Correct Power Only", power_only_cor])
                                writer.writerow([N, QBER, delta_p, sigma] + cascade_params + ["Incorrect Power Only", power_only_inc])
                            
    
                            


def exp_vary_QBER():

    Nrange = [250]#,500]#,1000,2000]
    #QBER_range = np.arange(0.01,0.111,0.01)
    QBER_range = np.array([0.11])
    repeats = 5
    max_exhaust = 10
    sigma_range = np.arange(0.1,1.1,0.1)
    algorithms = ["original"]#, "yanetal", "biconf", "option3", "option4", "option7", "option8"]
    QBER_type = "bernoulli"

    with open("results/tester_majority.csv","w") as csvfile:
        writer = csv.writer(csvfile,delimiter=",")
        header = ["N","Algorithm","Error T ype","Sigma","Result Type"] + QBER_range.tolist()
        writer.writerow(header)
        
        for N in Nrange:
            for alg in algorithms:
                cascade_params = [alg,QBER_type]
                for i in range(repeats):
                    for sigma in sigma_range:
                        correct = []
                        incorrect = []
                        correct_power_only = []
                        incorrect_power_only = []
                        for QBER in QBER_range:
                            cor, inc, power_only_cor, power_only_inc = run_EC_power_flip(N,QBER,sigma,2,cascade_params,initial_only=False)
                            correct.append(cor/N)
                            incorrect.append(inc/N)
                            correct_power_only.append(power_only_cor/N)
                            incorrect_power_only.append(power_only_inc/N)
                        writer.writerow([N] + cascade_params + [sigma,"Correct"] + correct)
                        writer.writerow([N] + cascade_params + [sigma,"Incorrect"] + incorrect)
                        writer.writerow([N] + cascade_params + [sigma,"Correct Power Only"] + correct_power_only)
                        writer.writerow([N] + cascade_params + [sigma,"Incorrect Power Only"] + incorrect_power_only)
      
import pandas as pd


def plot_exp_vary_del_p():
    plt.style.use("seaborn-v0_8-darkgrid")
    df = pd.read_csv("results/vary_delta_p_exp.csv")
    print(df)

    ##averaged results
    results = df.groupby(['N', 'Result Type', 'delta_p', 'Sigma', 'QBER'])["Result"].mean()
    results = results.reset_index()
    N_range = results['N'].unique()
    Result_types = results['Result Type'].unique()
    Sigma_r = results['Sigma'].unique()
    QBER_range = results['QBER'].unique()
    delta_p_range = results['delta_p'].unique()


    for N in N_range:
        for result_type in Result_types:
            for QBER in QBER_range:
                for delta_p in delta_p_range:
                    filtered = results[results["Result Type"] == result_type]
                    filtered = filtered[filtered["QBER"] == QBER]
                    filtered = filtered[filtered["N"] == N]
                    filtered = filtered[filtered["delta_p"] == delta_p]

                    if result_type == "Correct" or result_type == "Incorrect":
                        if QBER == 0.15:
                            plt.plot(filtered["Sigma"],filtered["Result"]/N,label="{0}: QBER={1:.2f},$\Delta_p$={2:.2f}".format(result_type, QBER,delta_p))
    plt.legend(bbox_to_anchor=(1.05,1.2))
  #  plt.show()
    plt.xlabel("Noise, $\sigma$",fontsize=16)
    plt.ylabel("Revealed Rate",fontsize=16)
    plt.savefig("results/vary_delta_exp.png",bbox_inches="tight")
    k=7





def plot_exp_vary_sigma():
    plt.style.use("seaborn-v0_8-darkgrid")
    df = pd.read_csv("results/exp_vary_QBER_test_conservative_fixed.csv")

    columns_to_average = [str(col) for col in df.columns[5:] if 0.01 <= float(col) <= 0.11]


    # Step 3: Group the DataFrame by 'N' and calculate the mean for the selected columns
    average_results = df.groupby(['N','Result Type', 'Sigma'])[columns_to_average].mean()
    average_results = average_results.reset_index()
    print(average_results.to_string())
    N_range = average_results['N'].unique()
    Result_types = average_results['Result Type'].unique()
    Sigma_r = average_results['Sigma'].unique()
    QBER_range = [0.01,0.02,0.03,0.04,0.05,0.060000000000000005,0.06999999999999999,0.08,0.09,0.09999999999999999,0.11]

    plt.subplot(2,1,2)
    for N in N_range:
        for result_type in Result_types:
            # Filter data for the current N and Result Type
            filtered_data = average_results[average_results['Result Type'] == result_type]
            filtered_data = filtered_data[filtered_data['N'] == N]
            #print(filtered_data.to_string())
            #k=7
            if (N==250) and result_type == "Correct":
                for QBER in QBER_range:
                    plt.plot(Sigma_r,filtered_data[str(QBER)],label="QBER={0:.2f}".format(QBER))
    
    plt.legend(bbox_to_anchor=(1.05,1.2))
    plt.xlabel("Noise, $\sigma$",fontsize=12)
    plt.ylabel("Correct Key Portion",fontsize=12)

    plt.subplot(2,1,1)
    for N in N_range:
        for result_type in Result_types:
            # Filter data for the current N and Result Type
            filtered_data = average_results[average_results['Result Type'] == result_type]
            filtered_data = filtered_data[filtered_data['N'] == N]
            #print(filtered_data.to_string())
            #k=7
            if (N==250) and result_type == "Incorrect":
                for QBER in QBER_range:
                    plt.plot(Sigma_r,filtered_data[str(QBER)],label="QBER={0:.2f}".format(QBER))


    plt.xlabel("Noise, $\sigma$",fontsize=12)
    plt.ylabel("Incorrect Key Portion",fontsize=12)
    plt.tight_layout()
    plt.savefig("results/poster_vary_sigma.png",dpi=500)
    plt.show()
        
    k=7

def plot_exp_vary_QBER():
    plt.style.use("seaborn-v0_8-darkgrid")
    N_prev = 0
    df = pd.read_csv("results/exp_vary_QBER_test_conservative_fixed.csv")

    columns_to_average = [str(col) for col in df.columns[5:] if 0.01 <= float(col) <= 0.11]

# Step 3: Group the DataFrame by 'N' and calculate the mean for the selected columns
    average_results = df.groupby(['N','Result Type','Sigma'])[columns_to_average].mean()
    print(average_results.to_string())

    for (n, result_type, sigma), row in average_results.iterrows():   
        if (sigma == 0.1 or sigma == 0.9): 
            plt.plot(np.asarray(columns_to_average,dtype=float), row.values, label=f'N={n}, $\sigma$={sigma}, {result_type}')

    plt.xlabel("QBER")
    plt.ylabel("Key recovery rate")

    plt.legend()
    plt.title("Effect of noise on with \n STM32-F405RGT6 (ARM Cortex-M4)")#attack on general hardware") 
    plt.savefig("results/exp_vary_QBER_large_general_test.png",bbox_inches = "tight",dpi=500)
    plt.show()
    k=7





def plot_diff_alg():
    plt.style.use("seaborn-v0_8-darkgrid")
    N_prev = 0
    df = pd.read_csv("results/tester_majority.csv")

    columns_to_average = [str(col) for col in df.columns[5:] if 0.01 <= float(col) <= 0.11]

# Step 3: Group the DataFrame by 'N' and calculate the mean for the selected columns
    average_results = df.groupby(['N','Result Type','Sigma', 'Algorithm'])[columns_to_average].mean()
    print(average_results.to_string())
    average_results = average_results.reset_index()
    N_range = average_results['N'].unique()
    Result_types = average_results['Result Type'].unique()
    Algs = average_results['Algorithm'].unique()
    Sigma_r = average_results['Sigma'].unique()

    
    for N in N_range:
        for result_type in Result_types:
            for alg in Algs:
                # Filter data for the current N and Result Type
                filtered_data = average_results[average_results['Result Type'] == result_type]
                filtered_data = filtered_data[filtered_data['N'] == N]
                filtered_data = filtered_data[filtered_data['Algorithm'] == alg]
                #print(filtered_data.to_string())
                #k=7
                if (N==250):
                        plt.plot(Sigma_r,filtered_data['0.11'],label="{0} - {1}".format(alg,result_type))

    plt.xlabel("Noise, $\sigma$",fontsize=12)
    plt.ylabel("Key recovery rate",fontsize=12)

    plt.legend()
   # plt.title("Effect of noise on with \n STM32-F405RGT6 (ARM Cortex-M4)")#attack on general hardware") 
    plt.savefig("results/tester.png",bbox_inches = "tight",dpi=500)
    plt.show()
    k=7

def initial_discovery_exp():
    plt.style.use("seaborn-v0_8-darkgrid")
    algorithms = ["original", "yanetal", "biconf", "option3", "option4", "option7", "option8"]
    QBER_type = "bernoulli"

    with open("results/exp_initial_disc.csv","w") as csvfile:
        writer = csv.writer(csvfile,delimiter=",")
    
        for alg in algorithms:
            #N_range = [100,200]#,500]
            N_range = [250,500]
            sigma = 0.1
            repeats = 5 
            for N in N_range:
                revealed = []
                QBER_range = np.arange(0.01,0.111,0.01)
                for QBER in QBER_range:
                    revealed_qber = 0
                    for r in range(repeats):
                        cascade_params = [alg,QBER_type]
                        cor, inc = run_EC_power_flip(N,QBER,sigma,2,cascade_params,initial_only=True)
                        revealed_qber = revealed_qber + (cor/N)
                        writer.writerow([N,QBER,sigma,revealed_qber])
                    revealed.append(revealed_qber/repeats)
                plt.plot(QBER_range,revealed,label="{0} N={1}".format(alg,N))

    plt.legend(bbox_to_anchor=(1.05,1))
    plt.xlabel("QBER")
    plt.ylabel("Portion of key revealed")
   # plt.title("Portion of key revealed after reducing initial system")
    plt.tight_layout()

    plt.savefig("results/initial_key_bits_reveled_bruno.png",dpi=500)
    plt.show()
    k=7


def exp_reconcilliation_efficiency():
    N_range = [250,500,1000]
    QBER_range = np.arange(0,0.21,0.01)
    sigma = 0.01
    repeats = 10
    for N in N_range:
        f_ec = []
        for QBER in QBER_range:
            parities_revealed = []
            for _ in range(repeats):
                Y, X, Y_ec, Eves_traces, remaining_errors = run_reconciliation("original",N,"bernoulli",QBER,sigma)
                parities_revealed.append(len(Eves_traces))
            m = sum(parities_revealed)/len(parities_revealed)
    
            h_eps = -QBER*np.log2(QBER) - (1 - QBER)*np.log2(1 - QBER)
            R = 1 - (m/N)
            f_ec.append((1-R)/h_eps)
        plt.plot(QBER_range,f_ec,label="N={0}".format(N))
    
    plt.legend()
    plt.show()
    k=7





#initial_discovery_exp()
#exp_vary_QBER()
#plot_diff_alg()
#plot_exp_vary_QBER()
#plot_exp_vary_sigma()
#exp_reconcilliation_efficiency()
#exp_vary_del_p()
#plot_exp_vary_del_p()
#exp_delp_to_sig_ratio()
plot_del_p_to_sig_ratio()
#N,QBER,sigma,cascade_params = 250,0.15,0.8,["original","bernoulli"]
#run_EC_power_flip(N,QBER,sigma,2,cascade_params,initial_only=False)
