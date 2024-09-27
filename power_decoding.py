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
    power_samples = power_samples[power_samples[:,0].argsort()]
    
    return power_samples


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
    
    confidence = np.abs(prob_0-prob_1)

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
        

##runs error correction algorithm and returns how Eve has perfromed in terms of key recovery
def run_EC_power_flip(N,QBER,sigma,max_iter,cascade_params,initial_only=False):
    #X,Y = cascade.gen_sifted_keys(N,QBER)
    #Y_ec, Eves_traces = cascade.cascade_EC(X,Y,QBER,max_iter,sigma,1)

    algo = cascade_params[0]
    error_type = cascade_params[1]
    Y, X, Y_ec, Eves_traces, remaining_errors = run_reconciliation(algo,N,error_type,QBER,sigma)
    X = np.asarray([*X._bits.values()])
    Y = np.asarray([*Y._bits.values()])
    Y_ec = np.asarray([*Y_ec._bits.values()])

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

        confidence_thresh = 0.95
        sample_min = 2

        if initial_only == False:
            while finished == False:
                remaining_start = (partial_solution == -1).sum()
                parity_list = calculate_intermediate_parities(Eves_traces,partial_solution)
                power_samples = sort_power_samples(Eves_traces,parity_list,partial_solution,Y)


                start_power = -3
                end_power = 4
                power_resolution = 3
                a_base = 0
                b_base = 1
                mid = 0.5
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
                    confidence_thresh = confidence_thresh*5/100
                else:
                    confidence_thresh = 0.95

            cor, inc = check_solution(partial_solution,Y)

        ##if unkowns remain, find by exhaustively matching parity vectors
        #if (inc == 0 and N-cor < 10):


        return cor, inc
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
        return 0, 0

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




def exp1():

    Nrange = [100,200,500,1000]
    QBER_range = np.arange(0.01,0.111,0.01)
    repeats = 5
    max_exhaust = 10
    sigma = 0.5
    algorithms = ["original", "yanetal", "biconf", "option3", "option4", "option7", "option8"]
    QBER_type = "bernoulli"
    cascade_params = [algorithms[0],QBER_type]

    

    cor_avg = []
    inc_avg = []
    for N in Nrange:
        corN = []
        incN = []
        for QBER in QBER_range:
            corQBER = 0
            incQBER = 0
            SR_AVG = 0
            for i in range(repeats):
                cor, inc = run_EC_power_flip(N,QBER,sigma,2,cascade_params,initial_only=False)
                corQBER = corQBER + cor
                incQBER = incQBER + inc
            corN.append((corQBER/repeats)/N)
            incN.append((incQBER/repeats)/N)
        cor_avg.append(corN)
        inc_avg.append(incN)       
    

    with open("results/exp1.csv","w") as csvfile:
        writer = csv.writer(csvfile,delimiter=",")
        writer.writerow(QBER_range)
        for i in range(len(cor_avg)):
            writer.writerow([Nrange[i]] + cascade_params + [sigma,"Correct"] + cor_avg[i])
            writer.writerow([Nrange[i]] + cascade_params + [sigma,"Incorrect"] + inc_avg[i])

def plot_exp1():
    plt.style.use("seaborn-v0_8-darkgrid")
    with open("results/exp1.csv","r") as csvfile:
        reader = csv.reader(csvfile,delimiter=",")
        read_count = 0
        for row in reader:
            if read_count == 0:
                QBER_range = np.asarray(row,dtype="float")
                read_count = read_count + 1
            else:
                N = row[0]
                alg =row[1]
                err_type = row[2]
                sigma = row[3]
                results_type = row[4]
                results = np.asarray(row[5:],dtype="float")
                if results_type == "Correct":
                    plt.plot(QBER_range,results,label="{0}: N={1}, $\sigma$={2}".format(results_type,N,sigma))


    plt.xlabel("QBER")
    plt.ylabel("Full recovery rate")

    plt.legend()
    plt.title("Effect of noise on with \n STM32-F405RGT6 (ARM Cortex-M4)")#attack on general hardware")
    plt.show()
    k=7
   # plt.savefig("results/exp1_large_general.png",bbox_inches = "tight")



def initial_discovery_exp():
    plt.style.use("seaborn-v0_8-darkgrid")
    algorithms = ["original", "yanetal", "biconf", "option3", "option4", "option7", "option8"]
    QBER_type = "bernoulli"
    
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
                revealed.append(revealed_qber/repeats)
            plt.plot(QBER_range,revealed,label="{0} N={1}".format(alg,N))

    plt.legend()
    plt.xlabel("QBER")
    plt.ylabel("Portion of key revealed")
    plt.title("Portion of key revealed after reducing initial system")
    plt.savefig("results/initial_key_bits_reveled_bruno.png",bbox_inches="tight")

#initial_discovery_exp()
#exp1()
plot_exp1()