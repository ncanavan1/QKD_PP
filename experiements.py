import numpy as np
import matplotlib.pyplot as plt
import csv
import power_decoding
import pandas as pd
import sys
import json


###test_general_experiement setup
def exp_general(exp_params):

    Nrange = exp_params["N"]
    repeats = exp_params["repeats"]
    noise_range = np.arange(exp_params["sigma_start"],exp_params["sigma_end"],exp_params["sigma_resolution"])
    QBER_range = np.arange(exp_params["QBER_start"],exp_params["QBER_end"],exp_params["QBER_resolution"])
    delta_p_range = np.arange(exp_params["delta_p_start"],exp_params["delta_p_end"],exp_params["delta_p_resolution"])
    algorithms = exp_params["algorithms"]
    QBER_type = exp_params["QBER_type"]
    min_confidence = exp_params["min_confidence"]
    samples_min = exp_params["samples_min"]
    output_file = exp_params["output_file"]

    with open(output_file,"w") as csvfile:
        writer = csv.writer(csvfile,delimiter=",")
        header = ["N", "QBER", "delta_p", "noise", "dpn_ratio",  "Algorithm", "Error Type" , "min_confidence", "samples_min", "Result Type", "Result"]
        writer.writerow(header)
        for alg in algorithms:
            for err_type in QBER_type:
                for N in Nrange:
                    for delta_p in delta_p_range:
                        for QBER in QBER_range:
                            for sigma in noise_range:
                                for _ in range(repeats):
                                    cor, inc, power_only_cor, power_only_inc = power_decoding.run_EC_power_flip(N,QBER,sigma,[alg,err_type],min_confidence,samples_min,delta_p = delta_p,initial_only=False)
                                    dpn = delta_p/sigma
                                    writer.writerow([N, QBER, delta_p, sigma, dpn, alg, err_type, min_confidence, samples_min, "Correct", cor])
                                    writer.writerow([N, QBER, delta_p, sigma, dpn, alg, err_type, min_confidence, samples_min, "Incorrect", inc])
                                    writer.writerow([N, QBER, delta_p, sigma, dpn, alg, err_type, min_confidence, samples_min, "Correct Power Only", power_only_cor])
                                    writer.writerow([N, QBER, delta_p, sigma, dpn, alg, err_type, min_confidence, samples_min, "Incorrect Power Only", power_only_inc])




   


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
                                cor, inc, power_only_cor, power_only_inc = power_decoding.run_EC_power_flip(N,QBER,sigma,2,cascade_params, delta_p = delta_p, initial_only=False)
                                writer.writerow([N, QBER, delta_p, dpn] + cascade_params + ["Correct", cor])
                                writer.writerow([N, QBER, delta_p, dpn] + cascade_params + ["Incorrect", inc])
                                writer.writerow([N, QBER, delta_p, dpn] + cascade_params + ["Correct Power Only", power_only_cor])
                                writer.writerow([N, QBER, delta_p, dpn] + cascade_params + ["Incorrect Power Only", power_only_inc])


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
                                cor, inc, power_only_cor, power_only_inc = power_decoding.run_EC_power_flip(N,QBER,sigma,2,cascade_params, delta_p = delta_p, initial_only=False)
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
                            cor, inc, power_only_cor, power_only_inc = power_decoding.run_EC_power_flip(N,QBER,sigma,2,cascade_params,initial_only=False)
                            correct.append(cor/N)
                            incorrect.append(inc/N)
                            correct_power_only.append(power_only_cor/N)
                            incorrect_power_only.append(power_only_inc/N)
                        writer.writerow([N] + cascade_params + [sigma,"Correct"] + correct)
                        writer.writerow([N] + cascade_params + [sigma,"Incorrect"] + incorrect)
                        writer.writerow([N] + cascade_params + [sigma,"Correct Power Only"] + correct_power_only)
                        writer.writerow([N] + cascade_params + [sigma,"Incorrect Power Only"] + incorrect_power_only)


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


  

###plots

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


def main(exp_file):
    with open(exp_file) as json_file:
        exp_params = json.load(json_file)
    
    if(exp_params["experiement"] == "exp2"):
        exp_general(exp_params)

    print("Hello")


if __name__ == "__main__":
    main("exp_setup.json")
    #if sys.argv[1].endswith(".json"):
     #   main(sys.argv[1])
    #else:
     #   print("Invlaid Experiemt JSON")

#initial_discovery_exp()
#exp_vary_QBER()
#plot_diff_alg()
#plot_exp_vary_QBER()
#plot_exp_vary_sigma()
#exp_reconcilliation_efficiency()
#exp_vary_del_p()
#plot_exp_vary_del_p()
#exp_delp_to_sig_ratio()
#plot_del_p_to_sig_ratio()
#N,QBER,sigma,cascade_params = 250,0.15,0.8,["original","bernoulli"]
#run_EC_power_flip(N,QBER,sigma,2,cascade_params,initial_only=False)