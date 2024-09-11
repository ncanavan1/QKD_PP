import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

def simulate_xor_trace_stm32_kim(A, B, noise_std):
    # Calculate Hamming Weight of the XOR result
    xor_result = A ^ B
    hw = bin(xor_result).count('1')

    alpha = 0.4
    beta = -4.75
    
    # Base power consumption with Gaussian noise
    base_power = alpha*hw + beta
    noise = np.abs(np.random.normal(0, noise_std))
    
    # Simulated power trace value
    simulated_power = base_power + noise
    
    return simulated_power

def simulate_xor_trace_general(A, B, noise_std):
    # Calculate Hamming Weight of the XOR result
    xor_result = A ^ B
    hw = bin(xor_result).count('1')

    alpha = 1
    beta = -0
    
    # Base power consumption with Gaussian noise
    base_power = alpha*hw + beta
    noise = np.abs(np.random.normal(0, noise_std))
    
    # Simulated power trace value
    simulated_power = base_power + noise
    
    return simulated_power



"""
# Example usage
A = 1  # Example binary input A
B = 0  # Example binary input B
noise_std = 0.3  # Standard deviation of noise

trace = simulate_xor_trace(A, B, noise_std)
print(trace)
"""

def run_noise_model():
    
    run1 = [0,0]
    run2 = [0,1]
    run3 = [1,0]
    run4 = [1,1]
    runs = [run1,run2,run3,run4]
    noise_range = np.arange(0,1,0.01)
    repeats = 1
    for i in range(4):
        A = runs[i][0]
        B = runs[i][1]
        performance = []
        for sigma in noise_range:
            pwr_avg = 0
            for rep in range(repeats):
                #pwr = simulate_xor_trace_stm32_kim(A,B,sigma)
                pwr = simulate_xor_trace_general(A,B,sigma)
                pwr_avg = pwr_avg + pwr
            pwr_avg = pwr_avg/repeats
            performance.append(pwr_avg)
        plt.plot(noise_range,performance,label="A={0},B={1}".format(A,B))


    average_pwr = 0.5
    noisy_average = np.ones(noise_range.shape[0])*average_pwr
    for i in range(noise_range.shape[0]):
        noisy_average[i] = noisy_average[i] + 0.8*noise_range[i]

    plt.plot(noise_range,noisy_average,"k--")
    #plt.title("Power consumption of an XOR operation on \nSTM32-F405RGT6 (ARM Coretex-M4)")
    plt.title("Power consumption of a generalised XOR operation")

    plt.legend()
    plt.ylabel("Power, mW")
    plt.xlabel("Noise, $\sigma$")
    #plt.savefig("results/XOR_power_consumption_kim.png",bbox_inches = "tight")
    plt.savefig("results/XOR_power_consumption_general.png",bbox_inches = "tight")




def calc_cdfs(x,a_base,b_base,sigma):
    cdf_a = norm.cdf(x,a_base,sigma)
    cdf_b = norm.cdf(x,b_base,sigma)
    return cdf_a, cdf_b



def prob_y_x_unknown(y,sigma,cdf_a,cdf_b,power_resolution): #p(power=y)
    ##setting power resolution to integral limits

    y_below = round(y,power_resolution)
    y_below_pos = np.where(x == y_below)[0][0]
    y_above = y_below+10**(-power_resolution)
    y_above_pos = y_below_pos + 1

    p_a = cdf_a[y_above_pos] - cdf_a[y_below_pos]
    p_b = cdf_b[y_above_pos] - cdf_b[y_below_pos]
    p_unkown = 0.5*(p_a+p_b)
    print(p_a) 
    print(p_b) 
    print(p_unkown) 

    pr_x_given_y = (p_b * 0.5)/p_unkown
    print(pr_x_given_y)    
    
    plt.xlabel("Power")
    plt.ylabel("CDF")
    plt.plot(x,cdf_a,label="x=a")
    plt.plot(x,cdf_b,label="x=b")
    plt.legend()
    plt.grid(True)
  #  plt.show()


##set y to at least one magnitude higher resolution than power resolution
start_power = -3
end_power = 3
power_resolution = 2
a_base = -0.5
b_base = 0.5

sigma = 0.5
y = 0.7

x = np.arange(start_power,end_power,10**(-power_resolution))
x = np.round(x,decimals=power_resolution)
cdf_a, cdf_b = calc_cdfs(x,a_base,b_base,sigma)



prob_y_x_unknown(y,sigma,cdf_a,cdf_b,power_resolution)
#run_noise_model()
    
