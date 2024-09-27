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
    beta = 0
    
    # Base power consumption with Gaussian noise
    base_power = alpha*hw + beta
    #noise = np.abs(np.random.normal(0, noise_std))
    noise = np.random.normal(0,noise_std)
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


##retrurns:- base prob of power =y
## - prob of power = y given x=a
## - prob of power = y given x=b
def probs_given_y(x,y,cdf_a,cdf_b,power_resolution): #p(power=y)
    ##setting power resolution to integral limits

    y_below = round(y,power_resolution)
    y_below_pos = np.where(x == y_below)[0][0]
    y_above = y_below+10**(-power_resolution)
    y_above_pos = y_below_pos + 1

    p_a = cdf_a[y_above_pos] - cdf_a[y_below_pos]
    p_b = cdf_b[y_above_pos] - cdf_b[y_below_pos]
    cdf_u = 0.5*(cdf_a + cdf_b)
    #p_unkown = 0.5*(p_a+p_b)
    p_u = cdf_u[y_above_pos] - cdf_u[y_below_pos]
    #pr_y_given_a = (p_a * 0.5)/p_unkown
    """
    plt.style.use("seaborn-v0_8-darkgrid")

    plt.xlabel("Power")
    plt.ylabel("CDF")
    plt.plot(x,cdf_a,label=(r'XOR$\rightarrow$ a'))
    plt.plot(x,cdf_b,label=(r'XOR$\rightarrow$ b'))


    x_cord = [0,0.3]
    y_cord_a = []
    y_cord_b = []
    y_pos = np.where(x == x_cord[0])[0][0]
    y_cord_a.append(cdf_a[y_pos])
    y_cord_b.append(cdf_b[y_pos])
    y_pos = np.where(x == x_cord[1])[0][0]
    y_cord_a.append(cdf_a[y_pos])
    y_cord_b.append(cdf_b[y_pos])


    plt.axvline(x_cord[0],color="k",linestyle="--",ymax=y_cord_a[0])
    plt.axvline(x_cord[1],color="k",linestyle="--",ymax=y_cord_a[1])
    plt.axhline(y_cord_a[0],color="k",linestyle="--",xmax=3/7)
    plt.axhline(y_cord_a[1],color="k",linestyle="--",xmax=3.3/7)
    plt.axhline(y_cord_b[0],color="k",linestyle="--",xmax=3/7)
    plt.axhline(y_cord_b[1],color="k",linestyle="--",xmax=3.3/7)

    plt.text(0.19,-0.1,r'$\Delta_{pwr}$')
    plt.legend()
    plt.savefig("results/cdf.png")
    plt.show()
    """

    return p_u, p_a, p_b

def plot_prob_x_range(x,cdf_a,cdf_b,power_resolution):
    p_as = []
    p_bs = []
    p_unkowns = []

    for y in x[:-1]:
        pu, pa, pb = probs_given_y(x,y,cdf_a,cdf_b,power_resolution)
        p_as.append((pa*0.5)/pu)
        p_bs.append((pb*0.5)/pu)
        p_unkowns.append(pu)
    plt.plot(x[:-1],p_as)
    plt.plot(x[:-1],p_bs)

    plt.show() 


def prob_given_S(S,x,cdf_a,cdf_b,power_resolution):

    Scount = len(S)
    p_us = []
    p_as = []
    p_bs = []

    for i in range(Scount):
        p_u, p_a, p_b = probs_given_y(x,S[i],cdf_a,cdf_b,power_resolution)
        p_us.append(p_u)
        p_as.append(p_a)
        p_bs.append(p_b)

    p_y_given_a = 1
    p_y_given_b = 1
   # for i in range(Scount):
    #    p_y_given_a = p_y_given_a*(p_as[i])/(p_as[i] + p_bs[i])
     #   p_y_given_b = p_y_given_b*(p_bs[i])/(p_as[i] + p_bs[i])


    p_s_given_a = np.prod(np.asarray(p_as))
    p_s_given_b = np.prod(np.asarray(p_bs))
   # p_s_unknown = np.prod(np.asarray(p_us))

    p_y_given_a = p_s_given_a/(p_s_given_a+p_s_given_b)
    p_y_given_b = p_s_given_b/(p_s_given_a+p_s_given_b)


   # for i in range(Scount):
    #    p_y_given_a = p_y_given_a*(p_as[i]/p_us[i])*(0.5)
     #   p_y_given_b = p_y_given_b*(p_bs[i]/p_us[i])*(0.5)

    #plt.show()
   # print(p_y_given_a)
   # print(p_y_given_b)
   # print(p_y_given_a+p_y_given_b)
   # plt.hist(S)
   # plt.show()
    return p_y_given_a, p_y_given_b


def runner():

    ##set y to at least one magnitude higher resolution than power resolution
    start_power = -3
    end_power = 4
    power_resolution = 3
    a_base = 0
    b_base = 1

    sigma = 0.5

    x = np.arange(start_power,end_power,10**(-power_resolution))
    x = np.round(x,decimals=power_resolution)
    cdf_a, cdf_b = calc_cdfs(x,a_base,b_base,sigma)


    S = []
    A,B = 0, 0
    for i in range(4):
        sim_power = simulate_xor_trace_general(A,B,sigma)
        S.append(sim_power)

    prob_given_S(S,x,cdf_a,cdf_b,power_resolution)
    #p_any, p_a, p_b = probs_given_y(x,y,cdf_a,cdf_b,power_resolution)
    #plot_prob_x_range(x,cdf_a,cdf_b,power_resolution)
    #run_noise_model()
#runner()