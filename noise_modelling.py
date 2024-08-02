import numpy as np
import matplotlib.pyplot as plt

def simulate_xor_trace(A, B, noise_std):
    # Calculate Hamming Weight of the XOR result
    xor_result = A ^ B
    hw = bin(xor_result).count('1')

    alpha = 0.4
    beta = -4.75
    
    # Base power consumption with Gaussian noise
    base_power = alpha*hw + beta
    noise = np.random.normal(0, noise_std)
    
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
    noise_range = np.arange(0,0.5,0.01)
    repeats = 1
    for i in range(4):
        A = runs[i][0]
        B = runs[i][1]
        performance = []
        for sigma in noise_range:
            pwr_avg = 0
            for rep in range(repeats):
                pwr = simulate_xor_trace(A,B,sigma)
                pwr_avg = pwr_avg + pwr
            pwr_avg = pwr_avg/repeats
            performance.append(pwr_avg)
        plt.plot(noise_range,performance,label="A={0},B={1}".format(A,B))

    plt.plot(noise_range,np.ones(noise_range.shape[0])*(-4.55),"k--")
    plt.legend()
    plt.ylabel("Power (units?)")
    plt.xlabel("Noise")
    plt.show()

#run_noise_model()
    
