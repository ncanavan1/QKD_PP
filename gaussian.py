import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from mpl_toolkits import mplot3d

def get_gx_guassian(x,sigma,mu):
    #prob = (1/(np.sqrt(2*np.pi*sigma**2)))*np.exp((-1*(x)**2)/(2*sigma**2))
    #prob = np.random.normal(0,sigma)
    #prob = norm(0,sigma).cdf(x)
    prob = (1/(sigma*np.sqrt(2*np.pi))) * np.exp(-1* ((x-mu)/(2*sigma))**2)
    return prob


def plot_gaussian():
    plt.style.use("seaborn-v0_8-darkgrid")

    mu1 = 0#-0.5
    mu2 = 1#0.5
    sigma_r = np.arange(0.05,1,0.001)
    x = np.arange(-2,3,0.01)
    fx_r = []
    fx_1 = []
    fx_2 = []
    for sigma in sigma_r:
        fx1 = get_gx_guassian(x,sigma,mu1)
        fx2 = get_gx_guassian(x,sigma,mu2)
        fx_1.append(fx1)
        fx_2.append(fx2)
        fx = 0.5*(fx1+fx2)
        fx_r.append(fx)
    fx_r = np.asarray(fx_r)
    fx_1 = np.asarray(fx_1)
    fx_2 = np.asarray(fx_2)
    X,Y = np.meshgrid(x,sigma_r)
    fig = plt.figure()
    ax = plt.axes(projection="3d")
    #ax.plot_surface(X,Y,fx_1,label="x=a")
    #ax.plot_surface(X,Y,fx_2,label="x=b")
    ax.plot_surface(X,Y,fx_r,label="x={a,b}",color="g")
    ax.set_xlabel("Power")
    ax.set_ylabel("Noise, $\sigma$")
    ax.set_zlabel("Probability Density")
    ax.legend()

    

    #plt.show()


    plt.figure(2,figsize=(7,4))
    plt.subplot(1,2,1)
    plt.plot(x,fx_1[0],label=(r'XOR$\rightarrow$a'))
    plt.plot(x,fx_2[0],label=(r'XOR$\rightarrow$b'))
    plt.plot(x,fx_r[0],label=(r'XOR$\rightarrow${a,b}'))
    plt.xlabel("Power")
    plt.ylabel("Probability Density")
    plt.axvline(x=0.5,color="k",linestyle="--")
  #  plt.axvline(x=1,color="k",linestyle="--")
   # plt.legend(loc="upper right")
    plt.title("$\sigma$ = 0.01")

    plt.subplot(1,2,2)
    plt.plot(x,fx_1[-1],label=(r'XOR$\rightarrow$a'))
    plt.plot(x,fx_2[-1],label=(r'XOR$\rightarrow$b'))
    plt.plot(x,fx_r[-1],label=(r'XOR$\rightarrow${a,b}'))
    plt.xlabel("Power")
   # plt.ylabel("Probability Density")
    plt.axvline(x=0.5,color="k",linestyle="--")
   # plt.axvline(x=1,color="k",linestyle="--")

    plt.title("$\sigma$ = 1")
    plt.legend(bbox_to_anchor=(-0.1,1.25),loc="upper center")
    plt.savefig("results/proability_convergance.png",bbox_inches="tight")

plot_gaussian()