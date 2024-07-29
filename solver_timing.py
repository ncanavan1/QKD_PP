import cascade_EC as cascade
import linalgtools as my_solver
import time
import matplotlib.pylab as plt
import numpy as np

def build_matrix(N,traces):
    M = len(traces)
    H_E = np.zeros([M,N])
    P_A = np.zeros(M)
    Included_matrix = np.zeros([M,N])
    alpha = 0.4 ##emperic
    beta = 2.5


    check = 0
    for trace in traces:
        inter_parity = trace[0]
        pos = trace[1]
        out_parity = trace[2]
        P_A[check] = out_parity
        for i in reversed(range(0,pos.size)): ##this was 1,pos.size WHY??
            H_E[check,pos[i]] = inter_parity[i] ^ inter_parity[i-1]
            Included_matrix[check,pos[i]] = 1
        H_E[check,pos[0]] = inter_parity[0]
        check = check + 1
    
    return Included_matrix, P_A


def time_RREF(N,QBER,XOR_noise,max_iter):
    X,Y = cascade.gen_sifted_keys(N,QBER)
    Y_ec, Eves_traces = cascade.cascade_EC(X,Y,QBER,max_iter,XOR_noise,1)
    if (X == Y_ec).all():
        Included_matrix, P_A = build_matrix(N,Eves_traces)
        
        
        start = time.time()
        ####solve system in completely unknown
        row_ech, b, row_order = my_solver.row_echelon_form(Included_matrix.copy(),P_A.copy())
        partial_solution = my_solver.read_off_basic_variables(row_ech,b) 
        end = time.time()

        elapsed = end - start   
        print("Time for N={0}; {1}s".format(N,elapsed))
        return(elapsed)
    

Nlist = [128,256,512,1024,2048,4096]
tlist = []
for N in Nlist:
    tlist.append(time_RREF(N,0.1,0.05,3))
plt.plot(Nlist,tlist)
plt.xlabel("Key Size, N")
plt.ylabel("Execution Time, s")
plt.savefig("save_results/solvertiming.png")