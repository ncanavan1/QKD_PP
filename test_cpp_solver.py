import cascade_EC
import noisy_decoding
import LinAlg
import matplotlib.pyplot as plt
import numpy as np

N = 10
X,Y = cascade_EC.gen_sifted_keys(N,0.1)
Y_ec, Eves_traces = cascade_EC.cascade_EC(X,Y,0.1,3,0.1,1)

Included_matrix, P_A, H_E, Conf_vector, partial_key = noisy_decoding.gen_system_from_trace(Eves_traces,N)


linsys = LinAlg.mod2system_Solver()
P_AT = np.atleast_2d(P_A).transpose()
npsys = np.append(Included_matrix,P_AT,axis = 1)

npsys = np.asarray(
        [[1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1],
        [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1],
        [0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1],
        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0],
        [1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1],
        [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1],
        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0],
        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0],
        [1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1],
        [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1]])

linsys.setSystem(npsys.shape[0],npsys.shape[1]-1,npsys)
linsys.gauss_jordan_elimination_positional()

solutions = np.zeros([2,2])
solutions = linsys.find_solutions()
print(solutions)