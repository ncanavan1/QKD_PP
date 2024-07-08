import numpy as np
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
#plt.style.use('seaborn-poster')

res = 1
x = np.arange(0,8.1,res) ##number of correct
y = np.arange(0,8.1,res) ##number of incorrect
X, Y = np.meshgrid(x, y)
sigma = 0.1




P = ((1-sigma)**(X/(X+Y)))*(sigma**(Y/(X+Y)))
R_rat = (Y-X)/(Y+x)

Rc = Y/(Y+X) 
Ri = X/(Y+X) 



#conf = (Rc**(1-sigma) * Ri**(sigma)) - (Ri**(1-sigma) * Rc**(sigma))
#conf = (X**(1-sigma) * Y**(sigma)) - (Y**(1-sigma) * X**(sigma))
#conf[0,:] = np.log10(X[0,:] + 1)
#conf[:,0] = np.log10(Y[:,0] + 1)
#alpha = 0.4
#conf = ((1-sigma)**(X) * sigma**(Y))/((1-sigma)**(Y) * sigma**(X))
#conf = np.emath.logn(alpha,conf)

#diff = np.abs(X-Y)

#conf = (np.tanh(alpha*diff))**(1/(1-sigma)**2)

#for i in range(P.shape[0]):
#    for j in range(P.shape[1]):
#        if P[i,j] > 1:
#            P[i,j] = 1 


fig = plt.figure(figsize = (12,10))
ax = plt.axes(projection='3d')
#surf = ax.plot_surface(X, Y, P, cmap = plt.cm.cividis)
surf2 = ax.plot_surface(X,Y,P, cmap = plt.cm.magma)
ax.set_xlabel('#Correct', labelpad=20)
ax.set_ylabel('#Incorrect', labelpad=20)
ax.set_zlabel('Confidence', labelpad=20)
plt.show()