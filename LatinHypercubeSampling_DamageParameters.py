"""
@author: jbrondex

Description:
------------
This file uses the hypercube latin sampling method to get a set of combinations of parameters (Sigmath_ice, B, lambdah) for the damage model
that correspond to an optimal coverage of the parameter space

"""
import numpy as np
import matplotlib.pyplot as plt
from smt.sampling_methods import LHS

###Define parameter space
Sigmath_limits = [0.05, 0.2] ###Values deduced from Teterousse steady state with empty cavity (2010 geometry) compared to crevasses of August 2011
B_limits = [0.5, 2.0] ###Values of Krug 2014
###For lambdah we consider 3 discrete values
Lambdah_values = [0.0, 0.5, 1.0] ###0.4 in Pralong 2006
xlimits = np.array([Sigmath_limits, B_limits])
sampling = LHS(xlimits=xlimits)

### How many combinations (Sigmath, B) do we want to consider for each value of lambdah ?
num = 8
x = sampling(num)
print(x.shape)

### Create vectors X, Y, Z: X,Y are sampled values of Sigmath,B while Z are discrete values of lambdah
X=[]
Y=[]
Z=[]
for i,lambdah in enumerate(Lambdah_values):
    X.append(x[:, 0])
    Y.append(x[:, 1])
    z=lambdah*np.ones((1,num))
    Z.append(z)

### Show combinations in a Fig.
fig = plt.figure()
ax = fig.add_subplot(projection='3d')

ax.scatter(X, Y, Z,"o")
ax.set_xlabel(r'$\sigma_{th}~(\mathrm{MPa})$')
ax.set_ylabel(r'$B~\mathrm{(MPa^{-1}~a^{-1})}$')
ax.set_zlabel(r'$\lambda_h~\mathrm{(_)}$')
plt.show()