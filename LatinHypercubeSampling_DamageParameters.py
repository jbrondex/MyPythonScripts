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
Lambdah_limits = [0.0, 0.5] ###0.4 in Pralong 2006
xlimits = np.array([Sigmath_limits, B_limits, Lambdah_limits])
sampling = LHS(xlimits=xlimits)

### How many combinations do we want to consider ?
num = 20
x = sampling(num)
print(x.shape)

### Show combinations in a Fig.
fig = plt.figure()
ax = fig.add_subplot(projection='3d')

ax.scatter(x[:, 0], x[:, 1], x[:, 2],"o")
ax.set_xlabel(r'$\sigma_{th}~(\mathrm{MPa})$')
ax.set_ylabel(r'$B~\mathrm{(MPa^{-1}~a^{-1})}$')
ax.set_zlabel(r'$\lambda_h~\mathrm{(_)}$')
plt.show()