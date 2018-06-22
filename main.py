"""
Demo of OptProblems, plot 2d problems.

ALKYLATION
CAMEL
FUNC2D
GOLDPR
GOMEZ
HS23
HS66
KS224
KS250
KS346
NEWBRANIN
PRES
SPEEDREDUCER
"""

import matplotlib.pyplot as plt
import numpy as np
import OptProblems

# problem = OptProblems.Cons('CAMEL')
problem = OptProblems.NonCons('HARTMANN 3-D')
lb = problem.lb
ub = problem.ub
obj = problem.obj
cns = problem.cns

num = 40
x1 = np.linspace(lb[0], ub[0], num)
x2 = np.linspace(lb[1], ub[1], num)
x1_plot, x2_plot = np.meshgrid(x1, x2)
xs = np.hstack((x1_plot.reshape(-1,1), x2_plot.reshape(-1,1)))
fs = [obj(x) for x in xs]
f_plot = np.array(fs).reshape(num, num)

if cns is not None:
    gs = [cns(x) for x in xs]
    dim = len(gs[0]) if isinstance(gs[0], list) else 1
    g_plot = np.array(gs).reshape(num, num, dim)

plt.figure() # figsize=(12.80, 10.24)
plt.xlabel('$x_1$')
plt.ylabel('$x_2$')
plt.contour(x1_plot, x2_plot, f_plot, cmap=plt.cm.plasma) # plt.pcolormesh(x1_plot, x2_plot, f_plot)
plt.colorbar()

if cns is not None:
    for i in range(dim):
        plt.contour(x1_plot, x2_plot, g_plot[:,:,i], colors="k", levels=0, linestyles='dotted', lw=0.1)

plt.axis( [ lb[0], ub[0], lb[1], ub[1] ] )
plt.axis('equal')
plt.xlim(lb[0], ub[0])
plt.ylim(lb[1], ub[1])

plt.show()
