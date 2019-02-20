# OptProblem

Python package for optimization problem collections, created by NTU ME SOLab.

# How to Use

## Solution

```python
import OptProblem

name = '2.4 GOLDPR'
problem = OptProblemSet.Cons(name)

print('-- solution --')
print(problem.xopt)
print(problem.fopt)
```

## scipy.optimize

```python
import numpy as np
import OptProblem
import scipy.optimize

# -- problem setup

name = '2.4 GOLDPR'
problem = OptProblemSet.Cons(name)

def cns(x):
    g = -1.0*np.array(problem.cns(x))
    return g.tolist()

# -- start optimization

x0 = ((np.array(problem.lb) + np.array(problem.ub)) / 2.0).tolist()

bounds = []
for lb_i, ub_i in zip(problem.lb, problem.ub):
    bounds.append((lb_i, ub_i))

ineq_cons = {'type':'ineq', 'fun': cns}
method = 'SLSQP'
options = {'disp': True}

res = scipy.optimize.minimize(problem.obj, x0, method=method, bounds=bounds,
                              constraints=ineq_cons, options=options)
```

## SOLab DIRECT Algorithm

```python
import DIRECT
import OptProblem
import warnings; warnings.filterwarnings("ignore")

# -- problem setup

name = '2.5 GOMEZ'
problem = OptProblemSet.Cons(name)

# -- start optimization

options = {}
solver = DIRECT.Solver(problem, options=options)
solver.optimize()
```

# Reference

1. https://www.sfu.ca/~ssurjano/optimization.html
2. http://infinity77.net/global_optimization/test_functions.html#test-functions-index
3. http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO.htm