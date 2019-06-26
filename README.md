# opt-prob-collection

Optimization problem collections by NTU ME SOLab.

## Requirement

```
python 2.7 or 3
---
numpy
matplotlib
```

## How to Use

### Problem list

```python
import opt_prob

print(opt_prob.Cons.names)
print(opt_prob.NonCons.names)
```

```
>>>
-- Constrained Problem --
['1.4 G4 Problem', '1.6 G6 Problem', '1.7 G7 Problem', '1.8 G8 Problem', '1.9 G9 Problem', '1.10 G10 Problem', '2.1 ALKYLATION', '2.2 CAMEL', '2.3 FUNC2D', '2.4 GOLDPR', '2.5 GOMEZ', '2.6 HS23', '2.8 KS224', '2.9 KS250', '2.10 KS346', '2.11 NEWBRANIN', '2.12 PRES']
-- Non-Constrained Problem --
['1.1 Ackley Function', '1.2 Bukin Function N. 6', '1.3 Cross-in-Tray Function', '1.4 Drop-Wave Function', '1.5 Eggholder Function', '1.6 Gramacy and Lee (2012) Function', '1.7 Griewank Function', '1.8 Holder Table Function', '1.10 Levy Function', '1.11 Levy Function N. 13', '1.12 Rastrigin Function', '1.13 Schaffer Function N. 2', '1.14 Schaffer Function N. 4', '1.15 Schwefel Function', '1.16 Shubert Function', '2.1 Bohachevsky Function', '2.2 Perm Function', '2.3 Rotated Hyper-Ellipsoid Function', '2.4 Sphere Function Modified', '2.5 Sum of Different Powers Function', '2.6 Sum Squares Function', '2.7 Trid Function', '3.1 Booth Function', '3.2 Matyas Function', '3.3 McCormick Function', '3.5 Zakharov Function', '4.1 Three-Hump Camel Function', '4.2 Six-Hump Camel Function', '4.3 Dixon-Price Function', '4.4 Rosenbrock Function', '5.2 Easom Function', '5.3 Michalewicz Function', '6.1 Beale Function', '6.2 Branin Function', '6.3 Colville Function', '6.4 Forrester et al. (2008) Function', '6.5 Goldstein-Price Function', '6.6 Hartmann 3-D Function', '6.7 Hartmann 4-D Function', '6.8 Hartmann 6-D Function', '6.9 Perm Function', '6.11 Shekel Function 5', '6.12 Shekel Function 7', '6.13 Shekel Function 10', '6.14 Styblinski-Tang Function']
```

### Solution

```python
import opt_prob

name = '2.4 GOLDPR'
problem = opt_prob.Cons(name)

print('-- Problem --')
print(problem)
print('-- Solution --')
print(problem.xopt)
print(problem.fopt)
```

```
>>>
-- Problem --

2.4 GOLDPR

Dimensions: 2

Goldstein Price
L.Pronzato, E.Walter, A.Venot, and J.F.Lebruchec. "A general purpose global
optimizer: Implementation and applicaitons". Mathematics and Computers in Simulation,
26:412-422, 1984.


-- Solution --
[0.5955, -0.4045]
5.6694
```

### Ploting

```python
import opt_prob

name = '2.4 GOLDPR'
problem = opt_prob.Cons(name)

opt_prob.plot(problem)
```

![](./pic1.png)
![](./pic2.png)

### scipy.optimize

```python
import numpy as np
import opt_prob
import scipy.optimize

# -- problem setup

name = '2.4 GOLDPR'
problem = opt_prob.Cons(name)

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

### SOLab DIRECT Algorithm

```python
import direct
import opt_prob

# -- problem setup

name = '2.5 GOMEZ'
problem = opt_prob.Cons(name)

# -- start optimization

solver = direct.Solver(problem)
solver.optimize()
```

## Reference

1. https://www.sfu.ca/~ssurjano/optimization.html
2. http://infinity77.net/global_optimization/test_functions.html#test-functions-index
3. http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO.htm
