import numpy as np
import OptProblemSet
import scipy.optimize

# -- problem

name = '2.4 GOLDPR'
problem = OptProblemSet.Cons(name)
obj = problem.obj

def cns(x):
    g = -1.0*np.array(problem.cns(x))
    return g.tolist()

# --

x0 = ((np.array(problem.lb) + np.array(problem.ub)) / 2.0).tolist()

bounds = []
for lb_i, ub_i in zip(problem.lb, problem.ub):
    bounds.append((lb_i, ub_i))

ineq_cons = {'type':'ineq', 'fun': cns}
method = 'SLSQP'
options = {'disp': True}

res = scipy.optimize.minimize(obj, x0, method=method, bounds=bounds, constraints=ineq_cons, options=options)

print('-- optimal result --')
print('res.success = {}'.format(res.success))
print(res.x)
print(res.fun)

print('-- solution --')
print(problem.xopt)
print(problem.fopt)