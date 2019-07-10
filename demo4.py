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

print(res)