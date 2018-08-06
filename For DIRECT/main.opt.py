"""
Demo for using DIRECT optimization algorithm
"""

import DIRECT
import OptProblems
import warnings

warnings.filterwarnings("ignore")

# -- problem

name = '2.5 GOMEZ'
problem = OptProblems.Cons(name)
# name = "6.2 Branin Function"
# problem = OptProblems.NonCons(name)

# -- start optimization

options = {}
# options['fcount'] = 10000
# options['stop_criterion'] = 'fopt'
options['stop_criterion'] = ['fopt', 'fcount']

# options['debug_plot_fd'] = True
# options['debug_plot_fd_penalty'] = True
# options['debug_plot_aux'] = True
# options['debug_plot_2d'] = True
# options['debug_savefig'] = True
# options['debug_savetxt'] = True

options['surrogate_assisted'] = True
# options['surrogate_type'] = 'original'
# options['surrogate_type'] = 'linear'
# options['surrogate_strategy'] = 'grid'
# options['surrogate_strategy'] = 'once_dist'
# options['surrogate_check_vertex'] = False
# options['surrogate_fitting_sample'] = 10

# options['dist_method'] = 'longside'
# options['cns_method'] = 'penalty'

solver = DIRECT.Solver(problem, options=options)
solver.optimize()
