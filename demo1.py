import opt_prob

name = '2.4 GOLDPR'
problem = opt_prob.Cons(name)

print('-- Problem --')
print(problem)
print('-- Solution --')
print(problem.xopt)
print(problem.fopt)
