"""
Non-Constrained Optimization Problem Collections:
6-Hump-Camel
GOLDPR
"""

import math

class Setup:

    """Setup problem with name.

    Args:
        name (str): Problem's name.

    Attributes:
        obj (function): Objective function.
        lb (list): Lower bound of variables.
        ub (list): Upper bound of variables.
        xopt (list): Solution's variables.
        fopt (float): Solution's objective value.
    """

    def __init__(self, name):

        if name == '6-Hump-Camel':

            self.__doc__ = """6-Hump-Camel
            """
            
            def obj(x):
                f = (4.0 - 2.1*x[0]**2 + x[0]**(4.0/3.0))*x[0]**2 + x[0]*x[1] + (-4.0 + 4.0*x[1]**2)*x[1]**2
                return f

            self.obj = obj
            self.lb = [-3, -2]
            self.ub = [3, 2]
            self.xopt = []
            self.fopt = -1.0316

        elif name == 'GOLDPR':

            self.__doc__ = """
            Goldstein Price
            L.Pronzato, E.Walter, A.Venot, and J.F.Lebruchec. "A general purpose global
            optimizer: Implementation and applicaitons". Mathematics and Computers in Simulation,
            26:412-422, 1984.
            """

            def obj(x):
                g1 = 19.0 - 14.0*x[0] + 3.0*x[0]**2 - 14.0*x[1] + 6.0*x[0]*x[1] + 3.0*x[1]**2
                g2 = 18.0 - 32.0*x[0] + 12.0*x[0]**2 + 48.0*x[1] - 36.0*x[0]*x[1] + 27.0*x[1]**2
                f = (1.0 + ((x[0] + x[1] + 1.0)**2)*g1) * (30.0 + ((2.0*x[0] - 3.0*x[1])**2)*g2)
                return f

            self.obj = obj
            self.lb = [-2, -2]
            self.ub = [2, 2]
            self.xopt = []
            self.fopt = 3.0

        else:
            raise "Unkown problem name."

    def __str__(self):
        return self.__doc__


if __name__ == '__main__':

    problem = Setup('GOLDPR')
    x0 = [0.0, 0.0]

    print('obj(x0) = {}'.format(problem.obj(x0)))
    print('lb = {}'.format(problem.lb))
    print('ub = {}'.format(problem.ub))
    print('xopt = {}'.format(problem.xopt))
    print('fopt = {}'.format(problem.fopt))
