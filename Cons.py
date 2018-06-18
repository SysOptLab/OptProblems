"""
Constrained Optimization Problem Collections:
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

import math

class Setup:

    """Setup problem with name.

    Args:
        name (str): Problem's name.

    Attributes:
        obj (function): Objective function.
        cns (function): Constraint function.
        lb (list): Lower bound of variables.
        ub (list): Upper bound of variables.
        xopt (list): Solution's variables.
        fopt (float): Solution's objective value.
    """

    def __init__(self, name):

        if name == 'ALKYLATION':
            
            self.__doc__ = """ALKYLATION"""

            def obj(x):
                X1 = x[0]
                X2 = x[1]
                X3 = x[2]
                X4 = x[3]
                X5 = x[4]
                x5 = 1.22*X4 - X1
                f = -(0.063*X4*X5 - 5.04*X1 - 0.035*X2 - 10.0*X3 - 3.36*x5)
                return f

            def cns(x):
                X1 = x[0]
                X2 = x[1]
                X3 = x[2]
                X4 = x[3]
                X5 = x[4]
                X6 = x[5]
                X7 = x[6]
                x5 = 1.22*X4 - X1
                x6 = (98000*X3)/(X4*X6 + 1000.0*X3)
                x8 = (X2 + x5)/X1

                g1 = 0.99*X4 - (X1*(1.12 + 0.13167*x8 - 0.00667*x8**2))
                g2 = (X1*(1.12 + 0.13167*x8 - 0.00667*x8**2)) - (100.0/99.0)*X4
                g3 = 0.99*X5 - (86.35 + 1.098*x8 - 0.038*x8**2 + 0.325*(x6 - 89.0))
                g4 = (86.35 + 1.098*x8 - 0.038*x8**2 + 0.325*(x6-89.0)) - (100.0/99.0)*X5
                g5 = 0.9*X6 - (35.82 - 0.222*X7)
                g6 = (35.82 - 0.222*X7) - (10.0/9.0)*X6
                g7 = 0.99*X7-(-133+3*X5)
                g8 = (-133 + 3.0*X5) - (100.0/99.0)*X7
                g9 = x5 - 2000
                g10 = -x5
                g11 = x6 - 93.0
                g12 = 85.0 - x6
                g13 = x8 - 12.0
                g14 = 3.0 - x8
                return [g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, g11, g12, g13, g14]

            self.obj = obj
            self.cns = cns
            self.lb = [0, 0, 0, 0, 90, 0.01, 145]
            self.ub = [2000, 16000, 120, 5000, 95, 4, 162]
            self.xopt = [1698.1, 15819, 54.107, 3031.2, 95.000, 1.5618, 153.54]
            self.fopt = 1768.75

        elif name == 'CAMEL':

            self.__doc__ = """
            Three Hump Camel-back w/ added constraint
            J.W.Hardy. An implemented extension of Branin's method.  In L.C.W. Dixon and
            G.P. Szego (Eds.), Towards Global Optimization. pp 117-142.
            North-Holland, Amsterdam, 1975.
            """
            
            def obj(x):
                f = 2.0*x[0]**2 - 1.05*x[0]**4 + (x[0]**6)/6.0 - x[0]*x[1] + x[1]**2
                return f

            def cns(x):
                g = (3.0-x[0])**2 + (1-x[1])**2 - 3.0
                return g

            self.obj = obj
            self.cns = cns
            self.lb = [-3, -1.5]
            self.ub = [ 3,  1.5]
            self.xopt = [1.7476, 0.8738]
            self.fopt = 0.2986

        elif name == 'FUNC2D':

            self.__doc__ = """FUNC2D"""

            def obj(x):
                f = 2 + 0.01*(x[1]-x[0]**2)**2 + (1-x[0])**2 + 2*(2-x[1])**2 + 7*math.sin(0.5*x[0])*math.sin(0.7*x[1]*x[0])
                return f

            def cns(x):
                g = -math.sin(x[0] - x[1] - math.pi/8.0)
                return g

            self.obj = obj
            self.cns = cns
            self.lb = [0, 0]
            self.ub = [5, 5]
            self.xopt = [2.7450, 2.3523]
            self.fopt = -1.1743

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
                f = (1.0+((x[0]+x[1]+1.0)**2)*g1)*(30.0+((2.0*x[0]-3.0*x[1])**2)*g2)
                f = math.log(f)
                return f

            def cns(x):
                g1 = -(3.0*x[0]) + (-3.0*x[1])**3
                g2 = x[0] - x[1] - 1.0
                return [g1, g2]

            self.obj = obj
            self.cns = cns
            self.lb = [-2, -2]
            self.ub = [2, 2]
            self.xopt = [0.5955, -0.4045]
            self.fopt = 5.6694

        elif name == 'GOMEZ':
            
            self.__doc__ = """GOMEZ"""

            def obj(x):
                f = (4-2.1*(x[0]**2)+(x[0]**4)/3)*(x[0]**2) + x[0]*x[1] + (-4+4*(x[1]**2))*(x[1]**2)
                return f

            def cns(x):
                g = -math.sin(4*math.pi*x[0]) + 2*(math.sin(2*math.pi*x[1])**2)
                return g

            self.obj = obj
            self.cns = cns
            self.lb = [-1, -1]
            self.ub = [ 1,  1]
            self.xopt = [0.10925714458181, -0.62344776471809]
            self.fopt = -0.9711

        elif name == 'HS23':
            
            self.__doc__ = """HS23"""

            def obj(x):
                f = x[0]**2 + x[1]**2
                return f

            def cns(x):
                g1 = -x[0] - x[1] + 1.0
                g2 = -x[0]**2 - x[1]**2 + 1.0
                g3 = -9.0*x[0]**2 - x[1]**2 + 9.0
                g4 = -x[0]**2 + x[1]
                g5 = -x[1]**2 + x[0]
                return [g1, g2, g3, g4, g5]

            self.obj = obj
            self.cns = cns
            self.lb = [-50, -50]
            self.ub = [50, 50]
            self.xopt = [1.0, 1.0]
            self.fopt = 2.0

        elif name == 'HS66':
            
            self.__doc__ = """HS66"""

            def obj(x):
                f = 0.2*x[2] - 0.8*x[1]
                return f

            def cns(x):
                g1 = math.exp(x[0]) - x[1]
                g2 = math.exp(x[1]) - x[2]
                return [g1, g2]
            
            self.obj = obj
            self.cns = cns
            self.lb = [0, 0, 0]
            self.ub = [100, 100, 10]
            self.xopt = [0.1841, 1.2022, 3.3273]
            self.fopt = 0.5182

        elif name == 'KS224':
            
            self.__doc__ = """
            Klaus Schittkowski Problem Collection 224
            """

            def obj(x):
                f = 2.0*x[0]**2 + x[1]**2 - 48.0*x[0] - 40.0*x[1]
                return f

            def cns(x):
                g1 = -1.0*(x[0] + 3.0*x[1])
                g2 = -1.0*(18.0 - x[0] - 3*x[1])
                g3 = -1.0*(x[0] + x[1])
                g4 = -1.0*(8.0 - x[0] - x[1])
                return [g1, g2, g3, g4]

            self.obj = obj
            self.cns = cns
            self.lb = [0, 0]
            self.ub = [6, 6]
            self.xopt = [4.0, 4.0]
            self.fopt = 304.0

        elif name == 'KS250':
            
            self.__doc__ = """
            Klaus Schittkowski Problem Collection p.74
            """

            def obj(x):
                f = -x[0]*x[1]*x[2]
                return f

            def cns(x):
                g1 = - x[0] - 2*x[1] - 2*x[2]
                g2 = x[0] + 2*x[1] + 2*x[2] - 72.0
                return [g1, g2]

            self.obj = obj
            self.cns = cns
            self.lb = [0, 0, 0]
            self.ub = [20, 11, 42]
            self.xopt = [20.0, 11.0, 15.0]
            self.fopt = -3300.0

        elif name == 'KS346':

            self.__doc__ = """
            Klaus Schittkowski Problem Collection p.167
            """

            def obj(x):
                f = -(0.0201/1e7)*(x[0]**4)*x[1]*(x[2]**2)
                return f

            def cns(x):
                g1 = (x[0]**2)*x[1] - 675.0
                g2 = (x[0]**2)*(x[2]**2)/1e7 - 0.419
                return [g1, g2]

            self.obj = obj
            self.cns = cns
            self.lb = [0, 0, 0]
            self.ub = [36, 5, 125]
            self.xopt = [16.51, 2.477, 124]
            self.fopt = -5.68478

        elif name == 'NEWBRANIN':

            self.__doc__ = """
            Branin test function
            F.H.Branin. Widely convergent method for finding multiple solutions of simultaneous
            nonlinear equations.  IBM Journal of Research and Development, 16:504-522, 1972.
            """

            def obj(x):
                f = -(x[0]-10.0)**2 - (x[1]-15.0)**2
                return f

            def cns(x):
                a = 1.0
                b = 5.1/(4.0*(math.pi**2))
                c = 5.0/math.pi
                d = 6.0
                e = 10.0
                f = 1.0/(8.0*math.pi)
                branin = a*(x[1] - b*x[0]**2 + c*x[0] - d)**2 + e*(1-f)*math.cos(x[0]) + e
                g = branin - 5.0
                return g

            self.obj = obj
            self.cns = cns
            self.lb = [-5, 0]
            self.ub = [10, 15]
            self.xopt = [3.2730, 0.0489]
            self.fopt = -268.789

        elif name == 'PRES':

            self.__doc__ = """PRES"""

            def obj(x):
                f = x[0] + x[1]
                return f

            def cns(x):
                g1 = 20.0 - (x[0]**2)*x[1]
                g2 = x[0]**2 + 8.0*x[1] - 75.0
                g3 = 1.0 - ((x[0] + x[1] - 5.0)**2)/30.0 - ((x[0] - x[1] - 12.0)**2)/120.0
                return [g1, g2, g3]

            self.obj = obj
            self.cns = cns
            self.lb = [0, 0]
            self.ub = [10, 10]
            self.xopt = [3.1139, 2.0627]
            self.fopt = 5.1765
        
        elif name == 'SPEEDREDUCER':

            self.__doc__ = """
            Floudas & Pardalos Ex 11.3 p. 153
            Weight Min of Speed Reducer
            originally from
            Chew, C.S., and Q. Zheng. Integral Global Optimization, vol 298 of Lecture
            Notes in Economics and Mathematical Systems. Springer-Verlag, 1988.
            """

            def obj(x):
                f = 0.7854*x[0]*(x[1]**2)*(3.3333*(x[2]**2)+14.9334*x[2]-43.0934) - \
                    1.508*x[0]*((x[5]**2)+(x[6]**2)) + \
                    7.477*((x[5]**3)+(x[6]**3)) + \
                    0.7854*(x[3]*(x[5]**2)+x[4]*(x[6]**2))
                return f

            def cns(x):
                A1 = math.sqrt((745*x[3]/x[1]/x[2])**2 + (16.911)**6)
                A2 = math.sqrt((745*x[4]/x[1]/x[2])**2 + (157.51)**6)
                B1 = 0.1*x[5]**3
                B2 = 0.1*x[6]**3

                g1 = 27.0 - x[0]*(x[1]**2)*x[2]
                g2 = 397.5 - x[0]*(x[1]**2)*(x[2]**2)
                g3 = 1.93 - x[1]*(x[5]**4)*x[2]/(x[3]**3)
                g4 = 1.93 - x[1]*(x[6]**4)*x[2]/(x[4]**3)
                g5 = A1*B1 - 1101.0
                g6 = A2/B2 - 850.0
                g7 = x[1]*x[2] - 40.0
                g8 = 5.0 - x[0]/x[1]
                g9 = x[0]/x[1] - 12.0
                g10 = 1.5*x[5] - x[3] + 1.9
                g11 = 1.5*x[6] - x[4] + 1.9
                return [g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, g11]

            self.obj = obj
            self.cns = cns
            self.lb = [2.6, 0.7, 17, 7.3, 7.3, 2.9, 5.0]
            self.ub = [3.6, 0.8, 28, 8.3, 8.3, 3.9, 5.5]
            self.xopt = [3.5, 0.7, 17, 7.3, 7.71, 3.35, 5.287]
            self.fopt = 2994.47

        else:
            raise "Unkown problem name."

    def __str__(self):
        return self.__doc__


if __name__ == '__main__':

    problem = Setup('GOLDPR')
    x0 = [0.0, 0.0]

    print('obj(x0) = {}'.format(problem.obj(x0)))
    print('cns(x0) = {}'.format(problem.cns(x0)))
    print('lb = {}'.format(problem.lb))
    print('ub = {}'.format(problem.ub))
    print('xopt = {}'.format(problem.xopt))
    print('fopt = {}'.format(problem.fopt))
