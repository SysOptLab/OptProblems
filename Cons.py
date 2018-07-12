"""
Constrained Optimization Problem Collections:

1. http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO.htm

1.4 G4 Problem
1.6 G6 Problem
1.7 G7 Problem
1.8 G8 Problem
1.9 G9 Problem
1.10 G10 Problem

2. Others

2.1 ALKYLATION
2.2 CAMEL
2.3 FUNC2D
2.4 GOLDPR
2.5 GOMEZ
2.6 HS23
2.8 KS224
2.9 KS250
2.10 KS346
2.11 NEWBRANIN
2.12 PRES

"""

import math
import numpy as np

if __package__:
    from . import debug_plot
else:
    import debug_plot

class Cons:

    """Constrained optimization problem.

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

    def __init__(self, name, dimensions=2):

        if name == '1.4 G4 Problem':
            
            self.__doc__ = """
            1.4 G4 Problem

            Dimensions: 5
            """

            def obj(x):
                y = 5.3578547*x[2]**2+0.8356891*x[0]*x[4]+37.293239*x[0]-40792.141
                return y

            def cns(x):
                u = 85.334407+0.0056858*x[1]*x[4]+0.0006262*x[0]*x[3]-0.0022053*x[2]*x[4]
                g1 = -u
                g2 = u-92
                v = 80.51249+0.0071317*x[1]*x[4]+0.0029955*x[0]*x[1]+0.0021813*x[2]**2
                g3 = -v+90
                g4 = v-110
                w = 9.300961+0.0047026*x[2]*x[4]+0.0012547*x[0]*x[2]+0.0019085*x[2]*x[3]
                g5 = -w+20
                g6 = w-25
                return [g1, g2, g3, g4, g5, g6]
            
            self.obj = obj
            self.cns = cns
            self.lb = [78,33,27,27,27]
            self.ub = [102,45,45,45,45]
            self.xopt = [78,33,29.995,45,36.7758]
            self.fopt = -30665.539

        elif name == '1.6 G6 Problem':
            
            self.__doc__ = """
            1.6 G6 Problem

            Dimensions: 2
            """

            def obj(x):
                y = (x[0]-10.0)**3+(x[1]-20.0)**3
                return y

            def cns(x):
                g1 = -(x[0]-5)**2-(x[1]-5)**2+100
                g2 = (x[0]-6)**2+(x[1]-5)**2-82.81
                return [g1, g2]
            
            self.obj = obj
            self.cns = cns
            self.lb = [13, 0]
            self.ub = [100, 100]
            self.xopt = [14.095,0.84296]
            self.fopt = -6961.81388

        elif name == '1.7 G7 Problem':
            
            self.__doc__ = """
            1.7 G7 Problem

            Dimensions: 10
            """

            def obj(x):
                y = x[0]**2+x[1]**2+x[0]*x[1]-14*x[0]-16*x[1]+(x[2]-10)**2+ \
                    4*(x[3]-5)**2+(x[4]-3)**2+2*(x[5]-1)**2+5*x[6]**2+ \
                    7*(x[7]-11)**2+2*(x[8]-10)**2+(x[9]-7)**2+45
                return y

            def cns(x):
                g1 = 4*x[0]+5*x[1]-3*x[6]+9*x[7]-105
                g2 = 10*x[0]-8*x[1]-17*x[6]+2*x[7]
                g3 = -8*x[0]+2*x[1]+5*x[8]-2*x[9]-12
                g4 = 3*(x[0]-2)**2+4*(x[1]-3)**2+2*x[2]**2-7*x[3]-120
                g5 = 5*x[0]**2+8*x[1]+(x[2]-6)**2-2*x[3]-40
                g6 = 0.5*(x[0]-8)**2+2*(x[1]-4)**2+3*x[4]**2-x[5]-30
                g7 = x[0]**2+2*(x[1]-2)**2-2*x[0]*x[1]+14*x[4]-6*x[5]
                g8 = -3*x[0]+6*x[1]+12*(x[8]-8)**2-7*x[9]
                return [g1, g2, g3, g4, g5, g6, g7, g8]
            
            self.obj = obj
            self.cns = cns
            self.lb = (np.ones(10)*-10.0).tolist()
            self.ub = (np.ones(10)*10.0).tolist()
            self.xopt = [2.171996, 2.363683, 8.773926, 5.095984, 0.9906548, 1.430574,1.321644, 9.828726, 8.280092, 8.375927]
            self.fopt = 24.3062091

        elif name == '1.8 G8 Problem':
            
            self.__doc__ = """
            1.8 G8 Problem

            Dimensions: 2
            """

            def obj(x):
                y = -(np.sin(2*np.pi*x[0])**3*np.sin(2*np.pi*x[1]))/(x[0]**3*(x[0]+x[1]))
                return y

            def cns(x):
                g1 = x[0]**2-x[1]+1
                g2 = 1-x[0]+(x[1]-4)**2
                return [g1, g2]
            
            self.obj = obj
            self.cns = cns
            self.lb = [0, 0]
            self.ub = [10, 10]
            self.xopt = [1.2279713, 4.2453733]
            self.fopt = -0.095825

        elif name == '1.9 G9 Problem':
            
            self.__doc__ = """
            1.9 G9 Problem

            Dimensions: 7
            """

            def obj(x):
                y = (x[0]-10)**2+5*(x[1]-12)**2+x[2]**4+3*(x[3]-11)**2+ \
                    10*x[4]**6+7*x[5]**2+x[6]**4-4*x[5]*x[6]-10*x[5]-8*x[6]
                return y

            def cns(x):
                v1 = 2*x[0]**2;
                v2 = x[1]**2;
                g1 = v1+3*v2**2+x[2]+4*x[3]**2+5*x[4]-127;
                g2 = 7*x[0]+3*x[1]+10*x[2]**2+x[3]-x[4]-282;
                g3 = 23*x[0]+v2+6*x[5]**2-8*x[6]-196;
                g4 = 2*v1+v2-3*x[0]*x[1]+2*x[2]**2+5*x[5]-11*x[6];
                return [g1, g2, g3, g4]
            
            self.obj = obj
            self.cns = cns
            self.lb = (np.ones(7)*-10.0).tolist()
            self.ub = (np.ones(7)*10.0).tolist()
            self.xopt = [2.330499, 1.951372, -0.4775414, 4.365726, -0.6244870, 1.038131, 1.594227]
            self.fopt = 680.6300573

        elif name == '1.10 G10 Problem':
            
            self.__doc__ = """
            1.10 G10 Problem

            Dimensions: 8
            """

            def obj(x):
                y = x[0]+x[1]+x[2]
                return y

            def cns(x):
                g1 = -1+0.0025*(x[3]+x[5])
                g2 = -1+0.0025*(-x[3]+x[4]+x[6])
                g3 = -1+0.01*(-x[4]+x[7])
                g4 = 100*x[0]-x[0]*x[5]+833.33252*x[3]-83333.333
                g5 = x[1]*x[3]-x[1]*x[6]-1250*x[3]+1250*x[4]
                g6 = x[2]*x[4]-x[2]*x[7]-2500*x[4]+1250000
                return [g1, g2, g3, g4, g5, g6]
            
            self.obj = obj
            self.cns = cns
            self.lb = [100, 1000, 1000, 10, 10, 10, 10, 10]
            self.ub = [10000, 10000, 10000, 1000, 1000, 1000, 1000, 1000]
            self.xopt = [579.3167, 1359.943, 5110.071, 182.0174, 295.5985, 217.9799, 286.4162, 395.5979]
            self.fopt = 7049.3307

        elif name == '2.1 ALKYLATION':
            
            self.__doc__ = """
            2.1 ALKYLATION

            Dimensions: 5
            """

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
            self.fopt = -1768.75

        elif name == '2.2 CAMEL':

            self.__doc__ = """
            2.2 CAMEL

            Dimensions: 2

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
            self.fopt = 0.29861

        elif name == '2.3 FUNC2D':

            self.__doc__ = """
            2.3 FUNC2D

            Dimensions: 2
            """

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

        elif name == '2.4 GOLDPR':

            self.__doc__ = """
            2.4 GOLDPR

            Dimensions: 2

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

        elif name == '2.5 GOMEZ':
            
            self.__doc__ = """
            2.5 GOMEZ

            Dimensions: 2
            """

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

        elif name == '2.6 HS23':
            
            self.__doc__ = """
            2.6 HS23

            Dimensions: 2
            """

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

        elif name == '2.8 KS224':
            
            self.__doc__ = """
            2.8 KS224

            Dimensions: 2

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
            self.fopt = -304.0

        elif name == '2.9 KS250':
            
            self.__doc__ = """
            2.9 KS250

            Dimensions: 3

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

        elif name == '2.10 KS346':

            self.__doc__ = """
            2.10 KS346

            Dimensions: 3

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

        elif name == '2.11 NEWBRANIN':

            self.__doc__ = """
            2.11 NEWBRANIN

            Dimensions: 2

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

        elif name == '2.12 PRES':

            self.__doc__ = """
            2.12 PRES
            
            Dimensions: 2
            """

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
            self.fopt = 5.1766
        
        elif name == '':
            
            self.__doc__ = """
            """

            def obj(x):
                return y

            def cns(x):
                return []
            
            self.obj = obj
            self.cns = cns
            self.lb = []
            self.ub = []
            self.xopt = None
            self.fopt = None

        else:
            raise "Unkown problem name."

    def plot(self):
        debug_plot.plot(self)

    def __str__(self):
        string = ''
        for line in self.__doc__.split('\n'):
            string += line.lstrip() + '\n'
        return string


if __name__ == '__main__':

    problem = Cons('2.4 GOLDPR')
    x0 = [0.0, 0.0]

    print('obj(x0) = {}'.format(problem.obj(x0)))
    print('cns(x0) = {}'.format(problem.cns(x0)))
    print('lb = {}'.format(problem.lb))
    print('ub = {}'.format(problem.ub))
    print('xopt = {}'.format(problem.xopt))
    print('fopt = {}'.format(problem.fopt))
