"""
Non-constrained optimization problem collections, based on:
https://www.sfu.ca/~ssurjano/optimization.html

1. Many Local Minima:
Ackley
Bukin N. 6
Cross-in-Tray
Drop-Wave
Eggholder
Gramacy & Lee 2012
Griewank
Holder Table
Langermann
Levy
Levy N. 13
Rastrigin
Schaffer N. 2
Schaffer N. 4
Schwefel
Shubert

4. Valley-Shaped
Six-Hump Camel

6. Other
Goldstein-Price
HARTMANN 3-D

"""

import math
import numpy as np

if __package__:
    from . import debug_plot
else:
    import debug_plot

class NonCons:

    """Non-constrained optimization problem.

    Args:
        name (str): Problem's name.

    Attributes:
        obj (function): Objective function.
        cns (None):
        lb (list): Lower bound of variables.
        ub (list): Upper bound of variables.
        xopt (list): Solution's variables.
        fopt (float): Solution's objective value.
    """

    def __init__(self, name):

        if name == 'Ackley':

            self.__doc__ = """
            The Ackley function is widely used for testing optimization algorithms.
            In its two-dimensional form, as shown in the plot above, it is characterized
            by a nearly flat outer region, and a large hole at the centre. The function
            poses a risk for optimization algorithms, particularly hillclimbing
            algorithms, to be trapped in one of its many local minima.

            1. Adorio, E. P., & Diliman, U. P. MVF - Multivariate Test Functions Library
            in C for Unconstrained Global Optimization (2005). Retrieved June 2013,
            from http://http://www.geocities.ws/eadorio/mvf.pdf.
            2. Molga, M., & Smutnicki, C. Test functions for optimization needs (2005).
            Retrieved June 2013, from http://www.zsd.ict.pwr.wroc.pl/files/docs/functions.pdf.
            3. Back, T. (1996). Evolutionary algorithms in theory and practice: evolution
            strategies, evolutionary programming, genetic algorithms. Oxford University
            Press on Demand.
            """

            def obj(x):
                d = 2.0
                c = 2*math.pi
                b = 0.2
                a = 20
                sum1 = 0
                sum2 = 0
                for xi in x:
                    sum1 = sum1 + xi**2.0
                    sum2 = sum2 + math.cos(c*xi)
                term1 = - a * math.exp(-b*math.sqrt(sum1/d))
                term2 = - math.exp(sum2/d)
                y = term1 + term2 + a + math.exp(1)
                return y

            self.obj = obj
            self.cns = None
            self.lb = [-40.0, -40.0]
            self.ub = [40.0, 40.0]
            self.xopt = [0.0, 0.0]
            self.fopt = 0.0

        elif name == 'Bukin N. 6':

            self.__doc__ = """
            The sixth Bukin function has many local minima, all of which lie in a ridge.

            Global Optimization Test Functions Index. Retrieved June 2013,
            from http://infinity77.net/global_optimization/test_functions.html#test-functions-index.
            """

            def obj(x):
                x1 = x[0]
                x2 = x[1]
                term1 = 100.0 * math.sqrt(abs(x2 - 0.01*x1**2))
                term2 = 0.01 * abs(x1 + 10.0)
                y = term1 + term2
                return y

            self.obj = obj
            self.cns = None
            self.lb = [-15, -3]
            self.ub = [-5, 3]
            self.xopt = [-10, 1]
            self.fopt = 0.0

        elif name == 'Cross-in-Tray':

            self.__doc__ = """
            The Cross-in-Tray function has multiple global minima. It is shown here
            with a smaller domain in the second plot, so that its characteristic
            "cross" will be visible. 

            Test functions for optimization. In Wikipedia. Retrieved June 2013,
            from https://en.wikipedia.org/wiki/Test_functions_for_optimization.
            """

            def obj(x):
                x1 = x[0]
                x2 = x[1]
                fact1 = math.sin(x1)*math.sin(x2)
                fact2 = math.exp(abs(100.0 - math.sqrt(x1**2+x2**2)/math.pi))
                y = - 0.0001 * (abs(fact1*fact2)+1.0)**0.1
                return y

            self.obj = obj
            self.cns = None
            self.lb = [-10, -10]
            self.ub = [10, 10]
            self.xopt = [1.3491, -1.3491] # others [1.3491, 1.3491] [-1.3491, 1.3491] [-1.3491, -1.3491]
            self.fopt = -2.06261

        elif name == 'Drop-Wave':

            self.__doc__ = """
            The Drop-Wave function is multimodal and highly complex. The second
            plot above shows the function on a smaller input domain, to
            illustrate its characteristic features. 

            Global Optimization Test Functions Index. Retrieved June 2013,
            from http://infinity77.net/global_optimization/test_functions.html#test-functions-index.
            """

            def obj(x):
                x1 = x[0]
                x2 = x[1]
                frac1 = 1.0 + math.cos(12.0*math.sqrt(x1**2+x2**2))
                frac2 = 0.5*(x1**2+x2**2) + 2.0
                y = - frac1/frac2
                return y

            self.obj = obj
            self.cns = None
            self.lb = [-2, -2]
            self.ub = [2, 2]
            self.xopt = [0, 0]
            self.fopt = -1.0

        elif name == 'Eggholder':

            self.__doc__ = """
            The Eggholder function is a difficult function to optimize,
            because of the large number of local minima. 

            Global Optimization Test Functions Index. Retrieved June 2013,
            from http://infinity77.net/global_optimization/test_functions.html#test-functions-index.
            """

            def obj(x):
                x1 = x[0]
                x2 = x[1]
                term1 = -(x2+47.0) * math.sin(math.sqrt(abs(x2+x1/2.0+47.0)))
                term2 = -x1 * math.sin(math.sqrt(abs(x1-(x2+47.0))))
                y = term1 + term2
                return y

            self.obj = obj
            self.cns = None
            self.lb = [-600, -600]
            self.ub = [600, 600]
            self.xopt = [512, 404.2319]
            self.fopt = -959.6407

        elif name == 'Gramacy & Lee 2012':

            self.__doc__ = """
            This is a simple one-dimensional test function. 

            1. Gramacy, R. B., & Lee, H. K. (2012). Cases for the nugget in
            modeling computer experiments. Statistics and Computing, 22(3), 713-722.
            2. Ranjan, P. (2013). Comment: EI Criteria for Noisy Computer Simulators.
            Technometrics, 55(1), 24-28.
            """

            def obj(x):
                term1 = math.sin(10.0*math.pi*x) / (2.0*x)
                term2 = (x-1.0)**4
                f = term1 + term2
                return f

            self.obj = obj
            self.cns = None
            self.lb = [0.5]
            self.ub = [2.5]
            self.xopt = None
            self.fopt = None

        elif name == 'Griewank':

            self.__doc__ = """
            The Griewank function has many widespread local minima, which are regularly
            distributed. The complexity is shown in the zoomed-in plots.

            1. Global Optimization Test Problems. Retrieved June 2013, from
            http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO.htm.
            2. Molga, M., & Smutnicki, C. Test functions for optimization needs (2005).
            Retrieved June 2013, from http://www.zsd.ict.pwr.wroc.pl/files/docs/functions.pdf.
            """

            def obj(x):
                d = 2
                total = 0
                prod = 1
                for ii in range(1, 3):
                    xi = x[ii-1]
                    total = total + (xi**2)/4000.0
                    prod = prod * math.cos(xi/math.sqrt(ii))
                f = total - prod + 1
                return f

            self.obj = obj
            self.cns = None
            self.lb = [-10, -10]
            self.ub = [10, 10]
            self.xopt = [0, 0]
            self.fopt = 0.0

        elif name == 'Holder Table':

            self.__doc__ = """
            The Holder Table function has many local minima, with four global minima. 

            1. Global Optimization Test Functions Index. Retrieved June 2013, from
            http://infinity77.net/global_optimization/test_functions.html#test-functions-index.
            2. Test functions for optimization. In Wikipedia. Retrieved June 2013, from
            https://en.wikipedia.org/wiki/Test_functions_for_optimization.
            """

            def obj(x):
                x1 = x[0]
                x2 = x[1]
                fact1 = math.sin(x1)*math.cos(x2)
                fact2 = math.exp(abs(1.0 - math.sqrt(x1**2+x2**2)/math.pi))
                f = - abs(fact1*fact2)
                return f

            self.obj = obj
            self.cns = None
            self.lb = [-10, -10]
            self.ub = [10, 10]
            self.xopt = [8.05502, 9.66459] # others [8.05502, -9.66459] [-8.05502, 9.66459] [-8.05502, -9.66459]
            self.fopt = -19.2085

        elif name == 'Langermann':

            self.__doc__ = """
            The Langermann function is multimodal, with many unevenly distributed local minima.

            1. Adorio, E. P., & Diliman, U. P. MVF - Multivariate Test Functions Library in C for
            Unconstrained Global Optimization (2005). Retrieved June 2013, from
            http://http://www.geocities.ws/eadorio/mvf.pdf.
            2. Molga, M., & Smutnicki, C. Test functions for optimization needs (2005). Retrieved
            June 2013, from http://www.zsd.ict.pwr.wroc.pl/files/docs/functions.pdf.
            """

            def obj(x):
                d = 2.0
                m = 5
                c = [1, 2, 5, 2, 3]
                A = np.array([[3, 5], [5, 2], [2, 1], [1, 4], [7, 9]])
                outer = 0
                for ii in range(5):
                    inner = 0
                    for jj in range(2):
                        xj = x[jj]
                        Aij = A[ii, jj]
                        inner = inner + (xj-Aij)**2
                    new = c[ii] * math.exp(-inner/math.pi) * math.cos(math.pi*inner)
                    outer = outer + new
                f = outer
                return f

            self.obj = obj
            self.cns = None
            self.lb = [0, 0]
            self.ub = [10, 10]
            self.xopt = None
            self.fopt = None

        elif name == 'Levy':

            self.__doc__ = """
            The function is usually evaluated on the hypercube.

            1. Global Optimization Test Functions Index. Retrieved June 2013, from
            http://infinity77.net/global_optimization/test_functions.html#test-functions-index.
            2. Laguna, M., & Marti, R. Experimental Testing of Advanced Scatter Search 
            Designs for Global Optimization of Multimodal Functions (2002). Retrieved June 
            2013, from http://www.uv.es/rmarti/paper/docs/global1.pdf.
            """

            def obj(x):
                d = 2
                w = []
                for ii in range(2):
                    w.append(1.0 + (x[ii] - 1.0)/4.0)
                term1 = (math.sin(math.pi*w[0]))**2
                term3 = (w[d-1]-1.0)**2 * (1.0+(math.sin(2*math.pi*w[d-1]))**2)
                total = 0
                for ii in range(1):
                    wi = w[ii]
                    new = (wi-1.0)**2 * (1.0+10.0*(math.sin(math.pi*wi+1))**2)
                    total = total + new
                f = term1 + total + term3
                return f

            self.obj = obj
            self.cns = None
            self.lb = [-10, -10]
            self.ub = [10, 10]
            self.xopt = [1, 1]
            self.fopt = 0.0

        elif name == 'Levy N. 13':

            self.__doc__ = """
            The function is usually evaluated on the square.

            Global Optimization Test Functions Index. Retrieved June 2013, from
            http://infinity77.net/global_optimization/test_functions.html#test-functions-index.
            """

            def obj(x):
                x1 = x[0]
                x2 = x[1]
                term1 = (math.sin(3*math.pi*x1))**2
                term2 = (x1-1.0)**2 * (1+(math.sin(3*math.pi*x2))**2)
                term3 = (x2-1.0)**2 * (1+(math.sin(2*math.pi*x2))**2)
                f = term1 + term2 + term3;
                return f

            self.obj = obj
            self.cns = None
            self.lb = [-10, -10]
            self.ub = [10, 10]
            self.xopt = [1, 1]
            self.fopt = 0.0

        elif name == 'Rastrigin':

            self.__doc__ = """
            Dimensions: d
            The Rastrigin function has several local minima. It is highly multimodal,
            but locations of the minima are regularly distributed. It is shown in the
            plot above in its two-dimensional form. 

            1. Global Optimization Test Problems. Retrieved June 2013, from
            http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO.htm.
            2. Pohlheim, H. GEATbx Examples: Examples of Objective Functions (2005).
            Retrieved June 2013, from http://www.geatbx.com/download/GEATbx_ObjFunExpl_v37.pdf.
            """

            def obj(x):
                d = 2
                total = 0
                for ii in range(2):
                    xi = x[ii]
                    total = total + (xi**2 - 10.0*math.cos(2.0*math.pi*xi))
                f = 10.0*d + total
                return f

            self.obj = obj
            self.cns = None
            self.lb = [-5, -5]
            self.ub = [5, 5]
            self.xopt = [0, 0]
            self.fopt = 0.0

        elif name == 'Schaffer N. 2':

            self.__doc__ = """
            The second Schaffer function.

            Test functions for optimization. In Wikipedia. Retrieved June 2013, from
            https://en.wikipedia.org/wiki/Test_functions_for_optimization.
            """

            def obj(x):
                x1 = x[0]
                x2 = x[1]
                fact1 = (math.sin(x1**2-x2**2))**2 - 0.5
                fact2 = (1.0 + 0.001*(x1**2+x2**2))**2
                f = 0.5 + fact1/fact2;
                return f

            self.obj = obj
            self.cns = None
            self.lb = [-5, -5]
            self.ub = [5, 5]
            self.xopt = [0, 0]
            self.fopt = 0.0

        elif name == 'Schaffer N. 4':

            self.__doc__ = """
            The fourth Schaffer function.

            Test functions for optimization. In Wikipedia. Retrieved June 2013, from
            https://en.wikipedia.org/wiki/Test_functions_for_optimization.
            """

            def obj(x):
                x1 = x[0]
                x2 = x[1]
                fact1 = math.cos(math.sin(abs(x1**2-x2**2))) - 0.5
                fact2 = (1.0 + 0.001*(x1**2+x2**2))**2
                f = 0.5 + fact1/fact2
                return f

            self.obj = obj
            self.cns = None
            self.lb = [-50, -50]
            self.ub = [50, 50]
            self.xopt = [0, 0]
            self.fopt = 0.0

        elif name == 'Schwefel':

            self.__doc__ = """
            Dimensions: d 
            The Schwefel function is complex, with many local minima.

            1. GEATbx: Examples of Objective Functions. Retrieved September 2014, from
            http://www.pg.gda.pl/~mkwies/dyd/geadocu/fcnfun7.html.
            2. Global Optimization Test Functions Index. Retrieved June 2013, from
            http://infinity77.net/global_optimization/test_functions.html#test-functions-index.
            3. Global Optimization Test Problems. Retrieved June 2013, from
            http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO.htm.
            4. Laguna, M., & Marti, R. Experimental Testing of Advanced Scatter Search Designs
            for Global Optimization of Multimodal Functions (2002). Retrieved June 2013, from
            http://www.uv.es/rmarti/paper/docs/global1.pdf.
            """

            def obj(x):
                d = 2
                total = 0
                for ii in range(2):
                    xi = x[ii]
                    total = total + xi*math.sin(math.sqrt(abs(xi)))
                f = 418.9829*d - total
                return f

            self.obj = obj
            self.cns = None
            self.lb = [-500, -500]
            self.ub = [500, 500]
            self.xopt = [420.9687, 420.9687]
            self.fopt = 0

        elif name == 'Shubert':
            
            self.__doc__ = """
            The Shubert function has several local minima and many global minima.

            Global Optimization Test Problems. Retrieved June 2013, from
            http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO.htm
            """

            def obj(x):
                x1 = x[0]
                x2 = x[1]
                sum1 = 0
                sum2 = 0
                for ii in range(1,6):
                    new1 = ii * math.cos((ii+1)*x1+ii)
                    new2 = ii * math.cos((ii+1)*x2+ii)
                    sum1 = sum1 + new1
                    sum2 = sum2 + new2
                y = sum1 * sum2
                return y

            self.obj = obj
            self.cns = None
            self.lb = [-10, -10]
            self.ub = [10, 10]
            self.xopt = None
            self.fopt = -186.7309

        elif name == '':
            
            self.__doc__ = """
            """

            def obj(x):
                return f

            self.obj = obj
            self.cns = None
            self.lb = []
            self.ub = []
            self.xopt = None
            self.fopt = None

        elif name == 'Six-Hump Camel':

            self.__doc__ = """
            6-Hump-Camel
            """

            def obj(x):
                f = (4.0-2.1*x[0]**2 + (x[0]**4.0)/3.0)*x[0]**2 + x[0]*x[1] + (-4.0 + 4.0*x[1]**2)*x[1]**2
                return f

            self.obj = obj
            self.cns = None
            self.lb = [-3, -2]
            self.ub = [3, 2]
            self.xopt = [0.0898, -0.7126] # anohter = [-0.0989, 0.7126]
            self.fopt = -1.0316

        elif name == 'Goldstein-Price':

            self.__doc__ = """
            The Goldstein-Price function has several local minima. 

            1. Dixon, L. C. W., & Szego, G. P. (1978). The global optimization
            problem: an introduction. Towards global optimization, 2, 1-15.
            2. Molga, M., & Smutnicki, C. Test functions for optimization needs (2005).
            Retrieved June 2013, from http://www.zsd.ict.pwr.wroc.pl/files/docs/functions.pdf
            3. Picheny, V., Wagner, T., & Ginsbourger, D. (2012). A benchmark of
            kriging-based infill criteria for noisy optimization.
            """

            def obj(x):
                g1 = 19.0 - 14.0*x[0] + 3.0*x[0]**2 - 14.0*x[1] + 6.0*x[0]*x[1] + 3.0*x[1]**2
                g2 = 18.0 - 32.0*x[0] + 12.0*x[0]**2 + 48.0*x[1] - 36.0*x[0]*x[1] + 27.0*x[1]**2
                f = (1.0 + ((x[0] + x[1] + 1.0)**2)*g1) * (30.0 + ((2.0*x[0] - 3.0*x[1])**2)*g2)
                f = math.log(f)
                return f

            self.obj = obj
            self.cns = None
            self.lb = [-2, -2]
            self.ub = [2, 2]
            self.xopt = [0, -1]
            self.fopt = 3.0

        elif name == 'HARTMANN 3-D':

            self.__doc__ = """
            The 3-dimensional Hartmann function has 4 local minima.

            Global Optimization Test Problems. Retrieved June 2013, from
            http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO.htm
            """

            def obj(x):
                alpha = [1.0, 1.2, 3.0, 3.2]
                A = np.array([
                    [3.0, 10, 30],
                    [0.1, 10, 35],
                    [3.0, 10, 30],
                    [0.1, 10, 35]])
                P = 10.0**(-4.0) * np.array([
                    [3689, 1170, 2673],
                    [4699, 4387, 7470],
                    [1091, 8732, 5547],
                    [381, 5743, 8828]])
                outer = 0.0
                for ii in range(4):
                    inner = 0
                    for jj in range(3):
                        xj = x[jj]
                        Aij = A[ii, jj]
                        Pij = P[ii, jj]
                        inner = inner + Aij*(xj-Pij)**2
                    new = alpha[ii] * math.exp(-inner)
                    outer = outer + new
                f = - outer
                return f

            self.obj = obj
            self.cns = None
            self.lb = [0, 0, 0]
            self.ub = [1, 1, 1]
            self.xopt = [0.114614, 0.555649, 0.852547]
            self.fopt = -3.86278

        else:
            raise "Unkown problem name."

    def plot(self):
        debug_plot.plot(self)

    def __str__(self):
        return self.__doc__

# --

if __name__ == '__main__':

    problem = NonCons('GOLDPR')
    x0 = [0.0, 0.0]

    print('obj(x0) = {}'.format(problem.obj(x0)))
    print('lb = {}'.format(problem.lb))
    print('ub = {}'.format(problem.ub))
    print('xopt = {}'.format(problem.xopt))
    print('fopt = {}'.format(problem.fopt))
