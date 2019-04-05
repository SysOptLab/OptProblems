"""
Non-constrained optimization problem collections.

1. Many Local Minima

1.1 Ackley Function
1.2 Bukin Function N. 6
1.3 Cross-in-Tray Function
1.4 Drop-Wave Function
1.5 Eggholder Function
1.6 Gramacy and Lee (2012) Function
1.7 Griewank Function
1.8 Holder Table Function
1.10 Levy Function
1.11 Levy Function N. 13
1.12 Rastrigin Function
1.13 Schaffer Function N. 2
1.14 Schaffer Function N. 4
1.15 Schwefel Function
1.16 Shubert Function

2. Bowl-Shaped

2.1 Bohachevsky Function
2.2 Perm Function
2.3 Rotated Hyper-Ellipsoid Function
2.4 Sphere Function Modified
2.5 Sum of Different Powers Function
2.6 Sum Squares Function
2.7 Trid Function

3. Plate-Shaped

3.1 Booth Function
3.2 Matyas Function
3.3 McCormick Function
3.5 Zakharov Function

4. Valley-Shaped

4.1 Three-Hump Camel Function
4.2 Six-Hump Camel Function
4.3 Dixon-Price Function
4.4 Rosenbrock Function

5. Steep Ridges/Drops

5.2 Easom Function
5.3 Michalewicz Function

6. Other

6.1 Beale Function
6.2 Branin Function
6.3 Colville Function
6.4 Forrester et al. (2008) Function
6.5 Goldstein-Price Function
6.6 Hartmann 3-D Function
6.7 Hartmann 4-D Function
6.8 Hartmann 6-D Function
6.9 Perm Function
6.11 Shekel Function 5
6.12 Shekel Function 7
6.13 Shekel Function 10
6.14 Styblinski-Tang Function

"""

import numpy as np

if __package__:
    from . import debug_plot
else:
    import debug_plot

class NonCons:

    """Non-constrained optimization problem

    Args:
        name (str): problem's name

    Attributes:
        obj (func): objfunction
        cns (None):
        lb (list): lower bound of variables
        ub (list): upper bound of variables
        xopt (list): solution's variables
        fopt (float): solution's obj value
    """

    def __init__(self, name, dimensions=2):

        if name == '1.1 Ackley Function':

            self.__doc__ = """
            1.1 Ackley Function

            Dimensions: d

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
                d = dimensions
                a = 20.0
                b = 0.2
                c = 2.0*np.pi
                sum1 = 0.0
                sum2 = 0.0
                for xi in x:
                    sum1 = sum1 + xi**2.0
                    sum2 = sum2 + np.cos(c*xi)
                term1 = - a * np.exp(-b*np.sqrt(sum1/d))
                term2 = - np.exp(sum2/d)
                y = term1 + term2 + a + np.exp(1)
                return y

            self.obj = obj
            self.cns = None
            self.lb = [-39.0, -39.0] # original bound = -40
            self.ub = [40.0, 40.0]
            self.xopt = [0.0, 0.0]
            self.fopt = 0.0

        elif name == '1.2 Bukin Function N. 6':

            self.__doc__ = """
            1.2 Bukin Function N. 6

            Dimensions: 2

            The sixth Bukin function has many local minima, all of which lie in a ridge.

            Global Optimization Test Functions Index. Retrieved June 2013, from
            http://infinity77.net/global_optimization/test_functions.html#test-functions-index.
            """

            def obj(x):
                x1 = x[0]
                x2 = x[1]
                term1 = 100.0 * np.sqrt(abs(x2 - 0.01*x1**2))
                term2 = 0.01 * abs(x1 + 10.0)
                y = term1 + term2
                return y

            self.obj = obj
            self.cns = None
            self.lb = [-15, -3]
            self.ub = [-5, 3]
            self.xopt = [-10, 1]
            self.fopt = 0.0

        elif name == '1.3 Cross-in-Tray Function':

            self.__doc__ = """
            1.3 Cross-in-Tray Function

            Dimensions: 2

            The Cross-in-Tray function has multiple global minima. It is shown here
            with a smaller domain in the second plot, so that its characteristic
            "cross" will be visible. 

            Test functions for optimization. In Wikipedia. Retrieved June 2013,
            from https://en.wikipedia.org/wiki/Test_functions_for_optimization.
            """

            def obj(x):
                x1 = x[0]
                x2 = x[1]
                fact1 = np.sin(x1)*np.sin(x2)
                fact2 = np.exp(abs(100.0 - np.sqrt(x1**2+x2**2)/np.pi))
                y = - 0.0001 * (abs(fact1*fact2)+1.0)**0.1
                return y

            self.obj = obj
            self.cns = None
            self.lb = [-10, -10]
            self.ub = [10, 10]
            self.xopt = [[1.3491, -1.3491], [1.3491, 1.3491], [-1.3491, 1.3491], [-1.3491, -1.3491]]
            self.fopt = -2.06261

        elif name == '1.4 Drop-Wave Function':

            self.__doc__ = """
            1.4 Drop-Wave Function

            Dimensions: 2

            The Drop-Wave function is multimodal and highly complex. The second
            plot above shows the function on a smaller input domain, to
            illustrate its characteristic features. 

            Global Optimization Test Functions Index. Retrieved June 2013,
            from http://infinity77.net/global_optimization/test_functions.html#test-functions-index.
            """

            def obj(x):
                x1 = x[0]
                x2 = x[1]
                frac1 = 1.0 + np.cos(12.0*np.sqrt(x1**2+x2**2))
                frac2 = 0.5*(x1**2+x2**2) + 2.0
                y = - frac1/frac2
                return y

            self.obj = obj
            self.cns = None
            self.lb = [-1.9, -1.9] # original bound = -2.0
            self.ub = [2, 2]
            self.xopt = [0, 0]
            self.fopt = -1.0

        elif name == '1.5 Eggholder Function':

            self.__doc__ = """
            1.5 Eggholder Function

            Dimensions: 2

            The Eggholder function is a difficult function to optimize,
            because of the large number of local minima. 

            Global Optimization Test Functions Index. Retrieved June 2013,
            from http://infinity77.net/global_optimization/test_functions.html#test-functions-index.
            """

            def obj(x):
                x1 = x[0]
                x2 = x[1]
                term1 = -(x2+47.0) * np.sin(np.sqrt(abs(x2+x1/2.0+47.0)))
                term2 = -x1 * np.sin(np.sqrt(abs(x1-(x2+47.0))))
                y = term1 + term2
                return y

            self.obj = obj
            self.cns = None
            self.lb = [-600, -600]
            self.ub = [600, 600]
            self.xopt = [512, 404.2319]
            self.fopt = -959.6407

        elif name == '1.6 Gramacy and Lee (2012) Function':

            self.__doc__ = """
            1.6 Gramacy and Lee (2012) Function
            
            Dimensions: 1

            This is a simple one-dimensional test function. 

            1. Gramacy, R. B., & Lee, H. K. (2012). Cases for the nugget in
            modeling computer experiments. Statistics and Computing, 22(3), 713-722.
            2. Ranjan, P. (2013). Comment: EI Criteria for Noisy Computer Simulators.
            Technometrics, 55(1), 24-28.
            """

            def obj(x):
                x = np.array(x)
                term1 = np.sin(10.0*np.pi*x) / (2.0*x)
                term2 = (x-1.0)**4
                f = term1[0] + term2[0]
                return f

            self.obj = obj
            self.cns = None
            self.lb = [0.5]
            self.ub = [2.5]
            self.xopt = [0.54856368]
            self.fopt = -0.8690111349647177

        elif name == '1.7 Griewank Function':

            self.__doc__ = """
            1.7 Griewank Function
            
            Dimensions: d

            The Griewank function has many widespread local minima, which are regularly
            distributed. The complexity is shown in the zoomed-in plots.

            1. Global Optimization Test Problems. Retrieved June 2013, from
            http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO.htm.
            2. Molga, M., & Smutnicki, C. Test functions for optimization needs (2005).
            Retrieved June 2013, from http://www.zsd.ict.pwr.wroc.pl/files/docs/functions.pdf.
            """

            def obj(x):
                d = dimensions
                total = 0
                prod = 1
                for ii in range(1, d+1):
                    xi = x[ii-1]
                    total = total + (xi**2)/4000.0
                    prod = prod * np.cos(xi/np.sqrt(ii))
                f = total - prod + 1
                return f

            self.obj = obj
            self.cns = None
            self.lb = [-9, -9] # original bound = -10.0
            self.ub = [10, 10]
            self.xopt = [0, 0]
            self.fopt = 0.0

        elif name == '1.8 Holder Table Function':

            self.__doc__ = """
            1.8 Holder Table Function
            
            Dimensions: 2

            The Holder Table function has many local minima, with four global minima. 

            1. Global Optimization Test Functions Index. Retrieved June 2013, from
            http://infinity77.net/global_optimization/test_functions.html#test-functions-index.
            2. Test functions for optimization. In Wikipedia. Retrieved June 2013, from
            https://en.wikipedia.org/wiki/Test_functions_for_optimization.
            """

            def obj(x):
                x1 = x[0]
                x2 = x[1]
                fact1 = np.sin(x1)*np.cos(x2)
                fact2 = np.exp(abs(1.0 - np.sqrt(x1**2+x2**2)/np.pi))
                f = - abs(fact1*fact2)
                return f

            self.obj = obj
            self.cns = None
            self.lb = [-10, -10]
            self.ub = [10, 10]
            self.xopt = [[8.05502, 9.66459], [8.05502, -9.66459], [-8.05502, 9.66459], [-8.05502, -9.66459]]
            self.fopt = -19.2085

        elif name == '1.10 Levy Function':

            self.__doc__ = """
            1.10 Levy Function
            
            Dimensions: d

            The function is usually evaluated on the hypercube.

            1. Global Optimization Test Functions Index. Retrieved June 2013, from
            http://infinity77.net/global_optimization/test_functions.html#test-functions-index.
            2. Laguna, M., & Marti, R. Experimental Testing of Advanced Scatter Search 
            Designs for Global Optimization of Multimodal Functions (2002). Retrieved June 
            2013, from http://www.uv.es/rmarti/paper/docs/global1.pdf.
            """

            def obj(x):
                d = dimensions
                w = []
                for ii in range(2):
                    w.append(1.0 + (x[ii] - 1.0)/4.0)
                term1 = (np.sin(np.pi*w[0]))**2
                term3 = (w[d-1]-1.0)**2 * (1.0+(np.sin(2*np.pi*w[d-1]))**2)
                total = 0
                for ii in range(1):
                    wi = w[ii]
                    new = (wi-1.0)**2 * (1.0+10.0*(np.sin(np.pi*wi+1))**2)
                    total = total + new
                f = term1 + total + term3
                return f

            self.obj = obj
            self.cns = None
            self.lb = [-9, -9] # origin bound = -10.0
            self.ub = [10, 10]
            self.xopt = [1, 1]
            self.fopt = 0.0

        elif name == '1.11 Levy Function N. 13':

            self.__doc__ = """
            1.11 Levy Function N. 13

            Dimensions: 2 

            The function is usually evaluated on the square.

            Global Optimization Test Functions Index. Retrieved June 2013, from
            http://infinity77.net/global_optimization/test_functions.html#test-functions-index.
            """

            def obj(x):
                x1 = x[0]
                x2 = x[1]
                term1 = (np.sin(3*np.pi*x1))**2
                term2 = (x1-1.0)**2 * (1+(np.sin(3*np.pi*x2))**2)
                term3 = (x2-1.0)**2 * (1+(np.sin(2*np.pi*x2))**2)
                f = term1 + term2 + term3;
                return f

            self.obj = obj
            self.cns = None
            self.lb = [-10, -10]
            self.ub = [10, 10]
            self.xopt = [1, 1]
            self.fopt = 0.0

        elif name == '1.12 Rastrigin Function':

            self.__doc__ = """
            1.12 Rastrigin Function

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
                d = dimensions
                total = 0
                for xi in x:
                    total = total + (xi**2 - 10.0*np.cos(2.0*np.pi*xi))
                f = 10.0*d + total
                return f

            self.obj = obj
            self.cns = None
            self.lb = [-4, -4] # original bound = -5.0
            self.ub = [5, 5]
            self.xopt = [0, 0]
            self.fopt = 0.0

        elif name == '1.13 Schaffer Function N. 2':

            self.__doc__ = """
            1.13 Schaffer Function N. 2

            Dimensions: 2

            The second Schaffer function.

            Test functions for optimization. In Wikipedia. Retrieved June 2013, from
            https://en.wikipedia.org/wiki/Test_functions_for_optimization.
            """

            def obj(x):
                x1 = x[0]
                x2 = x[1]
                fact1 = (np.sin(x1**2-x2**2))**2 - 0.5
                fact2 = (1.0 + 0.001*(x1**2+x2**2))**2
                f = 0.5 + fact1/fact2;
                return f

            self.obj = obj
            self.cns = None
            self.lb = [-4, -4] # original bound = -5.0
            self.ub = [5, 5]
            self.xopt = [0, 0]
            self.fopt = 0.0

        elif name == '1.14 Schaffer Function N. 4':

            self.__doc__ = """
            1.14 Schaffer Function N. 4

            Dimensions: 2

            The fourth Schaffer function.

            Test functions for optimization. In Wikipedia. Retrieved June 2013, from
            https://en.wikipedia.org/wiki/Test_functions_for_optimization.
            """

            def obj(x):
                x1 = x[0]
                x2 = x[1]
                fact1 = np.cos(np.sin(abs(x1**2-x2**2))) - 0.5
                fact2 = (1.0 + 0.001*(x1**2+x2**2))**2
                f = 0.5 + fact1/fact2
                return f

            self.obj = obj
            self.cns = None
            self.lb = [-50, -50]
            self.ub = [50, 50]
            self.xopt = [0, 0]
            self.fopt = 0.0

        elif name == '1.15 Schwefel Function':

            self.__doc__ = """
            1.15 Schwefel Function
            
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
                d = dimensions
                total = 0
                for ii in range(d):
                    xi = x[ii]
                    total = total + xi*np.sin(np.sqrt(abs(xi)))
                f = 418.9829*d - total
                return f

            self.obj = obj
            self.cns = None
            self.lb = [-500, -500]
            self.ub = [500, 500]
            self.xopt = [420.9687, 420.9687]
            self.fopt = 0

        elif name == '1.16 Shubert Function':
            
            self.__doc__ = """
            1.16 Shubert Function
            
            Dimensions: 2

            The Shubert function has several local minima and many global minima.

            Global Optimization Test Problems. Retrieved June 2013, from
            http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO.htm
            """

            def obj(x):
                x1 = x[0]
                x2 = x[1]
                sum1 = 0
                sum2 = 0
                for ii in range(1, 6):
                    new1 = ii * np.cos((ii+1)*x1+ii)
                    new2 = ii * np.cos((ii+1)*x2+ii)
                    sum1 = sum1 + new1
                    sum2 = sum2 + new2
                y = sum1 * sum2
                return y

            self.obj = obj
            self.cns = None
            self.lb = [-10, -10]
            self.ub = [10, 10]
            self.xopt = [-1.425128, -0.800273]
            self.fopt = -186.7309

        elif name == '2.1 Bohachevsky Function':
            
            self.__doc__ = """
            2.1 Bohachevsky Function
            
            Dimensions: 2 

            The Bohachevsky functions all have the same similar bowl shape. The one shown above is
            the first function.

            Global Optimization Test Problems. Retrieved June 2013, from
            http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO.htm.
            """

            def obj(x):
                x1 = x[0]
                x2 = x[1]

                term1 = x1**2
                term2 = 2*x2**2
                term3 = -0.3 * np.cos(3*np.pi*x1)
                term4 = -0.4 * np.cos(4*np.pi*x2)

                y = term1 + term2 + term3 + term4 + 0.7
                return y

            self.obj = obj
            self.cns = None
            self.lb = [-99.0, -99.0] # original bound = -100.0
            self.ub = [100.0, 100.0]
            self.xopt = [0.0, 0.0]
            self.fopt = 0.0

        elif name == '2.2 Perm Function':
            
            self.__doc__ = """
            2.2 Perm Function
            
            Dimensions: d

            The function is usually evaluated on the hypercube xi [-d, d], for all i = 1, ..., d. 

            Global Optimization Test Problems. Retrieved June 2013, from
            http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO.htm.
            """

            def obj(x):
                b = 10
                d = dimensions
                outer = 0

                for ii in range(1, d+1):
                    inner = 0
                    for jj in range(1, d+1):
                        xj = x[jj-1]
                        inner = inner + (jj+b)*(xj**ii-(1.0/jj)**ii)
                    outer = outer + inner**2

                y = outer;
                return y

            self.obj = obj
            self.cns = None
            self.lb = (np.ones(dimensions)*-2.0).tolist()
            self.ub = (np.ones(dimensions)*2.0).tolist()
            self.xopt = [1, 1.0/2.0] # [..., 1/3, 1/4 ..., 1/d]
            self.fopt = 0.0

        elif name == '2.3 Rotated Hyper-Ellipsoid Function':
            
            self.__doc__ = """
            2.3 Rotated Hyper-Ellipsoid Function
            
            Dimensions: d 

            The Rotated Hyper-Ellipsoid function is continuous, convex and unimodal. It is an
            extension of the Axis Parallel Hyper-Ellipsoid function, also referred to as the
            Sum Squares function. The plot shows its two-dimensional form.

            Molga, M., & Smutnicki, C. Test functions for optimization needs (2005). Retrieved
            June 2013, from http://www.zsd.ict.pwr.wroc.pl/files/docs/functions.pdf.
            """

            def obj(x):
                d = dimensions
                outer = 0

                for ii in range(1, d+1):
                    inner = 0
                    for jj in range(1, ii+1):
                        xj = x[jj-1]
                        inner = inner + xj**2
                    outer = outer + inner

                y = outer
                return y

            self.obj = obj
            self.cns = None
            self.lb = (np.ones(dimensions)*-59.0).tolist()
            self.ub = (np.ones(dimensions)*60.0).tolist()
            self.xopt = [0.0, 0.0] # [..., 0.0]
            self.fopt = 0.0

        elif name == '2.4 Sphere Function Modified':
            
            self.__doc__ = """
            2.4 Sphere Function Modified
            
            Dimensions: 6

            The Sphere function has d local minima except for the global one. It is continuous,
            convex and unimodal. The plot shows its two-dimensional form.

            1. Dixon, L. C. W., & Szego, G. P. (1978). The global optimization problem: an
            introduction. Towards global optimization, 2, 1-15.
            2. Molga, M., & Smutnicki, C. Test functions for optimization needs (2005). Retrieved
            June 2013, from http://www.zsd.ict.pwr.wroc.pl/files/docs/functions.pdf.
            3. Picheny, V., Wagner, T., & Ginsbourger, D. (2012). A benchmark of kriging-based
            infill criteria for noisy optimization.
            """

            def obj(x):
                d = 6
                sum = 0;
                for ii in range(1, d+1):
                    xi = x[ii-1]
                    sum = sum + (xi**2)*(2**ii)
                y = (sum - 1745.0) / 899.0
                return y

            self.obj = obj
            self.cns = None
            self.lb = (np.ones(6)*-0.9).tolist() # original bound = -1.0
            self.ub = (np.ones(6)*1.0).tolist()
            self.xopt = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            self.fopt = 0.0

        elif name == '2.5 Sum of Different Powers Function':
            
            self.__doc__ = """
            2.5 Sum of Different Powers Function
            
            Dimensions: d 

            The Sum of Different Powers function is unimodal. It is shown here in its
            two-dimensional form.

            Molga, M., & Smutnicki, C. Test functions for optimization needs (2005). Retrieved
            June 2013, from http://www.zsd.ict.pwr.wroc.pl/files/docs/functions.pdf.
            """

            def obj(x):
                d = dimensions
                sum = 0

                for ii in range(1, d+1):
                    xi = x[ii-1]
                    new = (abs(xi))**(ii+1.0)
                    sum = sum + new

                y = sum
                return y 

            self.obj = obj
            self.cns = None
            self.lb = (np.ones(dimensions)*-0.9).tolist() # original bound = -1.0
            self.ub = (np.ones(dimensions)*1.0).tolist()
            self.xopt = [0.0, 0.0] # [..., 0.0]
            self.fopt = 0.0

        elif name == '2.6 Sum Squares Function':
            
            self.__doc__ = """
            2.6 Sum Squares Function
            
            Dimensions: d 

            The Sum Squares function, also referred to as the Axis Parallel Hyper-Ellipsoid
            function, has no local minimum except the global one. It is continuous, convex
            and unimodal. It is shown here in its two-dimensional form. 

            1. Global Optimization Test Problems. Retrieved June 2013, from
            http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO.htm.

            2. Molga, M., & Smutnicki, C. Test functions for optimization needs (2005).
            Retrieved June 2013, from http://www.zsd.ict.pwr.wroc.pl/files/docs/functions.pdf.
            """

            def obj(x):
                d = dimensions
                sum = 0
                for ii in range(1, d+1):
                    xi = x[ii-1]
                    sum = sum + ii*xi**2
                y = sum
                return y

            self.obj = obj
            self.cns = None
            self.lb = (np.ones(dimensions)*-9.0).tolist() # original bound = -10.0
            self.ub = (np.ones(dimensions)*10.0).tolist()
            self.xopt = [0.0, 0.0] # [..., 0.0]
            self.fopt = 0.0

        elif name == '2.7 Trid Function':
            
            self.__doc__ = """
            2.7 Trid Function
            
            Dimensions: d 

            The Trid function has no local minimum except the global one. It is shown here in its
            two-dimensional form.

            1. Adorio, E. P., & Diliman, U. P. MVF - Multivariate Test Functions Library in C for
            Unconstrained Global Optimization (2005). Retrieved August 2017, from
            http://www.geocities.ws/eadorio/mvf.pdf.
            2. Global Optimization Test Problems. Retrieved June 2013, from
            http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO.htm.
            """

            def obj(x):
                d = dimensions
                sum1 = (x[0]-1.0)**2
                sum2 = 0

                for ii in range(2, d+1):
                    xi = x[ii-1]
                    xold = x[ii-2]
                    sum1 = sum1 + (xi-1.0)**2
                    sum2 = sum2 + xi*xold

                y = sum1 - sum2
                return y

            self.obj = obj
            self.cns = None
            self.lb = (np.ones(dimensions)*-0.9*dimensions**2).tolist() # original bound = -dimensions
            self.ub = (np.ones(dimensions)*dimensions**2).tolist()
            self.xopt = [1*(dimensions+1-1), 2*(dimensions+1-2)]
            self.fopt = - dimensions*(dimensions+4.0)*(dimensions-1.0)/6.0

        elif name == '3.1 Booth Function':
            
            self.__doc__ = """
            3.1 Booth Function
            
            Dimensions: 2 

            Global Optimization Test Problems. Retrieved June 2013, from
            http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO.htm.
            """

            def obj(x):
                x1 = x[0]
                x2 = x[1]

                term1 = (x1 + 2.0*x2 - 7)**2
                term2 = (2.0*x1 + x2 - 5)**2

                y = term1 + term2
                return y

            self.obj = obj
            self.cns = None
            self.lb = [-10.0, -10.0]
            self.ub = [10.0, 10.0]
            self.xopt = [1.0, 3.0]
            self.fopt = 0.0

        elif name == '3.2 Matyas Function':
            
            self.__doc__ = """
            3.2 Matyas Function
            
            Dimensions: 2 

            The Matyas function has no local minima except the global one.

            Global Optimization Test Problems. Retrieved June 2013, from
            http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO.htm.
            """

            def obj(x):
                x1 = x[0]
                x2 = x[1]
                term1 = 0.26 * (x1**2 + x2**2)
                term2 = -0.48*x1*x2
                y = term1 + term2
                return y

            self.obj = obj
            self.cns = None
            self.lb = [-9.0, -9.0] # original bound = -10.0
            self.ub = [10.0, 10.0]
            self.xopt = [0.0, 0.0]
            self.fopt = 0.0

        elif name == '3.3 McCormick Function':
            
            self.__doc__ = """
            3.3 McCormick Function
            
            Dimensions: 2 

            Adorio, E. P., & Diliman, U. P. MVF - Multivariate Test Functions Library in C for
            Unconstrained Global Optimization (2005). Retrieved June 2013, from
            http://http://www.geocities.ws/eadorio/mvf.pdf.
            """

            def obj(x):
                x1 = x[0]
                x2 = x[1]
                term1 = np.sin(x1 + x2)
                term2 = (x1 - x2)**2
                term3 = -1.5*x1
                term4 = 2.5*x2
                y = term1 + term2 + term3 + term4 + 1.0
                return y

            self.obj = obj
            self.cns = None
            self.lb = [-1.5, -3.0]
            self.ub = [4.0, 4.0]
            self.xopt = [-0.54719, -1.54719]
            self.fopt = -1.9133

        elif name == '3.5 Zakharov Function':
            
            self.__doc__ = """
            3.5 Zakharov Function
            
            Dimensions: d 

            The Zakharov function has no local minima except the global one. It is shown here
            in its two-dimensional form. 

            Global Optimization Test Problems. Retrieved June 2013, from
            http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO.htm.
            """

            def obj(x):
                d = dimensions
                sum1 = 0
                sum2 = 0

                for ii in range(1, d+1):
                    xi = x[ii-1]
                    sum1 = sum1 + xi**2
                    sum2 = sum2 + 0.5*ii*xi

                y = sum1 + sum2**2 + sum2**4
                return y

            self.obj = obj
            self.cns = None
            self.lb = (np.ones(dimensions)*-5.0).tolist()
            self.ub = (np.ones(dimensions)*10.0).tolist()
            self.xopt = [0.0, 0.0]
            self.fopt = 0.0

        elif name == '4.1 Three-Hump Camel Function':
            
            self.__doc__ = """
            4.1 Three-Hump Camel Function
            
            Dimensions: 2 

            The plot on the left shows the three-hump Camel function on its recommended input
            domain, and the plot on the right shows only a portion of this domain, to allow
            for easier viewing of the function's key characteristics. The function has three
            local minima. 

            Test functions for optimization. In Wikipedia. Retrieved June 2013, from
            https://en.wikipedia.org/wiki/Test_functions_for_optimization.
            """

            def obj(x):
                x1 = x[0]
                x2 = x[1]
                term1 = 2*x1**2
                term2 = -1.05*x1**4
                term3 = x1**6 / 6
                term4 = x1*x2
                term5 = x2**2
                y = term1 + term2 + term3 + term4 + term5
                return y

            self.obj = obj
            self.cns = None
            self.lb = [-4.0, -4.0] # original bound = -5.0
            self.ub = [5.0, 5.0]
            self.xopt = [0.0, 0.0]
            self.fopt = 0.0

        elif name == '4.2 Six-Hump Camel Function':

            self.__doc__ = """
            4.2 Six-Hump Camel Function
            
            Dimensions: 2 

            The plot on the left shows the six-hump Camel function on its recommended input domain,
            and the plot on the right shows only a portion of this domain, to allow for easier
            viewing of the function's key characteristics. The function has six local minima, two
            of which are global. 

            Molga, M., & Smutnicki, C. Test functions for optimization needs (2005). Retrieved June
            2013, from http://www.zsd.ict.pwr.wroc.pl/files/docs/functions.pdf.
            """

            def obj(x):
                y = (4.0-2.1*x[0]**2 + (x[0]**4.0)/3.0)*x[0]**2 + x[0]*x[1] + (-4.0 + 4.0*x[1]**2)*x[1]**2
                return y

            self.obj = obj
            self.cns = None
            self.lb = [-3, -2]
            self.ub = [3, 2]
            self.xopt = [[0.0898, -0.7126], [-0.0898, 0.7126]]
            self.fopt = -1.0316

        elif name == '4.3 Dixon-Price Function':
            
            self.__doc__ = """
            4.3 Dixon-Price Function
            
            Dimensions: d 

            Global Optimization Test Functions Index. Retrieved June 2013, from
            http://infinity77.net/global_optimization/test_functions.html#test-functions-index.
            """

            def obj(x):
                x1 = x[0]
                d = dimensions
                term1 = (x1-1)**2
                sum = 0

                for ii in range(2, d+1):
                    xi = x[ii-1]
                    xold = x[ii-2]
                    new = ii * (2*xi**2 - xold)**2
                    sum = sum + new

                y = term1 + sum
                return y

            self.obj = obj
            self.cns = None
            self.lb = (np.ones(dimensions)*-9.0).tolist() # original bound = -10.0
            self.ub = (np.ones(dimensions)*10.0).tolist()
            self.xopt = [2.0**(-(2**1-2)/(2**1)), 2.0**(-(2**2-2)/(2**2))] # [..., 2.0**(-(2**d-2)/(2**d))]
            self.fopt = 0.0

        elif name == '4.4 Rosenbrock Function':
            
            self.__doc__ = """
            4.4 Rosenbrock Function
            
            Dimensions: d 

            The Rosenbrock function, also referred to as the Valley or Banana function, is a
            popular test problem for gradient-based optimization algorithms. It is shown in
            the plot above in its two-dimensional form.

            The function is unimodal, and the global minimum lies in a narrow, parabolic valley.
            However, even though this valley is easy to find, convergence to the minimum is
            difficult (Picheny et al., 2012).

            1. Dixon, L. C. W., & Szego, G. P. (1978). The global optimization problem: an 
            introduction. Towards global optimization, 2, 1-15.
            2. Global Optimization Test Problems. Retrieved June 2013, from
            http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO.htm.
            3. Molga, M., & Smutnicki, C. Test functions for optimization needs (2005). 
            Retrieved June 2013, from http://www.zsd.ict.pwr.wroc.pl/files/docs/functions.pdf.
            4. Picheny, V., Wagner, T., & Ginsbourger, D. (2012). A benchmark of kriging-based
            infill criteria for noisy optimization.
            """

            def obj(x):
                d = dimensions
                sum = 0

                for ii in range(1, d):
                    xi = x[ii-1]
                    xnext = x[ii]
                    new = 100.0*(xnext-xi**2)**2 + (xi-1.0)**2
                    sum = sum + new

                y = sum
                return y

            self.obj = obj
            self.cns = None
            self.lb = (np.ones(dimensions)*-5.0).tolist()
            self.ub = (np.ones(dimensions)*10.0).tolist()
            self.xopt = (np.ones(dimensions)*1.0).tolist()
            self.fopt = 0.0

        elif name == '5.2 Easom Function':
            
            self.__doc__ = """
            5.2 Easom Function
            
            Dimensions: 2 

            The Easom function has several local minima. It is unimodal, and the global minimum has
            a small area relative to the search space. 

            Global Optimization Test Problems. Retrieved June 2013, from
            http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO.htm.
            """

            def obj(x):
                x1 = x[0]
                x2 = x[1]
                fact1 = -np.cos(x1)*np.cos(x2)
                fact2 = np.exp(-(x1-np.pi)**2-(x2-np.pi)**2)
                y = fact1*fact2
                return y

            self.obj = obj
            self.cns = None
            self.lb = [-100.0, -100.0]
            self.ub = [100.0, 100.0]
            self.xopt = [np.pi, np.pi]
            self.fopt = -1.0

        elif name == '5.3 Michalewicz Function':
            
            self.__doc__ = """
            5.3 Michalewicz Function
            
            Dimensions: d 

            The Michalewicz function has d! local minima, and it is multimodal. The parameter m defines
            the steepness of they valleys and ridges; a larger m leads to a more difficult search. The
            recommended value of m is m = 10.

            1. Global Optimization Test Functions Index. Retrieved June 2013, from
            http://infinity77.net/global_optimization/test_functions.html#test-functions-index.
            2. Global Optimization Test Problems. Retrieved June 2013, from
            http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO.htm.
            3. Molga, M., & Smutnicki, C. Test functions for optimization needs (2005).
            Retrieved June 2013, from http://www.zsd.ict.pwr.wroc.pl/files/docs/functions.pdf.
            """

            def obj(x):
                m = 10
                d = dimensions
                sum = 0

                for ii in range(1, d+1):
                    xi = x[ii-1]
                    new = np.sin(xi) * (np.sin(ii*xi**2/np.pi))**(2*m)
                    sum  = sum + new

                y = -sum
                return y

            self.obj = obj
            self.cns = None
            self.lb = np.zeros(dimensions).tolist()
            self.ub = (np.ones(dimensions)*np.pi).tolist()
            self.xopt = [2.20, 1.57]
            self.fopt = -1.8013

        elif name == '6.1 Beale Function':
            
            self.__doc__ = """
            6.1 Beale Function
            
            Dimensions: 2 

            The Beale function is multimodal, with sharp peaks at the corners of the input domain. 

            Global Optimization Test Problems. Retrieved June 2013, from
            http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO.htm.
            """

            def obj(x):
                x1 = x[0]
                x2 = x[1]
                term1 = (1.5 - x1 + x1*x2)**2
                term2 = (2.25 - x1 + x1*x2**2)**2
                term3 = (2.625 - x1 + x1*x2**3)**2
                y = term1 + term2 + term3
                return y

            self.obj = obj
            self.cns = None
            self.lb = [-4.5, -4.5]
            self.ub = [4.5, 4.5]
            self.xopt = [3.0, 0.5]
            self.fopt = 0.0

        elif name == '6.2 Branin Function':
            
            self.__doc__ = """
            6.2 Branin Function

            Dimensions: 2 

            The Branin, or Branin-Hoo, function has three global minima. The recommended values
            of a, b, c, r, s and t are: a = 1, b = 5.1 / (4pi2), c = 5 / pi, r = 6, s = 10 and
            t = 1 ~ (8pi).

            1. Dixon, L. C. W., & Szego, G. P. (1978). The global optimization problem: an introduction.
             Towards global optimization, 2, 1-15.
            2. Forrester, A., Sobester, A., & Keane, A. (2008). Engineering design via surrogate 
            modelling: a practical guide. Wiley.
            3. Global Optimization Test Problems. Retrieved June 2013, from
            http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO.htm.
            4. Molga, M., & Smutnicki, C. Test functions for optimization needs (2005). Retrieved June 
            2013, from http://www.zsd.ict.pwr.wroc.pl/files/docs/functions.pdf.
            5. Picheny, V., Wagner, T., & Ginsbourger, D. (2012). A benchmark of kriging-based infill 
            criteria for noisy optimization.
            """

            def obj(x):
                y = (x[1]-(5.1/(4*np.pi**2))*x[0]**2+5*x[0]/np.pi-6)**2+10*(1-1/(8*np.pi))*np.cos(x[0])+10;
                return y

            self.obj = obj
            self.cns = None
            self.lb = [-5, 0]
            self.ub = [10, 15]
            self.xopt = [[-np.pi, 12.275], [np.pi, 2.275], [9.42478, 2.475]]
            self.fopt = 0.3979

        elif name == '6.3 Colville Function':
            
            self.__doc__ = """
            6.3 Colville Function

            Dimensions: 4

            Global Optimization Test Problems. Retrieved June 2013, from
            http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO.htm.
            """

            def obj(x):
                x1 = x[0]
                x2 = x[1]
                x3 = x[2]
                x4 = x[3]
                term1 = 100 * (x1**2-x2)**2
                term2 = (x1-1)**2
                term3 = (x3-1)**2
                term4 = 90.0 * (x3**2-x4)**2
                term5 = 10.1 * ((x2-1)**2 + (x4-1)**2)
                term6 = 19.8*(x2-1)*(x4-1)
                y = term1 + term2 + term3 + term4 + term5 + term6
                return y

            self.obj = obj
            self.cns = None
            self.lb = (np.ones(4)*-9.0).tolist() # original bound = -10.0
            self.ub = (np.ones(4)*10.0).tolist()
            self.xopt = [1.0, 1.0, 1.0, 1.0]
            self.fopt = 0.0

        elif name == '6.4 Forrester et al. (2008) Function':
            
            self.__doc__ = """
            6.4 Forrester et al. (2008) Function

            Dimensions: 1 

            This function is a simple one-dimensional test function. It is multimodal, with one global
            minimum, one local minimum and a zero-gradient inflection point.
            
            Forrester, A., Sobester, A., & Keane, A. (2008). Engineering design via surrogate
            modelling: a practical guide. Wiley.
            """

            def obj(x):
                x = np.array(x)
                fact1 = (6*x - 2)**2;
                fact2 = np.sin(12*x - 4);
                y = fact1 * fact2
                return y[0]

            self.obj = obj
            self.cns = None
            self.lb = [0]
            self.ub = [1]
            self.xopt = [0.75724768]
            self.fopt = -6.0207400551464705

        elif name == '6.5 Goldstein-Price Function':

            self.__doc__ = """
            6.5 Goldstein-Price Function

            Dimensions: 2 

            The Goldstein-Price function has several local minima. 

            1. Dixon, L. C. W., & Szego, G. P. (1978). The global optimization
            problem: an introduction. Towards global optimization, 2, 1-15.
            2. Molga, M., & Smutnicki, C. Test functions for optimization needs (2005).
            Retrieved June 2013, from http://www.zsd.ict.pwr.wroc.pl/files/docs/functions.pdf
            3. Picheny, V., Wagner, T., & Ginsbourger, D. (2012). A benchmark of
            kriging-based infill criteria for noisy optimization.
            """

            def obj(x):
                x1 = x[0]
                x2 = x[1]
                fact1a = (x1 + x2 + 1)**2
                fact1b = 19 - 14*x1 + 3*x1**2 - 14*x2 + 6*x1*x2 + 3*x2**2
                fact1 = 1 + fact1a*fact1b
                fact2a = (2*x1 - 3*x2)**2
                fact2b = 18 - 32*x1 + 12*x1**2 + 48*x2 - 36*x1*x2 + 27*x2**2
                fact2 = 30 + fact2a*fact2b
                y = fact1*fact2
                return y

            self.obj = obj
            self.cns = None
            self.lb = [-2, -2]
            self.ub = [2, 2]
            self.xopt = [0, -1]
            self.fopt = 3.0

        elif name == '6.6 Hartmann 3-D Function':

            self.__doc__ = """
            6.6 Hartmann 3-D Function

            Dimensions: 3 

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
                    new = alpha[ii] * np.exp(-inner)
                    outer = outer + new
                f = - outer
                return f

            self.obj = obj
            self.cns = None
            self.lb = [0, 0, 0]
            self.ub = [1, 1, 1]
            self.xopt = [0.114614, 0.555649, 0.852547]
            self.fopt = -3.8628

        elif name == '6.7 Hartmann 4-D Function':
            
            self.__doc__ = """
            6.7 Hartmann 4-D Function

            Dimensions: 4

            The 4-dimensional Hartmann function is multimodal. It is given here in the form of
            Picheny et al. (2012), having a mean of zero and a variance of one. The authors also
            add a small Gaussian error term to the output. 

            1. Dixon, L. C. W., & Szego, G. P. (1978). The global optimization problem: an
            introduction. Towards global optimization, 2, 1-15.
            2. Global Optimization Test Problems. Retrieved June 2013, from
            http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO.htm.
            3. Picheny, V., Wagner, T., & Ginsbourger, D. (2012). A benchmark of kriging-based
            infill criteria for noisy optimization.
            """

            def obj(x):

                a = np.empty([4, 6])
                a[0,0]=10.0;	a[0,1]=3.0;		a[0,2]=17.0;	a[0,3]=3.5;		a[0,4]=1.7;		a[0,5]=8.0
                a[1,0]=0.05;	a[1,1]=10.0;	a[1,2]=17.0;	a[1,3]=0.1;		a[1,4]=8.0;		a[1,5]=14.0
                a[2,0]=3.0;		a[2,1]=3.5;		a[2,2]=1.7;		a[2,3]=10.0;	a[2,4]=17.0;	a[2,5]=8.0
                a[3,0]=17.0;	a[3,1]=8.0;		a[3,2]=0.05;	a[3,3]=10.0;	a[3,4]=0.1;		a[3,5]=14.0

                c = np.empty([4])
                c[0]=1.0;   c[1]=1.2;   c[2]=3.0;   c[3]=3.2

                p = np.empty([4, 6])
                p[0,0]=0.1312;	p[0,1]=0.1696;	p[0,2]=0.5569;	p[0,3]=0.0124;	p[0,4]=0.8283;	p[0,5]=0.5886
                p[1,0]=0.2329;	p[1,1]=0.4135;	p[1,2]=0.8307;	p[1,3]=0.3736;	p[1,4]=0.1004;	p[1,5]=0.9991
                p[2,0]=0.2348;	p[2,1]=0.1451;	p[2,2]=0.3522;	p[2,3]=0.2883;	p[2,4]=0.3047;	p[2,5]=0.6650
                p[3,0]=0.4047;	p[3,1]=0.8828;	p[3,2]=0.8732;	p[3,3]=0.5743;	p[3,4]=0.1091;	p[3,5]=0.0381

                s = 0

                for i in range(1, 5):
                    sm = 0
                    for j in range(1, 5):
                        sm = sm+a[i-1,j-1]*(x[j-1]-p[i-1,j-1])**2
                    s = s+c[i-1]*np.exp(-sm)
                s = 1.0/0.839 * (1.1 - s)
                y = s
                return y

            self.obj = obj
            self.cns = None
            self.lb = np.zeros(4).tolist()
            self.ub = np.ones(4).tolist()
            self.xopt = None
            self.fopt = -3.1335

        elif name == '6.8 Hartmann 6-D Function':
            
            self.__doc__ = """
            6.8 Hartmann 6-D Function

            Dimensions: 6 

            The 6-dimensional Hartmann function has 6 local minima. 

            1. Dixon, L. C. W., & Szego, G. P. (1978). The global optimization problem: an
            introduction. Towards global optimization, 2, 1-15.
            2. Global Optimization Test Problems. Retrieved June 2013, from
            http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO.htm.
            3. Picheny, V., Wagner, T., & Ginsbourger, D. (2012). A benchmark of kriging-based
            infill criteria for noisy optimization.
            """

            def obj(x):

                a = np.empty([4, 6])
                a[0,0]=10.0;	a[0,1]=3.0;		a[0,2]=17.0;	a[0,3]=3.5;		a[0,4]=1.7;		a[0,5]=8.0
                a[1,0]=0.05;	a[1,1]=10.0;	a[1,2]=17.0;	a[1,3]=0.1;		a[1,4]=8.0;		a[1,5]=14.0
                a[2,0]=3.0;		a[2,1]=3.5;		a[2,2]=1.7;		a[2,3]=10.0;	a[2,4]=17.0;	a[2,5]=8.0
                a[3,0]=17.0;	a[3,1]=8.0;		a[3,2]=0.05;	a[3,3]=10.0;	a[3,4]=0.1;		a[3,5]=14.0

                c = np.empty([4])
                c[0]=1.0;   c[1]=1.2;   c[2]=3.0;   c[3]=3.2

                p = np.empty([4, 6])
                p[0,0]=0.1312;	p[0,1]=0.1696;	p[0,2]=0.5569;	p[0,3]=0.0124;	p[0,4]=0.8283;	p[0,5]=0.5886
                p[1,0]=0.2329;	p[1,1]=0.4135;	p[1,2]=0.8307;	p[1,3]=0.3736;	p[1,4]=0.1004;	p[1,5]=0.9991
                p[2,0]=0.2348;	p[2,1]=0.1451;	p[2,2]=0.3522;	p[2,3]=0.2883;	p[2,4]=0.3047;	p[2,5]=0.6650
                p[3,0]=0.4047;	p[3,1]=0.8828;	p[3,2]=0.8732;	p[3,3]=0.5743;	p[3,4]=0.1091;	p[3,5]=0.0381

                s = 0

                for i in range(1, 5):
                    sm = 0
                    for j in range(1, 7):
                        sm = sm+a[i-1,j-1]*(x[j-1]-p[i-1,j-1])**2
                    s = s+c[i-1]*np.exp(-sm)

                y = -s
                return y

            self.obj = obj
            self.cns = None
            self.lb = np.zeros(6).tolist()
            self.ub = np.ones(6).tolist()
            self.xopt = [0.20169, 0.150011, 0.476874, 0.275332, 0.311652, 0.6573]
            self.fopt = -3.3224

        elif name == '6.9 Perm Function':
            
            self.__doc__ = """
            6.9 Perm Function

            Dimensions: d 

            The Perm d, beta function.

            Global Optimization Test Problems. Retrieved June 2013, from
            http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO.htm.
            """

            def obj(x):
                b = 0.5;
                d = dimensions
                outer = 0

                for ii in range(1, d+1):
                    inner = 0
                    for jj in range(1, d+1):
                        xj = x[jj-1]
                        inner = inner + (jj**ii+b)*((xj/jj)**ii-1.0)
                    outer = outer + inner**2

                y = outer
                return y

            self.obj = obj
            self.cns = None
            self.lb = (np.ones(dimensions)*-dimensions).tolist()
            self.ub = (np.ones(dimensions)*dimensions).tolist()
            self.xopt = [1, 2] # [3, 4, ...]
            self.fopt = 0.0

        elif name == '6.11 Shekel Function 5':
            
            self.__doc__ = """
            6.11 Shekel Function 5

            Dimensions: 4 

            The Shekel function has m local minima. Above are the recommended values of m, 
            the beta-vector and the C-matrix; beta is an m-dimensional vector, and C is a 4
            -by-m-dimensional matrix
            """

            def obj(x):
                m = 5
                b = 0.1 * np.array([1, 2, 2, 4, 4, 6, 3, 7, 5, 5])
                C = np.array([
                    [4.0, 1.0, 8.0, 6.0, 3.0, 2.0, 5.0, 8.0, 6.0, 7.0],
                    [4.0, 1.0, 8.0, 6.0, 7.0, 9.0, 5.0, 1.0, 2.0, 3.6],
                    [4.0, 1.0, 8.0, 6.0, 3.0, 2.0, 3.0, 8.0, 6.0, 7.0],
                    [4.0, 1.0, 8.0, 6.0, 7.0, 9.0, 3.0, 1.0, 2.0, 3.6]])
                outer = 0

                for ii in range(1, m+1):
                    bi = b[ii-1]
                    inner = 0
                    for jj in range(1, 5):
                        xj = x[jj-1]
                        Cji = C[jj-1, ii-1]
                        inner = inner + (xj-Cji)**2
                    outer = outer + 1/(inner+bi)

                y = -outer
                return y

            self.obj = obj
            self.cns = None
            self.lb = [0, 0, 0, 0]
            self.ub = [10, 10, 10, 10]
            self.xopt = [4, 4, 4, 4]
            self.fopt = -10.1532

        elif name == '6.12 Shekel Function 7':
            
            self.__doc__ = """
            6.12 Shekel Function 7

            Dimensions: 4 

            The Shekel function has m local minima. Above are the recommended values of m, 
            the beta-vector and the C-matrix; beta is an m-dimensional vector, and C is a 4
            -by-m-dimensional matrix

            1. Global Optimization Test Problems. Retrieved June 2013, from
            http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO.htm.

            2. Molga, M., & Smutnicki, C. Test functions for optimization needs (2005). Retrieved
            June 2013, from http://www.zsd.ict.pwr.wroc.pl/files/docs/functions.pdf.
            """

            def obj(x):
                m = 7
                b = 0.1 * np.array([1, 2, 2, 4, 4, 6, 3, 7, 5, 5])
                C = np.array([
                    [4.0, 1.0, 8.0, 6.0, 3.0, 2.0, 5.0, 8.0, 6.0, 7.0],
                    [4.0, 1.0, 8.0, 6.0, 7.0, 9.0, 5.0, 1.0, 2.0, 3.6],
                    [4.0, 1.0, 8.0, 6.0, 3.0, 2.0, 3.0, 8.0, 6.0, 7.0],
                    [4.0, 1.0, 8.0, 6.0, 7.0, 9.0, 3.0, 1.0, 2.0, 3.6]])
                outer = 0

                for ii in range(1, m+1):
                    bi = b[ii-1]
                    inner = 0
                    for jj in range(1, 5):
                        xj = x[jj-1]
                        Cji = C[jj-1, ii-1]
                        inner = inner + (xj-Cji)**2
                    outer = outer + 1/(inner+bi)

                y = -outer
                return y

            self.obj = obj
            self.cns = None
            self.lb = [0, 0, 0, 0]
            self.ub = [10, 10, 10, 10]
            self.xopt = [4, 4, 4, 4]
            self.fopt = -10.4029

        elif name == '6.13 Shekel Function 10':
            
            self.__doc__ = """
            6.13 Shekel Function 10

            Dimensions: 4 

            The Shekel function has m local minima. Above are the recommended values of m, 
            the beta-vector and the C-matrix; beta is an m-dimensional vector, and C is a 4
            -by-m-dimensional matrix
            """

            def obj(x):
                m = 10
                b = 0.1 * np.array([1, 2, 2, 4, 4, 6, 3, 7, 5, 5])
                C = np.array([
                    [4.0, 1.0, 8.0, 6.0, 3.0, 2.0, 5.0, 8.0, 6.0, 7.0],
                    [4.0, 1.0, 8.0, 6.0, 7.0, 9.0, 5.0, 1.0, 2.0, 3.6],
                    [4.0, 1.0, 8.0, 6.0, 3.0, 2.0, 3.0, 8.0, 6.0, 7.0],
                    [4.0, 1.0, 8.0, 6.0, 7.0, 9.0, 3.0, 1.0, 2.0, 3.6]])
                outer = 0

                for ii in range(1, m+1):
                    bi = b[ii-1]
                    inner = 0
                    for jj in range(1, 5):
                        xj = x[jj-1]
                        Cji = C[jj-1, ii-1]
                        inner = inner + (xj-Cji)**2
                    outer = outer + 1/(inner+bi)

                y = -outer
                return y

            self.obj = obj
            self.cns = None
            self.lb = [0, 0, 0, 0]
            self.ub = [10, 10, 10, 10]
            self.xopt = [4, 4, 4, 4]
            self.fopt = -10.5364

        elif name == '6.14 Styblinski-Tang Function':
            
            self.__doc__ = """
            6.14 Styblinski-Tang Function

            Dimensions: d 

            The Styblinski-Tang function is shown here in its two-dimensional form. 

            """

            def obj(x):
                d = dimensions
                sum = 0
                for ii in range(1, d+1):
                    xi = x[ii-1]
                    new = xi**4 - 16*xi**2 + 5*xi
                    sum = sum + new

                y = sum/2.0
                return y

            self.obj = obj
            self.cns = None
            self.lb = (np.ones(dimensions)*-5.0).tolist()
            self.ub = (np.ones(dimensions)*5.0).tolist()
            self.xopt = (np.ones(dimensions)*-2.903534).tolist()
            self.fopt = -39.16599*dimensions

        else:
            raise "Unkown problem name."

    def plot(self):
        debug_plot.plot(self)

    def __str__(self):
        string = ''
        for line in self.__doc__.split('\n'):
            string += line.lstrip() + '\n'
        return string

# --

if __name__ == '__main__':

    name = '6.5 Goldstein-Price Function'
    problem = NonCons(name)
    problem.plot()
