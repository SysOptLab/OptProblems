from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

# -- problem

def plot(problem):

    lb = problem.lb
    ub = problem.ub
    obj = problem.obj
    cns = problem.cns
    xopt = problem.xopt
    fopt = problem.fopt

    # -- verification

    try:
        print('obj(xopt) = {}'.format(obj(xopt)))
    except:
        print('obj(xopt) = {}'.format(obj(xopt[0])))

    try:
        print('xopt = {}'.format(xopt))
    except:
        print('xopt = {}'.format(xopt[0]))

    # -- plot

    if len(lb) == 1:

        # -- 1d line

        num = 100
        x_plot = np.linspace(lb, ub, num)
        fs = [obj(x) for x in x_plot]
        f_plot = np.array(fs).reshape(num)

        plt.figure()
        plt.plot(x_plot, f_plot)
        plt.xlabel('$x$')
        plt.ylabel('$f$')

    elif len(lb) == 2:

        # -- 2d coutour

        num = 40
        x1 = np.linspace(lb[0], ub[0], num)
        x2 = np.linspace(lb[1], ub[1], num)
        x1_plot, x2_plot = np.meshgrid(x1, x2)
        xs = np.hstack((x1_plot.reshape(-1,1), x2_plot.reshape(-1,1)))
        fs = [obj(x) for x in xs]
        f_plot = np.array(fs).reshape(num, num)

        if cns is not None:
            gs = [cns(x) for x in xs]
            dim = len(gs[0]) if isinstance(gs[0], list) else 1
            g_plot = np.array(gs).reshape(num, num, dim)

        plt.figure() # figsize=(12.80, 10.24)
        plt.xlabel('$x_1$')
        plt.ylabel('$x_2$')
        plt.pcolor(x1_plot, x2_plot, f_plot, cmap=plt.cm.jet) # plt.pcolormesh(x1_plot, x2_plot, f_plot)
        plt.colorbar()

        if cns is not None:
            for i in range(dim):
                plt.contour(x1_plot, x2_plot, g_plot[:,:,i], colors="k", levels=0, linestyles='dotted', lw=0.1)

        if isinstance(xopt[0], list):
            for xopt_i in xopt:
                plt.plot(xopt_i[0], xopt_i[1], '*r')
        else:
            plt.plot(xopt[0], xopt[1], '*r')

        plt.axis( [ lb[0], ub[0], lb[1], ub[1] ] )
        plt.axis('equal')
        plt.xlim(lb[0], ub[0])
        plt.ylim(lb[1], ub[1])

        # -- 3d surface

        fig = plt.figure()
        ax = fig.gca(projection='3d')
        surf = ax.plot_surface(x1_plot, x2_plot, f_plot, cmap=cm.jet, linewidth=0, antialiased=False, alpha=0.4)
        if isinstance(xopt[0], list):
            for xopt_i in xopt:
                ax.scatter(xopt_i[0], xopt_i[1], obj(xopt_i), color="r", marker='*', s=50)
        else:
            ax.scatter(xopt[0], xopt[1], obj(xopt), color="r", marker='*', s=50)
        
        ax.set_xlabel('$x_1$')
        ax.set_ylabel('$x_2$')
        ax.set_zlabel('$f$')
        ax.view_init(50, 235)

    else:
        print('Problem\'s dimensions is greater than 2.')
        pass

    plt.show()