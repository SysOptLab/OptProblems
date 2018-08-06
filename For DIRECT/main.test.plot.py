import os
import shutil
import math
import matplotlib.pyplot as plt

# -- problem

problems_cons = [
    "1.4 G4 Problem",
    "1.6 G6 Problem",
    "1.8 G8 Problem",
    "2.2 CAMEL",
    "2.3 FUNC2D",
    "2.4 GOLDPR",
    "2.5 GOMEZ",
    "2.6 HS23",
    "2.8 KS224",
    "2.9 KS250",
    "2.10 KS346",
    "2.11 NEWBRANIN",
    "2.12 PRES"]

problems_noncons = [
    "1.1 Ackley Function",
    "1.3 Cross-in-Tray Function",
    "1.4 Drop-Wave Function",
    "1.5 Eggholder Function",
    "1.6 Gramacy and Lee (2012) Function",
    "1.7 Griewank Function",
    "1.8 Holder Table Function",
    "1.10 Levy Function",
    "1.11 Levy Function N. 13",
    "1.12 Rastrigin Function",
    "1.13 Schaffer Function N. 2",
    "1.14 Schaffer Function N. 4",
    "1.15 Schwefel Function",
    "1.16 Shubert Function",
    "2.1 Bohachevsky Function",
    "2.2 Perm Function",
    "2.3 Rotated Hyper-Ellipsoid Function",
    "2.5 Sum of Different Powers Function",
    "2.6 Sum Squares Function",
    "2.7 Trid Function",
    "3.1 Booth Function",
    "3.2 Matyas Function",
    "3.3 McCormick Function",
    "3.5 Zakharov Function",
    "4.1 Three-Hump Camel Function",
    "4.2 Six-Hump Camel Function",
    "4.3 Dixon-Price Function",
    "4.4 Rosenbrock Function",
    "5.2 Easom Function",
    "5.3 Michalewicz Function",
    "6.1 Beale Function",
    "6.3 Colville Function",
    "6.2 Branin Function",
    "6.4 Forrester et al. (2008) Function",
    "6.5 Goldstein-Price Function",
    "6.6 Hartmann 3-D Function",
    "6.7 Hartmann 4-D Function",
    "6.8 Hartmann 6-D Function",
    "6.9 Perm Function",
    "6.11 Shekel Function 5",
    "6.12 Shekel Function 7",
    "6.13 Shekel Function 10",
    "6.14 Styblinski-Tang Function"]

# -- optimization

settings = [
    "ordinary",
    "penalty-long",
    "surr-opt-org",
    "surr-opt-rbf",
    "penalty-long_surr-opt-rbf",
    "surr-one-rbf",
    "penalty-long_surr-one-rbf",
    ]

colors = ['b', 'g', 'r', 'm', 'y', 'k']
linestyles = [':', ':', ':', ':', ':', ':']

# -- plot result and compare


def get_result(path):
    fcount = []
    fbest = []
    with open(path, 'r') as fileIO:
        for line_string in fileIO:
            line_list = line_string.split()
            fcount.append(int(line_list[1]))
            fbest.append(float(line_list[3]))
    return fcount, fbest

for name in problems_cons:
    root = "./20180726/"
    
    plt.figure(figsize=(12.80, 10.24))

    for setting, color, linestyle in zip(settings, colors, linestyles):
        path = root + setting + "/cons-" + name + '.txt'
        fcount, fbest = get_result(path)
        plt.plot(fcount, fbest, color, label=setting, linestyle=linestyle)
        for i, f in enumerate(fbest):
            if f != math.inf:
                plt.scatter(fcount[i], fbest[i], c=color, marker='o', s=20)
                break

    plt.xlabel('fcount')
    plt.ylabel('fbest')
    # plt.xlim([0, 200])
    plt.title(name)
    plt.legend()

    mng = plt.get_current_fig_manager()
    mng.resize(*mng.window.maxsize())

    plt.savefig(fname=root + '/cons-{}.png'.format(name))
    plt.close()

for name in problems_noncons:
    root = "./20180726/"
    
    plt.figure(figsize=(12.80, 10.24))

    for setting, color, linestyle in zip(settings, colors, linestyles):
        path = root + setting + "/noncons-" + name + '.txt'
        fcount, fbest = get_result(path)
        plt.plot(fcount, fbest, color, label=setting, linestyle=linestyle)
        for i, f in enumerate(fbest):
            if f != math.inf:
                plt.scatter(fcount[i], fbest[i], c=color, marker='o', s=20)
                break

    plt.xlabel('fcount')
    plt.ylabel('fbest')
    # plt.xlim([0, 200])
    plt.title(name)
    plt.legend()

    mng = plt.get_current_fig_manager()
    mng.resize(*mng.window.maxsize())

    plt.savefig(fname=root + '/noncons-{}.png'.format(name))
    plt.close()