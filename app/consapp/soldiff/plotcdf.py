
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
import math

def cdf(x, plot=True, *args, **kwargs):
    x, y = sorted(x), np.arange(len(x)) / len(x)
    return [x, y]

def plot_cdf3_csv(sol_1_fname, label1, sol_2_fname, label2, sol_3_fname, label3, title):
    
    sol_1 = np.genfromtxt(sol_1_fname, delimiter=',', unpack=True)

    rms_1_2d = np.sqrt(sol_1[1]*sol_1[1] + sol_1[2]*sol_1[2])
    rms_1_3d = np.sqrt(sol_1[1]*sol_1[1] + sol_1[2]*sol_1[2]+ sol_1[3]*sol_1[3])
    average_1_2d = np.average(rms_1_2d)
    average_1_3d = np.average(rms_1_3d)
    cdf_1_2d = cdf(rms_1_2d)
    cdf_1_3d = cdf(rms_1_3d)

    sol_2 = np.genfromtxt(sol_2_fname, delimiter=',', unpack=True)

    rms_2_2d = np.sqrt(sol_2[1]*sol_2[1] + sol_2[2]*sol_2[2])
    rms_2_3d = np.sqrt(sol_2[1]*sol_2[1] + sol_2[2]*sol_2[2] + sol_2[3]*sol_2[3])
    average_2_2d = np.average(rms_2_2d)
    average_2_3d = np.average(rms_2_3d)
    cdf_2_2d = cdf(rms_2_2d)
    cdf_2_3d = cdf(rms_2_3d)

    sol_3 = np.genfromtxt(sol_3_fname, delimiter=',', unpack=True)

    rms_3_2d = np.sqrt(sol_3[1]*sol_3[1] + sol_3[2]*sol_3[2])
    rms_3_3d = np.sqrt(sol_3[1]*sol_3[1] + sol_3[2]*sol_3[2] + sol_3[3]*sol_3[3])
    average_3_2d = np.average(rms_3_2d)
    average_3_3d = np.average(rms_3_3d)
    cdf_3_2d = cdf(rms_3_2d)
    cdf_3_3d = cdf(rms_3_3d)

    print("%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f" %(average_1_2d, average_1_3d,average_2_2d,average_2_3d,average_3_2d,average_3_3d))

    plt.plot(cdf_1_2d[0], cdf_1_2d[1], '-g', label=label1+' 2D')
    plt.plot(cdf_2_2d[0], cdf_2_2d[1], '-r', label=label2+' 2D')
    plt.plot(cdf_3_2d[0], cdf_3_2d[1], '-b', label=label3+' 2D')
    plt.plot(cdf_1_3d[0], cdf_1_3d[1], '-.g', label=label1+' 3D')
    plt.plot(cdf_2_3d[0], cdf_2_3d[1], '-.r', label=label2+' 3D')
    plt.plot(cdf_3_3d[0], cdf_3_3d[1], '-.b', label=label3+' 3D')

    plt.grid(True)
    plt.xlim(0, 5);
    plt.title('CDF '+title)
    plt.legend(loc="lower right")
    plt.show()

def plot_cdf2_csv(sol_1_fname, label1, sol_2_fname, label2, title):
    
    sol_1 = np.genfromtxt(sol_1_fname, delimiter=',', unpack=True)

    rms_1_2d = np.sqrt(sol_1[1]*sol_1[1] + sol_1[2]*sol_1[2])
    rms_1_3d = np.sqrt(sol_1[1]*sol_1[1] + sol_1[2]*sol_1[2]+ sol_1[3]*sol_1[3])
    average_1_2d = np.average(rms_1_2d)
    average_1_3d = np.average(rms_1_3d)
    cdf_1_2d = cdf(rms_1_2d)
    cdf_1_3d = cdf(rms_1_3d)

    sol_2 = np.genfromtxt(sol_2_fname, delimiter=',', unpack=True)

    rms_2_2d = np.sqrt(sol_2[1]*sol_2[1] + sol_2[2]*sol_2[2])
    rms_2_3d = np.sqrt(sol_2[1]*sol_2[1] + sol_2[2]*sol_2[2] + sol_2[3]*sol_2[3])
    average_2_2d = np.average(rms_2_2d)
    average_2_3d = np.average(rms_2_3d)
    cdf_2_2d = cdf(rms_2_2d)
    cdf_2_3d = cdf(rms_2_3d)

    print("%10.3f,%10.3f,%10.3f,%10.3f"%(average_1_2d, average_1_3d,average_2_2d,average_2_3d))

    plt.plot(cdf_1_2d[0], cdf_1_2d[1], '-g', label=label1+' 2D')
    plt.plot(cdf_2_2d[0], cdf_2_2d[1], '-r', label=label2+' 2D')
    plt.plot(cdf_1_3d[0], cdf_1_3d[1], '-.g', label=label1+' 3D')
    plt.plot(cdf_2_3d[0], cdf_2_3d[1], '-.r', label=label2+' 3D')

    plt.grid(True)
    plt.xlim(0, 5);
    plt.title('CDF '+title)
    plt.legend(loc="lower right")
    plt.show()

def plot_cdf1_csv(sol_1_fname, label1, title):
    
    sol_1 = np.genfromtxt(sol_1_fname, delimiter=',', unpack=True)

    rms_1_2d = np.sqrt(sol_1[1]*sol_1[1] + sol_1[2]*sol_1[2])
    rms_1_3d = np.sqrt(sol_1[1]*sol_1[1] + sol_1[2]*sol_1[2]+ sol_1[3]*sol_1[3])
    average_1_2d = np.average(rms_1_2d)
    average_1_3d = np.average(rms_1_3d)
    cdf_1_2d = cdf(rms_1_2d)
    cdf_1_3d = cdf(rms_1_3d)

    print("%10.3f,%10.3f" % (average_1_2d, average_1_3d))

    plt.plot(cdf_1_2d[0], cdf_1_2d[1], '-g', label=label1+' 2D')
    plt.plot(cdf_1_3d[0], cdf_1_3d[1], '-.g', label=label1+' 3D')

    plt.grid(True)
    plt.xlim(0, 5);
    plt.title('CDF '+title)
    plt.legend(loc="lower right")
    plt.show()