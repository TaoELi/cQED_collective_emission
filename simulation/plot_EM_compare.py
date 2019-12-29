"""
Plot the real space distribution of E-field intensity as a function of time
"""

import numpy as np
import columnplots as clp
from matplotlib.pyplot import FuncFormatter
import sys

path = sys.argv[-1]
#path="EM_FDTD_Full/"

def get_data(filename):
    return np.loadtxt(filename)

def running_mean(x, N):
    return np.convolve(x, np.ones((N,))/N, mode='valid')


n = 11
Lcavity = 2.0 * np.pi


def format_func(value, tick_number):
    # find number of multiples of pi/2
    N = value
    if N == 0:
        return "0"
    else:
        return r"${0}\pi$".format(N)
def format_func_small(value, tick_number):
    # find number of multiples of pi/2
    N = value
    if N == 0:
        return "0"
    else:
        return r"${0}\pi$".format(int(1000*N)/1000.0)

fig, axes = clp.initialize(col=2, row=4, width=4.3*3, height=6.3*2,
                    return_fig_args=True, sharey=True, fontsize=12, LaTeX=True, labelthem=True, labelthemPosition=[0.08, 0.98])


N = 50 #average over 50 points

paths = ["EM_FDTD_Full/", "EM_FDTD_Truncated/", "EM_MODE_Full/", "EM_MODE_Truncated/"]
methods = ["FDTD 400 modes", "FDTD 100 modes", "XP 400 modes", "XP 100 modes"]

for j in range(4):
    data_QM = get_data(paths[j]+"Ez2_QM.txt")
    data_EhR = get_data(paths[j]+"Ez2_MultiEh.txt")
    spacing_y = np.max(np.max(data_QM[1:-1,:])) * 1.15
    for i in range(n):
        X = [data_EhR[:,0], data_QM[:, 0]]
        X = [x / np.pi for x in X]
        Y = [data_EhR[:, i+1], data_QM[:, i+1]]
        Y = [(y + spacing_y*i) / spacing_y for y in Y]
        X_avg = [X[0][int(N/2)-1:-int(N/2)]]
        Y_avg = [running_mean(Y[0], N)]
        X = X + X_avg
        Y = Y + Y_avg
        for k in range(2):
            clp.plotone(X, Y, axes[k, j], colors=['r', "k--", "c-."], labels=["MMST", "QM", "coarse-grained MMST"] if (i == 0 and j == 0) else None, alphaspacing=0.01, ylim=[-0.2, 12], showlegend=True if (j == 0 and k == 0) else False, xlim=None if k == 0 else [0.9622, 1.0378])
    axes[0, j].set_title(methods[j])
    axes[1, j].set_xlabel("Cavity position")
    if j is 0:
        for k in range(2):
            axes[k, j].set_ylabel("E-field Intensity (arbi. units)")
    axes[0, j].xaxis.set_major_formatter(FuncFormatter(format_func))
    axes[1, j].xaxis.set_major_formatter(FuncFormatter(format_func))

# Add arrow
axes[0, -1].text(2.25, 5.0, 'time increment', rotation=90, color='b', size=12)
#axes[1, -1].text(2.25, 5.0, 'time increment', rotation=90, color='b', size=12)
for k in range(2):
    axes[k, -1].annotate('', xy=(1.05, 0.0), xycoords='axes fraction', xytext=(1.05, 1.0),
            arrowprops=dict(arrowstyle="<-", color='b'))
#ax.grid(True, linestyle='--')
clp.adjust(tight_layout=True, savefile="EM_compare_MMST.pdf")
