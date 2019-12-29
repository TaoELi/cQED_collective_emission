"""
Plot the real space distribution of E-field intensity as a function of time
"""

import numpy as np
import columnplots as clp
from matplotlib.pyplot import FuncFormatter
import sys

#path = sys.argv[-1]
path="EM_FDTD_Full/"

def get_data(filename):
    return np.loadtxt(filename)

def running_mean(x, N):
    return np.convolve(x, np.ones((N,))/N, mode='valid')

data_QM = get_data(path+"Ez2_QM.txt")
data_EhR = get_data(path+"Ez2_MultiEh.txt")

n = 11
Lcavity = 2.0 * np.pi

spacing_y = np.max(np.max(data_QM[1:-1,:])) * 1.15

fig, ax = clp.initialize(col=1, row=1, width=4.3, height=6.8, sharex=True,
                    return_fig_args=True, sharey=True, fontsize=12, LaTeX=True)


N = 50 #average over 50 points
for i in range(n):
    X = [data_EhR[:,0], data_QM[:, 0]]
    X = [x / np.pi for x in X]
    Y = [data_EhR[:, i+1], data_QM[:, i+1]]
    Y = [(y + spacing_y*i) / spacing_y for y in Y]
    X_avg = [X[0][int(N/2)-1:-int(N/2)]]
    Y_avg = [running_mean(Y[0], N)]
    X = X + X_avg
    Y = Y + Y_avg
    clp.plotone(X, Y, ax, colors=['r', "k--", "c-."], labels=["MMST", "QM", "coarse-grained MMST"] if i == 0 else None, alphaspacing=0.01, ylim=[-0.2, 12])

ax.set_xlabel("Cavity position")
ax.set_ylabel("E-field Intensity (arbi. units)")

def format_func(value, tick_number):
    # find number of multiples of pi/2
    N = value
    if N == 0:
        return "0"
    else:
        return r"${0}\pi$".format(N)

ax.xaxis.set_major_formatter(FuncFormatter(format_func))

# Add arrow
ax.text(2.25, 5.0, 'time increment', rotation=90, color='b', size=12)
ax.annotate('', xy=(1.05, 0.0), xycoords='axes fraction', xytext=(1.05, 1.0),
            arrowprops=dict(arrowstyle="<-", color='b'))
#ax.grid(True, linestyle='--')
clp.adjust(tight_layout=True, savefile="SE_intensity_MMST.pdf")
