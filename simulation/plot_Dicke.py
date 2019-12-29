import numpy as np
import columnplots as clp
import glob
import sys

path=sys.argv[-1]

kFGR = np.pi

from matplotlib.pyplot import FuncFormatter

def get_popu_data(filename, N):
    data = np.loadtxt(filename)
    t, p1, p2 = data[:-1,0], data[:-1, 1], data[:-1, 2]
    for i in range(N-1):
        p2 += data[:-1, 2*(i+2)]
    return t, p2 / N


def gather_data(Nlst, path = "SE_Dicke/"):
    xs, ys = [], []
    for N in Nlst:
        subpath = path + "N_%d" %N
        t_MMST, p2_MMST = get_popu_data("%s/traj_MultiEh.txt" %subpath, N)
        xs.append(t_MMST)
        ys.append(p2_MMST)
    return xs, ys

Nlst = [1, 7, 35]
labels = ["N = %d" %N for N in Nlst]

axes = clp.initialize(col=2, row=1, width=4.3, height=4.3*0.618*2, fontsize=12, LaTeX=True, labelthem=True, labelthemPosition=[0.04, 0.96])

xs, ys = gather_data(Nlst)

clp.plotone([x/kFGR for x in xs], ys, axes[0], labels=labels, colors=["r--", "c", "b"], xlabel="time ($t$)", ylabel="Avg. excited state popu. $\\bar{\\rho}_{ee}(t)$", lw=2, alphaspacing=0.02)

# In the second subplot, we need to calculate the effective decay rate
Nlst = [7, 9, 11, 15, 19, 25, 35, 45]
xs, ys = gather_data(Nlst)


from scipy.optimize import curve_fit
def exponenial_func(x, k):
    return np.exp(-k*x)

y_fit = []
k_lst = []
for i in range(len(Nlst)):
    popt, pcov = curve_fit(exponenial_func, xs[i][0:2500], ys[i][0:2500], p0=(1.0))
    y_fit.append(exponenial_func(xs[i], *popt))
    k_lst.append(popt[0])

xs2 = [np.array(Nlst)]
ys2 = [np.array(k_lst)]
clp.plotone(xs2, ys2, axes[1], labels=["MMST"], colors = ['r-o'], xlabel="number of TLSs ($N$)", ylabel="Fitted decay rate ($k_{fit}$)", lw=2, markersize=8)


axes[1].legend(loc="lower right")

print(ys2)
# Plot x-axis as a function of PI
def format_func(value, tick_number):
    # find number of multiples of pi/2
    N = value
    if N == 0:
        return "0"
    else:
        return r"${0}\pi$".format(N)

axes[0].xaxis.set_major_formatter(FuncFormatter(format_func))

# output figure
clp.adjust(tight_layout=True, savefile="SE_Dicke.pdf")
