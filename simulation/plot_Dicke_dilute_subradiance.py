import numpy as np
import columnplots as clp
import glob

kFGR = np.pi

from matplotlib.pyplot import FuncFormatter

def get_popu_data(filename, N):
    data = np.loadtxt(filename)
    t, p1, p2 = data[:-1,0], data[:-1, 1], data[:-1, 2]
    for i in range(N-1):
        p2 += data[:-1, 2*(i+2)]
    return t, p2 / N


def gather_data(Nlst, path):
    xs, ys = [], []
    for N in Nlst:
        subpath = path + "N_%d" %N
        t_MMST, p2_MMST = get_popu_data("%s/traj_MultiEh.txt" %subpath, N)
        xs.append(t_MMST)
        ys.append(p2_MMST)
    return xs, ys

def gather_data_same_N(N, path_lst):
    xs, ys = [], []
    for path in path_lst:
        subpath = path + "N_%d" %N
        t_MMST, p2_MMST = get_popu_data("%s/traj_MultiEh.txt" %subpath, N)
        xs.append(t_MMST)
        ys.append(p2_MMST)
    return xs, ys

# Plot x-axis as a function of PI
def format_func(value, tick_number):
    # find number of multiples of pi/2
    N = value
    if N == 0:
        return "0"
    else:
        return r"${0}\pi$".format(N)

axes = clp.initialize(col=2, row=1, width=4.3, height=4.3*0.618*2, fontsize=12, LaTeX=True, labelthem=True, labelthemPosition=[0.10, 0.95])

# Subplot 1: decay dynamics for different spacing with N = 35
Nlst = [1, 7, 35]
xs, ys = gather_data(Nlst, path="SE_Dicke_dx_13/")
labels = ["N=%d, $a = \lambda / 4$" %N for N in Nlst]
labels[0] = "N=1"
colors = ["r--", "c", "b"]
clp.plotone([x/kFGR for x in xs], ys, axes[0], labels=labels, colors=colors, xlabel="time ($t$)", ylabel="Avg. excited state popu. $\\bar{\\rho}_{ee}(t)$", lw=2, xlim=[-0.01, 1.5])#, ylog=True, ylim= [6e-3, 2.0])
axes[0].xaxis.set_major_formatter(FuncFormatter(format_func))

# Subplot 2: effective decay rate as a function of spacing for N = 35
Nlst = [7, 9, 11, 15, 19, 25, 35, 45]
xs, ys = gather_data(Nlst, path="SE_Dicke_dx_13/")

from scipy.optimize import curve_fit
def biexponenial_func(x, ks, kf, A):
    return A * np.exp(-ks*x) + (1.0 - A) * np.exp(-kf*x)

y_fit, ks_lst, kf_lst = [], [], []
for i in range(len(Nlst)):
    popt, pcov = curve_fit(biexponenial_func, xs[i][0:2500], ys[i][0:2500], p0=(0.1, 10, 0.5), maxfev=10000)
    y_fit.append(biexponenial_func(xs[i], *popt))
    ks = min(popt[0:2])
    kf = max(popt[0:2])
    ks_lst.append(ks)
    kf_lst.append(kf)

xs2 = [np.array(Nlst)]
ys2 = [1.0 / np.array(ks_lst)]
clp.plotone(xs2, ys2, axes[1], labels=["MMST"], colors = ['r-o'], xlabel="number of TLSs ($N$)", ylabel="Fitted slow decay lifetime ($1/k_{s}$)", lw=2, markersize=8)
axes[1].legend(loc="lower right")

axes[0].plot(xs[0] / kFGR, y_fit[0], "k-.", alpha=0.5, label="biexponential fitting")
axes[0].plot(xs[6] / kFGR, y_fit[6], "k-.", alpha=0.5)

print("kslow = ", ks_lst)
print("kfast = ", kf_lst)

# output figure
clp.adjust(tight_layout=True, savefile="SE_Dicke_dilute_subradiance.pdf")
