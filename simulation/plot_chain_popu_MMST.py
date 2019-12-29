import numpy as np
import columnplots as clp
import glob
import sys

path=sys.argv[-1]

kFGR = np.pi
N=101

from matplotlib.pyplot import FuncFormatter

def get_popu_data(filename):
    data = np.loadtxt(filename)
    t, p1, p2 = data[:-1,0], data[:-1, N-3], data[:-1, N+1]
    return t, p1, p2


def gather_data(path):
    t_QM, p1_QM, p2_QM = get_popu_data("%s/traj_QM.txt" %path)
    t_MMST, p1_MMST, p2_MMST = get_popu_data("%s/traj_MultiEh.txt" %path)
    xs = [t_MMST /kFGR, t_QM / kFGR, t_QM / kFGR]
    #y1s = [p1_MMST, p1_QM]
    y2s = [p2_MMST, p2_QM, p2_QM[0]*np.exp(-kFGR*t_QM)]
    return xs, y2s

labels = ["MMST", "QM", "free space decay"]
colors2 = ["r", "k", "0.7"]

axes = clp.initialize(col=1, row=3, width=4.3*2, height=4.3*0.618, fontsize=12, LaTeX=True, sharey=True, labelthem=True, labelthemPosition=[0.13, 0.95])

xs, y2s = gather_data("SE_101TLS_FDTD_dx_25")
clp.plotone(xs, y2s, axes[0], labels=labels, colors=colors2, xlabel="time ($t$)", ylabel="Excited state population", alphaspacing=0.1, lw=2)
axes[0].set_title("101 TLSs, $a = \lambda/2$")
xs, y2s = gather_data("SE_101TLS_FDTD_dx_13")
clp.plotone(xs, y2s, axes[1], labels=labels, colors=colors2, xlabel="time ($t$)",  alphaspacing=0.1, lw=2, showlegend=False)
axes[1].set_title("101 TLSs, $a = \lambda/4$")
xs, y2s = gather_data("SE_101TLS_FDTD_dx_0")
clp.plotone(xs[::-1], y2s[::-1], axes[2], labels=labels[::-1], colors=colors2[::-1], xlabel="time ($t$)",  alphaspacing=0.1, lw=2, showlegend=False)
axes[2].set_title("101 TLSs, $a = 0$")

# Plot x-axis as a function of PI
def format_func(value, tick_number):
    # find number of multiples of pi/2
    N = value
    if N == 0:
        return "0"
    else:
        return r"${0}\pi$".format(N)

for i in range(3):
    axes[i].xaxis.set_major_formatter(FuncFormatter(format_func))

# output figure
clp.adjust(tight_layout=True, savefile="chain_atoms.pdf")
