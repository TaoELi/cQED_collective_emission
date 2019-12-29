import numpy as np
import columnplots as clp
import glob
import sys

from matplotlib.pyplot import FuncFormatter

path="SE_1TLS_FDTD/"

def get_popu_data(filename):
    data = np.loadtxt(filename)
    t, p1, p2 = data[:-1,0], data[:-1, 1], data[:-1, 2]
    return t, p1, p2

t_QM, p1_QM, p2_QM = get_popu_data(path+"traj_QM.txt")

t_MMST, p1_MMST, p2_MMST = get_popu_data(path+"traj_MultiEh.txt")

xs = [t_MMST /np.pi, t_QM / np.pi]
y1s = [p1_MMST, p1_QM]
y2s = [p2_MMST, p2_QM]

labels = ["MMST", "QM"]
colors1 = ["r--", "k--"]
colors2 = ["r", "k"]

ax = clp.initialize(col=1, row=1, width=4.3, fontsize=12, LaTeX=True)

clp.plotone(xs, y1s, ax, colors=colors1, alphaspacing=0.0, lw=2)
clp.plotone(xs, y2s, ax, labels=labels, colors=colors2, xlabel="time ($t$)",
            ylabel="Electronic population", alphaspacing=0.0, lw=2)

# Plot x-axis as a function of PI
def format_func(value, tick_number):
    # find number of multiples of pi/2
    N = value
    if N == 0:
        return "0"
    else:
        return r"${0}\pi$".format(N)

ax.xaxis.set_major_formatter(FuncFormatter(format_func))

# output figure
clp.adjust(tight_layout=True, savefile="SE_popu_MMST.pdf")
