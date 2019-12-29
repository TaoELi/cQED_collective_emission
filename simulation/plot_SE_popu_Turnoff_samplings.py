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

t_MMST, p1_MMST, p2_MMST = get_popu_data(path+"traj_MultiEh.txt")
t_QM, p1_QM, p2_QM = get_popu_data(path+"traj_QM.txt")

t_p, p1_p, p2_p = get_popu_data(path+"traj_PreBinSQC.txt")
t_sed, p1_sed, p2_sed = get_popu_data(path+"traj_StochasticED.txt")


xs = [t_p /np.pi, t_sed/np.pi, t_MMST / np.pi, t_QM / np.pi]
y1s = [p1_p, p1_sed,  p1_MMST, p1_QM]
y2s = [p2_p, p2_sed,  p2_MMST, p2_QM]


labels = ["Sampling electronic ZPE (self-interaction)", "Sampling photonic ZPE (vacuum fluctuations)", "Sampling both (MMST)", "QM"]
colors1 = ["b--", "g--", "r--", 'k--']
colors2 = ["b--", "g-.", "r", 'k--']

ax = clp.initialize(col=1, row=1, width=4.6, fontsize=12, LaTeX=True, sharex=True)

clp.plotone(xs, y2s, ax, labels=labels, colors=colors2, xlabel="time (t)",
           ylabel="Excited state population", alphaspacing=0.0, lw=2, ylim=[-0.51, 2.0])


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
clp.adjust(tight_layout=True, savefile="SE_popu_turnoff_samplings.pdf")
