import numpy as np
import columnplots as clp
import glob
import sys

from matplotlib.pyplot import FuncFormatter

path = "GroundState_Stability/"

def get_popu_data(filename):
    data = np.loadtxt(filename)
    t, p1, p2 = data[:-1,0], data[:-1, 1], data[:-1, 2]
    return t, p1, p2

t_MMST, p1_MMST, p2_MMST = get_popu_data(path+"traj_MultiEh.txt")
t_p, p1_p, p2_p = get_popu_data(path+"traj_PreBinSQC.txt")
t_sed, p1_sed, p2_sed = get_popu_data(path+"traj_StochasticED.txt")


xs = [t_p /np.pi, t_sed/np.pi, t_MMST / np.pi]
y1s = [p1_p, p1_sed,  p1_MMST]
y2s = [p2_p, p2_sed,  p2_MMST]

labels = ["Sampling electronic ZPE (self-interaction)", "Sampling photonic ZPE (vacuum fluctuations)", "Sampling both (MMST)"]
colors1 = ["b--", "g-.", "r"]
colors2 = ["b--", "g--", "r--"]

ax = clp.initialize(col=1, row=1, width=4.3, fontsize=10, LaTeX=True)

clp.plotone(xs, y2s, ax, labels=labels, colors=colors1, xlabel="time ($t$)",
            ylabel="Excited state population", alphaspacing=0.0, lw=2, xlim=[0, 1.2], ylim=[-0.6, 1.0])

# Plot x-axis as a function of PI
def format_func(value, tick_number):
    # find number of multiples of pi/2
    N = int(value * 10)
    if N == 0:
        return "0"
    else:
        return r"${0}\pi$".format(N/10)

ax.xaxis.set_major_formatter(FuncFormatter(format_func))

# output figure
clp.adjust(tight_layout=True, savefile="ground_state_popu.pdf")
