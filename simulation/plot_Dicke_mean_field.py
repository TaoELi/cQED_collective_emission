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

def calculate_derivative(xs, ys):
    xnew, ynew = [], []
    for i in range(len(xs)):
        x, y = xs[i], ys[i]
        nsize = np.size(y)
        dydt = (y[0:nsize-1] - y[1:]) / (x[1] - x[0])
        x2 = x[0:nsize-1]
        xnew.append(x2)
        ynew.append(dydt)
    return xnew, ynew

def calc_analytical(xs, ys, Nlst):
    xnew, ynew = xs, []
    for i in range(len(xs)):
        x, y = xs[i], ys[i]
        N = Nlst[i]
        td = x[np.argmax(y)]
        y_analytical = kFGR * N / 4.0 * np.cosh(kFGR*N / 2.0 * (x - td))**(-2)
        ynew.append(y_analytical)
    return xnew, ynew

Nlst = [7, 35]
labels = ["MMST N = %d" %N for N in Nlst]

ax = clp.initialize(col=1, row=1, width=4.3, height=4.3*0.618, fontsize=12, LaTeX=True)

xs, ys = gather_data(Nlst)
xnew, dydt = calculate_derivative(xs,ys)

clp.plotone(xnew, dydt, ax, labels=labels, colors=["c", "b"], xlabel="time ($t$)", ylabel="$d\\bar{\\rho}_{ee}/dt$", lw=2, xlim=[-0.005, 0.2])

xnew, dydt = calc_analytical(xnew,dydt, Nlst)
labels = ["mean-field N = %d" %N for N in Nlst]
clp.plotone(xnew, dydt, ax, labels=labels, colors=["k-.", "k--"], xlabel="time ($t$)",  lw=2, xlim=[-0.005, 0.2], ylim=[-2, 30])



# output figure
clp.adjust(tight_layout=True, savefile="SE_Dicke_delay_time.pdf")
