import numpy as np
import columnplots as clp
import glob
import sys

path="SE_Dicke/"#sys.argv[-1]

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

def gather_dpdt(Nlst, path = "SE_Dicke/"):
    x_avg, y_avg = gather_data(Nlst, path=path)
    x_new, y_new = calculate_derivative(x_avg, y_avg)
    return x_new, y_new

def calc_analytical(xs, ys, Nlst):
    xnew, ynew = xs, []
    for i in range(len(xs)):
        x, y = xs[i], ys[i]
        N = Nlst[i]
        td = x[np.argmax(y)]
        y_analytical = kFGR * N / 4.0 * np.cosh(kFGR*N / 2.0 * (x - td))**(-2)
        ynew.append(y_analytical)
    return xnew, ynew

def get_N_trajs(path="SE_Dicke/", Ntraj=10, Natoms=35):
    xs, ys = [], []
    for i in range(Ntraj):
        t, p2 = get_popu_data("%s/N_%d_statics/traj_MultiEh_%d.txt" %(path, Natoms, i+1), Natoms)
        nsize = np.size(p2)
        tnew = t[0:nsize-1]
        dpdt = (p2[0:nsize-1] - p2[1:]) / (t[1] - t[0])
        xs.append(tnew)
        ys.append(dpdt)
    return xs, ys

def get_td_lst(Nlst, path="SE_Dicke/"):
    xs, ys = gather_dpdt(Nlst, path=path)
    td_lst = []
    for i, N in enumerate(Nlst):
        t, dpdt = xs[i], ys[i]
        td = t[np.argmax(dpdt)]
        td_lst.append(td)
    return [np.array(Nlst)], [np.array(td_lst)]

def get_td_deviation_lst(Nlst, Ntraj=100, path="SE_Dicke/"):
    filename = path.split("/")[0] + "_td_statistics.dat"
    try:
        data = np.loadtxt(filename)
        return data[:,0], data[:,1], data[:,2], data[:,3], data[:,4]
    except:
        Delta_td_lst, Delta_pmax_lst = [], []
        td_mean_lst, pmax_mean_lst = [], []
        for N in Nlst:
            xs, ys = get_N_trajs(path=path, Ntraj=Ntraj, Natoms=N)
            td_lst, pmax_lst= [], []
            for i, x in enumerate(xs):
                td_lst.append(x[np.argmax(ys[i])])
                pmax_lst.append(np.max(ys[i]))
            td_lst = np.array(td_lst)
            pmax_lst = np.array(pmax_lst)
            Delta_td = np.sqrt(np.mean(td_lst**2) - np.mean(td_lst)**2)
            Delta_td_lst.append(Delta_td)
            td_mean_lst.append(np.mean(td_lst))
            Delta_pmax = np.sqrt(np.mean(pmax_lst**2) - np.mean(pmax_lst)**2)
            Delta_pmax_lst.append(Delta_pmax)
            pmax_mean_lst.append(np.mean(pmax_lst))
        data = np.zeros((len(Nlst), 5))
        data[:,0], data[:,1], data[:,2], data[:,3], data[:,4] = np.array(Nlst), np.array(td_mean_lst), np.array(Delta_td_lst), np.array(pmax_mean_lst), np.array(Delta_pmax_lst)
        np.savetxt(filename, data)
        return np.array(Nlst), np.array(td_mean_lst), np.array(Delta_td_lst), np.array(pmax_mean_lst), np.array(Delta_pmax_lst)

def get_td_deviation_QED(Nlst):
    Delta_td_lst = []
    td_mean_lst = []
    for N in Nlst:
        td, delta_td = 0.0, 0.0
        for i in range(N):
            td += 1.0 / N * (1.0 / (i+1))
            delta_td += (1.0 / N / (i+1))**2
        td_mean_lst.append(td / kFGR)
        Delta_td_lst.append(np.sqrt(delta_td) / kFGR)
    return np.array(Nlst), np.array(td_mean_lst), np.array(Delta_td_lst)

axes = clp.initialize(col=2, row=1, width=4.3, height=4.3*0.618*2, fontsize=12, LaTeX=True, labelthem=True, labelthemPosition=[0.04, 0.95])

# For the first plot, gather many trajs
xs, ys = get_N_trajs(path=path, Ntraj=5, Natoms=35)
x_avg, y_avg = gather_dpdt([35])

clp.plotone(xs, ys, axes[0],  labels=["No.%d MMST traj N=35" %(i+1) for i in range(len(xs))], xlabel="time ($t$)", ylabel="$d\\bar{\\rho}_{ee}/dt$", lw=1.5, xlim=[-0.005, 0.2], ylim=[-10, 105])
clp.plotone(x_avg, y_avg, axes[0], colors=["b"], labels=["MMST avg traj N=35"], xlabel="time ($t$)", ylabel="$d\\bar{\\rho}_{ee}/dt$", lw=3, xlim=[-0.005, 0.2], ylim=[-10, 105])

# For the second plot, we try to find the statics
Nlst = [7, 9, 11, 15, 19, 25, 35, 45]
x_td, y_td = get_td_lst(Nlst)
x_qm, y_qm, e_qm = get_td_deviation_QED(Nlst)
x, y, e, pmax, pmax_e = get_td_deviation_lst(Nlst=Nlst, Ntraj=1000, path=path)
axes[1].errorbar(x, y, yerr=e, linestyle='None', marker='o', fmt='r-o', capsize=5)
axes[1].errorbar(x_qm, y_qm, yerr=e_qm, linestyle='None', marker='^', fmt='k-^', alpha=0.9, capsize=5)
#clp.plotone([x, x_td[0], x_qm], [y, y_td[0], y_qm], axes[1],  labels=["statics of 1000 MMST trajectories", "MMST avg", "QM"], colors=['ro', "bs", "k^"], xlabel="number of TLSs ($N$)", ylabel="delay time ($t_D$)")
clp.plotone([x,  x_qm], [y,  y_qm], axes[1],  labels=["Histogram of 1000 MMST trajectories", "QM"], colors=['ro', "k^"], xlabel="number of TLSs ($N$)", ylabel="delay time ($t_D$)")

#axes[2].errorbar(x, pmax, yerr=pmax_e, linestyle='None', marker='o', fmt='r-o', capsize=5)
# plot the previous result
# output figure
clp.adjust(tight_layout=True, savefile="SE_Dicke_try.pdf")
