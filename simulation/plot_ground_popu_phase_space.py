import numpy as np
import glob
import matplotlib.pyplot as plt
from scipy import stats
import columnplots as clp

path = "GroundState_Stability/"

n1 = 0.5
n2 = 0.5
cutoff = 2.0

xmin = n1-cutoff
xmax = n1+cutoff
ymin = n2-cutoff
ymax = n2+cutoff

methods = ["Sampling both \n (MMST)", "Sampling electronic ZPE \n (self-interaction)", "Sampling photonic ZPE \n (vacuum fluctuations)"]
# begin to plot the contour plot
def density_estimation(m1, m2):
    X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    positions = np.vstack([X.ravel(), Y.ravel()])
    values = np.vstack([m1, m2])
    kernel = stats.gaussian_kde(values)
    Z = np.reshape(kernel(positions).T, X.shape)
    return X, Y, Z

n = 3

NSTART = 0

fig, axes = clp.initialize(col=3, row=n, width=6, height=6, sharex=True, sharey=True,
                    commonX=[0.5, 0.00, "ground state popu."], commonY=[0.01, 0.5, "excited state popu."],
                    return_fig_args=True, LaTeX=True, labelthem=True, labelthemPosition=[0.95, 0.95])
Zgood = 0

for k in range(3):
    if k == 0:
        data_tot = np.loadtxt(path+"ElectronicDist_MultiEh.txt")
    elif k == 1:
        data_tot = np.loadtxt(path+"ElectronicDist_PreBinSQC.txt")
    else:
        data_tot = np.loadtxt(path+"ElectronicDist_StochasticED.txt")
    for i in range(n):
        m1, m2  = data_tot[:, NSTART+i*2], data_tot[:, NSTART+i*2+1]
        if i is 2:
            m1, m2  = data_tot[:, NSTART+3*2], data_tot[:, NSTART+3*2+1]
        axes[k, i].plot(m1, m2, 'k.', markersize=0.5, alpha=0.5)
        try:
            if k is 2 and i is 0:
                im = axes[k, i].imshow(np.rot90(Z*0.0), cmap=plt.cm.gist_earth_r, extent=[xmin, xmax, ymin, ymax])
            else:
                X, Y, Z = density_estimation(m1, m2)
                im = axes[k, i].imshow(np.rot90(Z), cmap=plt.cm.gist_earth_r, extent=[xmin, xmax, ymin, ymax])
                if k is 0 and i is 0:
                    Zgood = Z
                    im_good = im
        except:
            im = axes[k, i].imshow(np.rot90(Zgood*0.0), cmap=plt.cm.gist_earth_r, extent=[xmin, xmax, ymin, ymax])
            print("Error for inv")
        m1_centroid = np.mean(m1)
        m2_centroid = np.mean(m2)
        print(m1_centroid, m2_centroid)
        axes[k, i].axvline(x=1.0, color='0.5', linestyle='--', alpha=0.5)
        axes[k, i].axhline(y=0.0, color='0.5', linestyle='--', alpha=0.5)
        axes[k, i].plot(m1_centroid, m2_centroid, 'ro', markersize=5, alpha=0.8)
        axes[k, i].set_xlim([xmin, xmax])
        axes[k, i].set_ylim([ymin, ymax])
        if i == 0:
            axes[k, i].text(xmin+0.2, ymax-0.5, "t = 0")
            axes[k, i].text(xmin+0.2, ymin+0.2, methods[k], color='b', fontsize=10)
        elif i == 1:
            axes[k, i].text(xmin+0.2, ymax-0.5, "t = $\pi$/3")
        elif i == 2:
            axes[k, i].text(xmin+0.2, ymax-0.5, "t = $\pi$")
        else:
            axes[k, i].text(xmin+0.2, ymax-0.5, "t = %d$\pi$/3" %i)

cb_ax = fig.add_axes([1.0, 0.1, 0.02, 0.8])
cbar = fig.colorbar(im_good, cax=cb_ax)

clp.adjust(tight_layout=True, savefile="ground_state_phase_space_dist.pdf")
