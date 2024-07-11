import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec
import seaborn as sns
import matplotlib.cm as cm2

def scientificTex(a, precision=2):
    expo=int(np.log(abs(a))/np.log(10));
    expo=expo if expo>0 else expo-1
    manti=np.round(float(a)/10**expo, precision)
    if expo<-1 or expo > 1:
        if manti == 1.0:
            return r'10^{{{}}}'.format(expo)
        else:
            return r'{}\cdot 10^{{{}}}'.format(manti,expo)
    else:
        return r'{}'.format(np.round(a,precision))

fig_width_pt = 320  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inch
# ratio = 2
ratio = 4./3                            # Sane ratio
# ratio = 2./(np.sqrt(5)-1.0)                # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width/ratio            # height in inches
fig_size =  [fig_width,fig_height]

sns.set_style("ticks")
params = {'backend': 'ps',
          'axes.labelsize': 10,
          'font.family': 'serif',
          'font.serif': 'Computer Modern Roman',
          'font.weight': 'normal',
          'legend.fontsize': 10,
          'xtick.labelsize': 10,
          'ytick.labelsize': 10,
          'text.usetex': True,
          'figure.figsize': fig_size}
plt.rcParams.update(params)

cm=sns.color_palette('Paired',8)
sns.set_palette(cm)

data = np.loadtxt('SEG_process.txt')
rp_process = data[1:,:]
times = np.squeeze(data[0,:])
Nx = 64
L = 1
dx = L/Nx
x = np.linspace(0,L-dx,Nx)
X, T = np.meshgrid(x,times)

fig, ax = plt.subplots()
extent = [times[0], times[-1], 0, L]
img = ax.imshow(rp_process, interpolation='bilinear', cmap=cm2.RdBu, 
               aspect= (extent[1]-extent[0])/(extent[3]-extent[2])/(4/3), extent=extent,
               origin='lower', vmax=rp_process.max(), vmin=rp_process.min())
plt.colorbar(img, ax=ax, label=r"$\rho$")
# ax.grid()
ax.set_xlabel(r"time")
ax.set_ylabel(r"$L$")
#sns.despine(left=True, bottom=True)

print(times[0])
print(times[-1])

# lgd=ax.legend(loc='lower right', ncol=1, frameon=True)
# lgd.get_frame().set_linewidth(0.0)

plt.tight_layout(pad=0.25)
plt.savefig('SEG_process.pdf')
