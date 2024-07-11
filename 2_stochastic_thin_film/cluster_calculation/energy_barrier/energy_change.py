import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec
import seaborn as sns
import matplotlib.cm as cm2
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

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

fig_width_pt = 446  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inch
# ratio = 2
# ratio = 4./3                            # Sane ratio
ratio = 2./(np.sqrt(5)-1.0)                # Aesthetic ratio
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

# for energy change
data = np.loadtxt('STFE_U_all_0.0005.txt') # energy data
U_s = data[0]
U_s_sim = data[1]
U_a = data[2]
ind_loss = int(data[3])
U_all = data[4:]
dt = 1.5665e-7
Nx = 128
L = 1
dx = L/Nx
x = np.linspace(0,L-dx,Nx)
Nt = int(5e5)
t = np.linspace(0,Nt-1,Nt)*dt
t = t-t[-1]
end_cut = 10
epsilon = 0.0005 # noise amplitude
# Ua = -1.290017893812943e4 # this is energy for flat profile normalised by epsilon
# Us = -1.289524373725464e4 # this is energy for saddle normalised by epsilon
del(data)

# for profile
data_profiles = np.loadtxt('STFE_profiles_0.0005.txt') # profile data
x = data_profiles[:,0]
profiles = data_profiles[:,1:]
del(data_profiles)


fig, ax = plt.subplots()
ax.plot(t[:Nt-end_cut],U_all[:Nt-end_cut]-(U_s_sim-U_s)-U_a,color=cm[1],label=r"$E$",linewidth=0.3)
ax.plot(t[:Nt-end_cut],np.ones(Nt-end_cut)*U_s-U_a,color=cm[2],label=r"$E[h_s]$",linewidth=1.0)
ax.plot(t[:Nt-end_cut],np.ones(Nt-end_cut)*U_a-U_a,color=cm[4],label=r"$E[h_0]$",linewidth=1.0)
# ax.plot(np.ones(10)*t[ind_loss],np.linspace(-2,6,10),"--",color="k",label=r"t_{saddle}",linewidth=0.4)
# ax.grid()
print(t[ind_loss])
ax.set_xlabel(r"$t-t_r$")
ax.set_ylabel(r"$(E-E[h_0])/ \varepsilon$")
#sns.despine(left=True, bottom=True)
ax.set_ylim([U_a/10*epsilon,U_s-U_a/10*epsilon-U_a])
ax.set_xlim([t[0],t[Nt-end_cut]])

lgd=ax.legend(bbox_to_anchor=(0.82,0.52), ncol=1)
# lgd.get_frame().set_linewidth(0.0)

# inset plot for profile
# positioning of inset plot, starting from (inset_dim[0],inset_dim[1]), ending at (inset_dim[2],inset_dim[3])
inset_dim = [0.1,0.4,0.5,0.45]
# inset_dim[3] = inset_dim[2]/ratio
print(inset_dim)
axins = inset_axes(ax, width="95%", height="95%",
                   bbox_to_anchor=(inset_dim[0],inset_dim[1],inset_dim[2],inset_dim[3]),
                   bbox_transform = ax.transAxes, loc=3)
# ax.add_patch(plt.Rectangle((inset_dim[0], inset_dim[1]), inset_dim[2], inset_dim[3], ls="--", ec="c", fc="none",
#                            transform=ax.transAxes)) # show inset box
axins.plot(x,profiles[:,0],color=cm[0],linewidth=1.0)
axins.plot(x,profiles[:,1],color=cm[0],linewidth=1.0)
axins.plot(x,profiles[:,2],color=cm[0],linewidth=1.0)
axins.plot(x,profiles[:,3],"--",color="k",label=r"$h_s$",linewidth=0.9)
axins.plot(x,profiles[:,4],color=cm[0],linewidth=1.0)
axins.set_xlabel(r"$x$")
axins.set_ylabel(r"$h$")
axins.set_xlim([0,1])

# plot arrows
tempxA = t[10]; tempxB = 0.1
tempyA = 0; tempyB = np.mean(profiles[:,0])
axins.annotate('', xy=(tempxA, tempyA),xycoords=ax.transData,
            xytext=(tempxB, tempyB),textcoords=axins.transData,
            arrowprops=dict(shrinkA=0,shrinkB=0,color='black',arrowstyle="<->",linewidth=1)
            )

tempxA = t[int(3e5)]; tempxB = x[64]
tempyA = U_all[int(3e5)]-(U_s_sim-U_s)-U_a; tempyB = profiles[64,1]
axins.annotate('', xy=(tempxA, tempyA),xycoords=ax.transData,
            xytext=(tempxB, tempyB),textcoords=axins.transData,
            arrowprops=dict(shrinkA=0,shrinkB=0,color='black',arrowstyle="<->",linewidth=1)
            )

tempxA = t[ind_loss]; tempxB = x[120]
tempyA = U_all[ind_loss]-(U_s_sim-U_s)-U_a; tempyB = profiles[120,2]
axins.annotate('', xy=(tempxA, tempyA),xycoords=ax.transData,
            xytext=(tempxB, tempyB),textcoords=axins.transData,
            arrowprops=dict(shrinkA=0,shrinkB=0,color='black',arrowstyle="<->",linewidth=1)
            )

tempxA = t[ind_loss+35000]; tempxB = x[86]
tempyA = U_all[ind_loss+35000]-(U_s_sim-U_s)-U_a; tempyB = profiles[86,4]
axins.annotate('', xy=(tempxA, tempyA),xycoords=ax.transData,
            xytext=(tempxB, tempyB),textcoords=axins.transData,
            arrowprops=dict(shrinkA=0,shrinkB=0,color='black',arrowstyle="<->",linewidth=1)
            )


plt.tight_layout(pad=0.25)
# plt.show()
plt.savefig('STFE_energy.pdf')
