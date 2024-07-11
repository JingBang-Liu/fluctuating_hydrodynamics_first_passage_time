import pylab
import numpy as np
from matplotlib.gridspec import GridSpec
import seaborn as sns

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
#ratio = 2
ratio = 4./3                            # Sane ratio
#ratio = 2./(pylab.sqrt(5)-1.0)                # Aesthetic ratio
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
pylab.rcParams.update(params)

cm=sns.color_palette('Paired',8)
sns.set_palette(cm)

data_theory = np.loadtxt("data/SEG_epsis_tau_theory.txt")
data_MC = np.loadtxt("data/SEG_epsis_tau_MC.txt")
epsis = data_theory[:16,0]
epsis_MC = data_MC[:16,0]
naive = data_theory[:16,1]
mob = data_theory[:16,2]
con = data_theory[:16,3]
MC = data_MC[:,1]

pylab.semilogy(1/epsis, naive, '--', color=cm[0], label='original Eyring-Kramers')
pylab.semilogy(1/epsis, mob, '--', color=cm[1], label='with mobility correction')
pylab.semilogy(1/epsis, con, '-', color=cm[5], label='with conserved quantity')
pylab.semilogy(1/epsis_MC, MC, 'k.', label='Monte-Carlo simulation')
pylab.xlabel(r'noise strength $1/\varepsilon$')
pylab.ylabel(r'mean first passage time $w_B(x_-)$')
pylab.grid()
#sns.despine(left=True, bottom=True)

lgd=pylab.legend(loc='lower right', ncol=1, frameon=True)
lgd.get_frame().set_linewidth(0.0)

pylab.tight_layout(pad=0.25)
pylab.savefig("SEG-MC.pdf")
