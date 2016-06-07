import numpy as np
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
from StoPLBRJerezChen import StoPLBRM
from matplotlib.pyplot import flag
from matplotlib import colors
from matplotlib import rcParams
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import multiprocessing as mp

def CreateFigures(LTM=10):
    urk = np.loadtxt('SolutionDet.txt')
    ussls = np.loadtxt('SolutionSto.txt')
    uas = np.loadtxt('SolutionStoAs.txt')
    t= ussls[:,0]
    #Pots Vs Time
    plt.close()
    w=1.0
    fig_width_pt =538.58*w                   # Get this from LaTeX using \showthe\columnwidth(190mm)
    inches_per_pt = 1.0/72.27                # Convert pt to inch
    golden_mean = (np.sqrt(5) - 1.0) / 2.0   # Aesthetic ratio
    fig_width = fig_width_pt * inches_per_pt # width in inches
    fig_height = fig_width*golden_mean       # height in inches
    fig_size =  [fig_width, fig_height]
    params = {'backend': 'ps',
        'axes.labelsize': 10,
        'text.fontsize': 10,
        'legend.fontsize': 8,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'text.usetex': True,
        'figure.figsize': fig_size}
    rcParams.update(params)
    
##
    fig1, ax1 = plt.subplots()
    #plt.title('Jerez-Chen Model Experimental Parameters')
    OC_line = mlines.Line2D([], [], color='#000000', linestyle='-', lw=2, marker='',  markersize=1)
    OB_line = mlines.Line2D([], [], color='#666666', linestyle='-', lw=2, marker='',  markersize=1)
    sOC_line = mlines.Line2D([], [], color='#666666', linestyle='-.', lw=4, marker='', markersize=1)
    sOB_line = mlines.Line2D([], [], color='#000000', linestyle='-.', lw=4, marker='', markersize=1)
#
    ax1.plot(urk[:, 0], urk[:, 1],
             color='#000000',
             linestyle='-',
             lw=2,
             marker='',
             label='osteoclast'
    )
    ax1.plot(ussls[:, 0], ussls[:, 1],
             color='#666666', 
             linestyle='-.', 
             lw=3,
             marker='',
             label='Sto osteoclast'
    )
    ax1.set_xlabel('time (days)')
    ax1.set_ylabel('OCs (u1)', color='k')
    ax2 = ax1.twinx()
    #
    ax2.plot(urk[:, 0], urk[:, 2],
             color='#666666',
             lw = 2,
             linestyle='-',
             marker='',
             label='osteoblast'
    )
    ax2.plot(ussls[:, 0], ussls[:, 2],
             color='#000000',
             linestyle='-.',
             lw=3,
             marker='',
             label='Sto osteoblast'
    )
    ax2.set_xlabel('time (days)')
    ax2.set_ylabel(r'OBs (u2)')
    ax2.grid(True)
    plt.legend([OC_line, OB_line, sOC_line, sOB_line],
               ["OCs", "OBs", "stoOCs", "stoOBs"],
               bbox_to_anchor=(0., 1.02, 1., .5), 
               loc=3,
               ncol=4, 
               mode="expand",
               borderaxespad=0.
    )
    plt.savefig("CellsTime.eps")
    plt._show()
    #   
    #
    #Phase Potrait
    fig1, ax3 = plt.subplots()
    #plt.title('Stochastic Power Law Bone Remodeling')
    ax3.set_xlim(-1, 13)
    ax3.set_ylim(-50, urk[:,2].max()*1.15)
    ax3.plot(urk[:, 1], urk[:, 2],
             color='#000000',
             linestyle='-',
             lw=3,
             label='Deterministic'
    )
    ax3.plot(StoPlbrmJC.Ubar[0], StoPlbrmJC.Ubar[1],
             'o',
             mfc='none',
             ms=8,
             #label=r'$\bar{u}$'
    )
    ax3.plot(StoPlbrmJC.Uzero[0], StoPlbrmJC.Uzero[1],
             's', 
             mfc='#666666',
             ms=8
    )
    ax3.plot(ussls[:, 1], ussls[:, 2],
             color='#666666', 
             ms=1,
             label='Stochastic Short Time'
    )
    ax3.plot(uas[:, 1], uas[:, 2],
             color='#666666',
             lw=1,
             linestyle='-.',
             label='Long Time Stochastic')
    #ax3.plot(uas[0,0], uas[0,1], 'x', mfc='green', ms=12)
    ax3.set_xlabel(r'$u_1$')
    ax3.set_ylabel(r'$u_2$')
    ax3.grid(True)
    ax3.legend(loc=0)
    plt.savefig("PhasePotrait.eps")
    plt.show()
   #
   #
    fig = plt.figure()
    #ax = fig.gca(projection='3d')
    ax = Axes3D(fig)
    ax.set_xlabel(r'$u_1$')
    ax.set_ylabel(r'$t$')
    ax.set_zlabel(r'$u_2$')
    ax.view_init(elev=20., azim=195)
    #
    ax.plot(urk[:, 1], t, urk[:, 2],
            color='#000000',
            lw=1,
            linestyle='-',
            label='Deterministic'
    )
    ax.plot(ussls[:, 1], t, ussls[:, 2],
            color='#666666',
            lw=1,
            linestyle='-',
            #label='Sto Short Time'
    )
    ax.plot(urk[:, 1], LTM*t[-1]+t, urk[:, 2],
            color='#000000',
            lw=1,
            linestyle='-',
            #label='Det Long Time'
    )
    ax.plot(uas[:, 1], LTM*t[-1]+t, uas[:, 2],
            color='#666666',
            lw=1,
            linestyle='-',
            label='Stochastic'
    )
    ax.legend()
    plt.savefig("PhasePotrait3d.eps")
    plt.show()

# Model parameters -------------------------------------------------------------------------------
a1 = 0.3
a2 = 0.1
b1 = 0.2
b2 = 0.02
gamma1 = -0.3
gamma2 = 0.5
sigma = np.array([.1*b1, 0.1*b2])
k1 = 0.03
k2 = 0.0017

#Stencil Parameters
U0=[11, 237]
k = 7
p = 1
r = p
T0 = 0.0
#T = 1300
T = 2600
LTM=1
M=100
#
StoPlbrmJC = StoPLBRM()
StoPlbrmJC.InitializeMesh(k, p, r, T0, T)
#StoPlbrmJC.NoiseUpdate()
#------------------------------------------------------------------------------
StoPlbrmJC.SetParametersStoPLBRM(a1, b1, a2, b2, 1.0, gamma2, gamma1, 1.0, k1, k2, sigma, U0)
#StoPlbrmJC.RK()
#StoPlbrmJC.SSLS()
#StoPlbrmJC.LongTimeBehavior(LTM)
#StoPlbrmJC.Moments(M, LTM)
#StoPlbrmJC.SaveData()
#CreateFigures(LTM)
