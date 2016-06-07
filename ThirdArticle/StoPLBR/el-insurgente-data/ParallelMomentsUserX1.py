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

def PathSSLS(Uzero, seed):
	Up = StoPlbrmJC.SSLS(seed, Uzero)
	return Up

def LongTimePaths(seed, K):
	LongTimePath=StoPlbrmJC.LongTimeBehavior(seed, K)	
	return LongTimePath

def Figures():
	meanU = np.load("MeanShortTime.npy")
	meanSU = np.load("MeanSquareShortTime.npy")
	LongmeanU = np.load("MeanLongTime.npy")
	LongmeanSU = np.load("MeanSquareLongTime.npy")
	LongU1 = np.load("U1Paths.npy")
	LongU2 = np.load("U2Paths.npy")
	plt.close()
	w=1.0
	fig_width_pt =538.58*w                   # Get this from LaTeX using \show the \columnwidth(190mm)
	inches_per_pt = 1.0/72.27                # Convert pt to inch
	golden_mean = (np.sqrt(5) - 1.0) / 2.0   # A esthetic ratio
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
	
	fig1, ax1 = plt.subplots()
	ax1.plot(meanU[:,0], meanU[:,1],
		color='#000000',
		linestyle='-',
		lw=2,
		marker='',
		label='Mean'    
    )
	ax1.plot(np.sqrt(meanSU[:,0]),np.sqrt(meanSU[:,1]),
		color='#666666',
		linestyle='-',
		lw=2,
		marker='',
		label='Mean square'
	)
	ax1.set_xlabel(r'$u_1$')
	ax1.set_ylabel(r'$u_2$')
	ax1.legend(loc=0)
	plt.savefig("ShortTmeMoments.eps")
	
	fig2, ax2 = plt.subplots()
	ax2.plot(LongmeanU[:,0], LongmeanU[:,1],
		color='#000000',
		linestyle='-',
		lw=2,
		marker='',
		label='Mean'
    )
	ax2.plot(np.sqrt(LongmeanSU[:,0]), np.sqrt(LongmeanSU[:,1]),
		color='#666666',
		linestyle='-',
		lw=2,
		marker='',
		label='Mean square'
	)
	ax2.set_xlabel(r'$u_1$')
	ax2.set_ylabel(r'$u_2$')
	ax2.legend(loc=0)
	plt.savefig("LongTimeMoments.eps")
	
	fig3, ax3 = plt.subplots()
	for i in range(1):
		ax3.plot(LongU1[i,:], LongU2[i,:],
		#color='#000000',
		linestyle='-',
		lw=2,
		marker='',
		label='Mean'
		)
	ax3.set_xlabel(r'$u_1$')
	ax3.set_ylabel(r'$u_2$')
	plt.savefig("LongTimePaths.png")

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
k = 5
p = 0
r = p
T0 = 0.0
#T = 1
T = 2600
LTM=12
M=2**14
output = mp.Queue()

#
StoPlbrmJC = StoPLBRM()
StoPlbrmJC.InitializeMesh(k, p, r, T0, T)
StoPlbrmJC.NoiseUpdate(123456789)

#---------------------------------------------------------------------------------------------------------------------
StoPlbrmJC.SetParametersStoPLBRM(a1, b1, a2, b2, 1.0, gamma2, gamma1, 1.0, k1, k2, sigma, U0)
ussls = StoPlbrmJC.SSLS(123456789)
#
#
pool = mp.Pool(processes=8)
#UpathsSymb = [pool.apply_async(PathSSLS, args=(U0, i)) for i in np.random.random_integers(1, 123456789, M)]
#Upaths=[p.get() for p in UpathsSymb]
#Upaths = np.array(Upaths)
UpathsSymb = [pool.apply_async(LongTimePaths, args=(i,LTM)) for i in np.random.random_integers(1, 123456789, M)]
Upaths = [p.get() for p in UpathsSymb]
Upaths = np.array(Upaths)
np.save("U1Paths.npy", Upaths[:,:,0])
np.save("U2Paths.npy", Upaths[:,:,1])
del Upaths
#
U1Paths = np.load("U1Paths.npy")
meanU1 = U1Paths.mean(axis=0)
MSU1 = np.sqrt((U1Paths**2).mean(axis=0))
del U1Paths
#
U2Paths = np.load("U2Paths.npy")
meanU2 = U2Paths.mean(axis=0)
MSU2 = np.sqrt((U2Paths**2).mean(axis=0))
del U2Paths
#
meanU=np.array([meanU1,meanU2])
meanU=np.transpose(meanU)
np.save("MeanLongTime.npy", meanU)
del meanU
#
MSU = np.array([MSU1,MSU2])
MSU=np.transpose(MSU)
np.save("MeanSquareLongTime.npy", MSU)
del MSU
