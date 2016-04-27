import numpy as np
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
from StoPLBRJerezChen import StoPLBRM
from matplotlib.pyplot import flag
from matplotlib import colors
from matplotlib import rcParams
from mpl_toolkits.mplot3d import Axes3D
# Model parameters -------------------------------------------------------------------------------
a1 = 0.3			#Estos son los parametros del modelo ai = \alpha_i
a2 = 0.18
b1 = 0.2
b2 = 0.02
ns = 0.1			#Amplitud del ruido
gamma1 = -0.9
gamma2 = 0.5
sigma = np.array([ns*b1, ns*b2])
k1 = 0.03
k2 = 0.0017
#
Ji1 = (( b1 + 0.5 * (sigma[0] ** 2) ) / a1) ** (1.0 / gamma1)
Ji2 = (( b2 + 0.5 * (sigma[1] ** 2) ) / a2) ** (1.0 / gamma2)
#Stencil Parameters
U0=[10, 0.700]
k = 5
p = 0
r = p
T0 = 0.0
T = 650*5
LTM=8
M=500
#
seed = 109730187
np.random.seed(seed)
#
StoPlbrmJC = StoPLBRM()
StoPlbrmJC.InitializeMesh(k, p, r, T0, T)
#------------------------------------------------------------------------------
StoPlbrmJC.SetParametersStoPLBRM(a1, b1, a2, b2, 1.0, gamma2, gamma1, 1.0, k1, k2, sigma, U0)
Urk = StoPlbrmJC.RK()
StoPlbrmJC.NoiseUpdate(seed,flag=1)
Ussls = StoPlbrmJC.SSLS(seed, [1,1], fn = 1.0)

Uas, detUas = StoPlbrmJC.LongTimeBehavior(seed, LTM)
StoPlbrmJC.SaveData()
t=StoPlbrmJC.t[0 : -1 : StoPlbrmJC.R].reshape([StoPlbrmJC.t[0 : - 1 : StoPlbrmJC.R].shape[0],1])
tas=StoPlbrmJC.LongTimeM * t[-1] + t
U1rk = Urk[:, 0]
U2rk = Urk[:, 1]
#
stoU1 = Ussls[:, 0]
stoU2 = Ussls[:, 1]
#
detU1as = detUas[:, 0]
detU2as = detUas[:, 1]
stoU1as = Uas[:, 0]
stoU2as = Uas[:, 1]
'''
np.savetxt('SolutionDet.txt',\
		np.transpose(np.array([t[:,0], U1rk[:] , U2rk[:]])),\
		fmt=['%1.8f','%1.8f','%1.8f'],\
		delimiter='\t'\
)
#
np.savetxt('OneLongPathSolutionSto.txt',\
		np.transpose(np.array([t[:,0], stoU1[:] , stoU2[:]])),\
		fmt=['%1.8f','%1.8f','%1.8f'],\
		delimiter='\t'\
)
#
np.savetxt('SolutionDetAs.txt',\
		np.transpose(np.array([tas[:,0], detU1as[:] , detU2as[:]])),\
		fmt=['%1.8f','%1.8f','%1.8f'],\
		delimiter='\t'\
)
#
np.savetxt('SolutionSto.txt',\
		np.transpose(np.array([t[:,0], stoU1[:] , stoU2[:]])),\
		fmt=['%1.8f','%1.8f','%1.8f'],\
		delimiter='\t'\
)
#
np.savetxt('SolutionStoAs.txt',\
		np.transpose(np.array([tas[:,0], stoU1as[:] , stoU2as[:]])),\
		fmt=['%1.8f','%1.8f','%1.8f'],\
		delimiter='\t'\
)
'''
np.save('SolutionDet.npy', np.transpose(np.array([t[:,0], U1rk[:] , U2rk[:]])))
np.save('SolutionDetAs.npy', np.transpose(np.array([tas[:,0], detU1as[:] , detU2as[:]])))
np.save('SolutionSto.npy', np.transpose(np.array([t[:,0], stoU1[:] , stoU2[:]])))
np.save('SolutionStoAs.npy', np.transpose(np.array([tas[:,0], stoU1as[:] , stoU2as[:]])))
#----------------------Long Path-----------------------------------------------------------------------

plt.figure()
plt.plot(t, stoU2)
plt.figure()
plt.plot(tas, stoU2as)
plt.show()

#
T = 50 * 650
k = 6
#seed = np.random.random_integers(1, 123456789)
seed = 109730187
np.random.seed(seed)
StoPlbrmJC.InitializeMesh(k, p, r, T0, T)
StoPlbrmJC.SetParametersStoPLBRM(a1, b1, a2, b2, 1.0, gamma2, gamma1, 1.0, k1, k2, sigma, U0)
t = StoPlbrmJC.t[0 : -1 : StoPlbrmJC.R].reshape([StoPlbrmJC.t[0 : - 1 : StoPlbrmJC.R].shape[0],1])
StoPlbrmJC.NoiseUpdate(seed, flag=1)
Urk = StoPlbrmJC.RK()
Ussls = StoPlbrmJC.SSLS(seed, U0, fn = 1.0)
U1 = Urk[:, 0]
U2 = Urk[:, 1]
stoU1 = Ussls[:, 0]
stoU2 = Ussls[:, 1]
#
plt.figure()
plt.plot(t, stoU1)
plt.figure()
plt.plot(t, stoU2)
plt.show()
np.save('OneLongPathSolutionSto.npy', np.transpose(np.array([t[:,0], stoU1[:] , stoU2[:]])))
np.save('OneLongPathSolutionDet.npy', np.transpose(np.array([t[:,0], U1[:] , U2[:]])))