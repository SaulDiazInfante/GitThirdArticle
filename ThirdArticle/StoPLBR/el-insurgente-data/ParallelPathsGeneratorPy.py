#=======================================================================================================================
# This script generates the paths at Short and Long Time of the SDE Bone remodeling modeling
# under the assumption that the deterministic generation period is almost 650 days.
#INPUT:
#	StoPLBRM 	class which have the numerical schemes and other methods.
#OUTPUT:
#	U1PathsShortTime.npy
#	U2PathsShortTime.npy
#	U1PathsLongTime.npy
#	U2PathsLongTime.npy
#	Binary files with np-arrays of paths. 
#=======================================================================================================================
import numpy as np
import matplotlib.pyplot as plt
from StoPLBRJerezChen import StoPLBRM
import multiprocessing as mp
import time
def PathSSLS(Uzero, seed):
	Up = StoPlbrmJC.SSLS(seed, Uzero)
	return Up

def LongTimePaths(seed, K):
	LongTimePath=StoPlbrmJC.LongTimeBehavior(seed, K)	
	return LongTimePath
#
a1 = 0.15
a2 = 0.16
b1 = 0.2
b2 = 0.02
gamma1 = -0.837
gamma2 = 0.88
sigma = np.array([.1*b1, 0.1*b2])
k1 = 0.03
k2 = 0.0017

#Stencil Parameters
U0=[10, 0.7]
k = 5
p = 0
r = p
T0 = 0.0
#T = 1
T = 650*8
LTM = 5
M = 2 ** 12

StoPlbrmJC = StoPLBRM()
StoPlbrmJC.InitializeMesh(k, p, r, T0, T)
StoPlbrmJC.SetParametersStoPLBRM(a1, b1, a2, b2, 1.0, gamma2, gamma1, 1.0, k1, k2, sigma, U0)
output = mp.Queue()
pool = mp.Pool(processes=4)
#=======================================================================================================================
# Begin the Short Time Parallelization Block
#=======================================================================================================================

t = time.time()
print t
UpathsSymb = [pool.apply_async(PathSSLS, args=(U0, i)) for i in np.random.random_integers(1, 123456789, M)]
Upaths=[p.get() for p in UpathsSymb]
Upaths = np.array(Upaths)
t = time.time()
t = time.ctime(t)
t = t.replace(' ', '')
file_name1 = 'U1PathsShortTime' + t + '.npy'
file_name2 = 'U2PathsShortTime' + t + '.npy'
np.save(file_name1, Upaths[:, :, 0])
np.save(file_name2, Upaths[:, :, 1])

del Upaths
del UpathsSymb
#=======================================================================================================================
# Begin the Long Time Parallelization Block
#=======================================================================================================================
UpathsSymb = [pool.apply_async(LongTimePaths, args=(i,LTM)) for i in np.random.random_integers(1, 123456789, M)]
Upaths = [p.get() for p in UpathsSymb]
Upaths = np.array(Upaths)
t = time.time()
t = time.ctime(t)
t = t.replace(' ', '')
file_name1 = 'U1PathsLongTime' + t + '.npy'
file_name2 = 'U2PathsLongTime' + t + '.npy'
np.save(file_name1, Upaths[:, :, 0])
np.save(file_name2, Upaths[:, :, 1])
