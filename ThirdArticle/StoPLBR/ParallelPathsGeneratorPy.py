# ===================================================================================
# This script generates the paths at Short and Long Time of the SDE Bone remodeling modeling
# under the assumption that the deterministic generation period is almost 650 days.
# INPUT:
#   StoPLBRM 	class which have the numerical schemes and other methods.
# OUTPUT:
#   U1PathsShortTime.npy
#   U2PathsShortTime.npy
#   U1PathsLongTime.npy
#   U2PathsLongTime.npy
#   Binary files with np-arrays of paths.
# =======================================================================================================================
import numpy as np
from StoPLBRJerezChen import StoPLBRM
import multiprocessing as mp
import time


def path_ssls(u_zero, seed):
    ppg.noise_update(seed)
    up = ppg.ssls(seed, u_zero)
    return up


def long_time_paths(seed, k):
    long_time_path = ppg.long_time_behavior(seed, k)
    return long_time_path


a1 = 0.15  # Estos son los parametros del modelo ai = \alpha_i
a2 = 0.16
b1 = 0.2
b2 = 0.02
ns = 0.1  # Amplitud del ruido
gamma1 = -0.837
gamma2 = 0.88

sigma = np.array([ns * b1, ns * b2])
k1 = 0.03
k2 = 0.0017
# Stencil Parameters
U0 = [10, 0.7]
k = 5
p = 0
r = p
T0 = 0.0
# T = 1
T = 650 * 8
LTM = 6
M = 2 ** 2

ppg = StoPLBRM()
ppg.initialize_mesh(k, p, r, T0, T)
ppg.set_parameters_sto_plbrm(a1, b1, a2, b2, 1.0, gamma2, gamma1, 1.0, k1, k2, sigma, U0)

output = mp.Queue()
pool = mp.Pool(processes=4)
# ================================================================================================
# Begin the Short Time Parallelization Block
# ================================================================================================
u_paths_symb = [pool.apply_async(path_ssls, args=(U0, i)) for i in np.random.random_integers(1, 123456789, M)]
u_paths = [p.get() for p in u_paths_symb]
u_paths = np.array([u_paths[i][0] for i in np.arange(M)])
t = time.time()
t = time.ctime(t)
t = t.replace(' ', '')
file_name1 = 'U1PathsShortTime' + t + '.npy'
file_name2 = 'U2PathsShortTime' + t + '.npy'
np.save(file_name1, u_paths[:, :, 0])
np.save(file_name2, u_paths[:, :, 1])
# del u_paths
# del u_paths_symb
'''
# =================================================================================================
# Begin the Long Time Parallelization Blockd
# =================================================================================================
u_paths_symb = [pool.apply_async(long_time_paths, args=(i, LTM)) for i in np.random.random_integers(1, 123456789, M)]
u_paths = [p.get() for p in u_paths_symb]
u_paths = np.array([u_paths[i][0] for i in np.arange(M)])
t = time.time()
t = time.ctime(t)
t = t.replace(' ', '')
file_name1 = 'U1PathsLongTime' + t + '.npy'
file_name2 = 'U2PathsLongTime' + t + '.npy'
np.save(file_name1, u_paths[:, :, 0])
np.save(file_name2, u_paths[:, :, 1])
'''