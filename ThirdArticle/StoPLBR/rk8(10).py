import numpy as np
import tables
import matplotlib.pyplot as plt
import scipy as sci
from scipy.sparse import bsr_matrix
from scipy.integrate import odeint

def f(x, t):
    x1 = x[0]
    x2 = x[1]
    r1 = x1 * (2.0 - x2)
    r2 = x2 * (x1 - 1.0)
    return np.array([r1, r2])

def rk8_10(y, t, K):
    k_0 = f(y, t)
    k_list = [k_0]
    y_n =np.zeros([1, 2])
    for k in np.arange(1, K):
        sum = 0.0
        for j in np.arange(k):
            t_jj = t + h * alpha[j]
            sum += beta[k, j] * f(k_list[j], t_jj)
        k_j = k_0 + h * sum
        k_list.append(k_j)
    k_list = np.array(k_list)
    sum = 0.0
    for k in np.arange(len(k_list)):
        sum += c[k] * f(k_list[k], t)
    y_n = k_0 +h * sum
    return y_n

alpha = np.loadtxt('rk8(10)_alpha_coefficients.txt', dtype=np.float64,
                   delimiter=' ', usecols=(1,))
alpha = alpha.reshape([17, 1])
c = np.loadtxt('rk8(10)_c_coefficients.txt', dtype=np.float64,
               delimiter=' ', usecols=(1,))
c = c.reshape([17, 1])
beta_index_ij = np.loadtxt('rk8(10)_beta_coefficients.txt',
                           dtype=np.uint, delimiter=' ', usecols=(0, 1))
beta_data = np.loadtxt('rk8(10)_beta_coefficients.txt',
                       dtype=np.float64, delimiter=' ', usecols=(2,))
beta = bsr_matrix((beta_data, (beta_index_ij[:, 0],
                               beta_index_ij[:, 1])),
                  dtype=np.float64).toarray()
t_f = 4.0
t_j, h = np.linspace(0, t_f, 10000, dtype=np.float64, retstep=True)
y0 = np.array([2.0, 5])
sol = odeint(f, y0, t_j)
u_rk = np.zeros(sol.shape)
u_rk[0] = y0
buffer = np.array([t_j[0], y0[0],y0[1]])
# open an extensible file to add data
hdf5_path = "rk8(10)_sol.hdf5"
hdf5_file = tables.openFile(hdf5_path, mode='w')
filters = tables.Filters(complevel=5, complib='blosc')
data_storage = hdf5_file.createEArray(hdf5_file.root, 'data',
                                      tables.Atom.from_dtype(
                                          buffer.dtype),
                                      shape=(0, buffer.shape[-1]),
                                      filters=filters,
                                      expectedrows=len(buffer))
data_storage.append(buffer[:][None])
hdf5_file.close()
extendable_hdf5_file = tables.openFile(hdf5_path, mode='a')
extendable_hdf5_data = extendable_hdf5_file.root.data

for j in np.arange(u_rk.shape[0] - 1):
    u_rk[j+1] = rk8_10(u_rk[j], t_j[j], 17)
    buffer = np.array([t_j[j+1], u_rk[j+1, 0],  u_rk[j+1, 1]])
    extendable_hdf5_data.append(buffer[:][None])
extendable_hdf5_file.close()
plt.plot(sol[0], sol[1], label='lsoda')
plt.plot(u_rk[:, 0], u_rk[:, 1], label='rk8(10)')
plt.show()