import numpy as np
import matplotlib.pyplot as plt
import time

# Deterministic model function


def F(u):
    x, y, z, w = u  # unpack current values of x
    derivs = np.array([a1 * x * y**g1 - b1 * x + s1 * x * z,
                       a2 * y * x**g2 - b2 * y + s2 * y * z, r * z *
                       (1 - z / K) - b3 * z + s3 * z * x ** g2 +
                       s4 * z * y ** g1, - k1 * max(x - y0, 0) + k2 *
                       max(y - x0, 0)])
    # list of dx/dt=f functions
    return derivs

# Stochastic model's functions
# Drift Term


def f(x, h):
    u, v, w = x
    increment = np.array([a1 * v ** g1 - b1 + s1 * w,
                          a2 * u ** g2 - b2 + s2 * w,
                          r * (1 - w / K) - b3 + s3 * u ** g2 + s4 *
                          v ** g1])
    derivs = np.exp(h * increment)
    return derivs
# Diffusion Term


def g(x):
    u, v, w = x
    derivs = np.array([d1 * u, d2 * v, d3 * w * (1 - w / K)])
    return derivs


# Bone Mass
def G(u, v):
    return -k1 * np.sqrt(max(u - i1, 0)) + k2 * np.sqrt(max(v - i2, 0))

# np.random.seed(10)
# Parameters

a1 = 0.3
a2 = 0.1
b1 = 0.2
b2 = 0.02
b3 = 0.05
g1 = -0.3
g2 = 0.5
s1 = 0.001
s2 = -0.00003
s3 = 0.005
s4 = 0.005
r = 0.045
K = 300.0
d1 = b1 * 0.1
d2 = b2 * 0.1
d3 = r * 0.1

k1 = 0.075
k2 = 0.00152

x0 = (b2 / a2) ** (1 / g2)
y0 = (b1 / a1) ** (1 / g1)

#
i1 = ((b2 + 0.5 * d2 ** 2)/a2) ** (1 / g2)
i2 = ((b1 + 0.5 * d1 ** 2)/a1) ** (1 / g1)

# cancer

l = r - b3 + s3*b2/a2 + s4 * b1 / a1
m = r / K + s2 * s3/a2 + s1 * s4/a1
j3 = ((K ** 2)/(d3 ** 2)) * (-m + (d3 ** 2) / K
      + np.sqrt(m ** 2 + (2 * m * (d3 ** 2) / (K ** 2)) * (l / m - K)
      + ((d3 ** 2) / K ** 2) * ((s3 * d2 ** 2) / a2 + (s4 * d1 **
                                                      2)/a1)))
j2 = ((b2 + 0.5 * d2 ** 2 - s2 * j3) / a2) ** (1/g2)
j1 = ((b1 + 0.5 * d1 ** 2 - s1 * j3) / a2) ** (1/g1)

T = 15000
N = 2 ** 21
dt = np.float(T) / np.float(N)
# L = T * N
dW = np.sqrt(dt) * np.random.randn(N, 3)
h = 0.02
n = int(T / h)
# Step size and nodes
t1 = np.linspace(0.0, T, n+1)
# Deterministic Model
det = np.zeros((n + 1, 4))
det[0] = [10.0, 5.0, 50.0, 40.0]
for j in np.arange(n):
    q1 = F(det[j]) * h
    q2 = F(det[j] + 0.5 * q1) * h
    q3 = F(det[j] + 0.5 * q2) * h
    q4 = F(det[j] + q3) * h
    det[j+1] = det[j] + 1.0 / 6.0 * (q1 + 2 * q2 + 2 * q3 + q4)

# Stochastic Model
Xsm = np.zeros((N, 3))
mo = np.zeros(N)
x_temp1 = [10.0, 5.0, 50.0]
Xsm[0, :] = x_temp1
mo[0] = 95.0
for j in np.arange(N - 1):
    # Steklov Method
    x_temp1 = f(Xsm[j, :], dt) * Xsm[j, :] + g(Xsm[j, :]) * dW[j]
    Xsm[j+1, :] = x_temp1
    K1 = G(Xsm[j][0], Xsm[j][1])*dt
    K2 = G(Xsm[j][0] + 0.5 * K1, Xsm[j][1] + 0.5 * K1) * dt
    K3 = G(Xsm[j][0] + 0.5 * K2, Xsm[j][1] + 0.5 * K2) * dt
    K4 = G(Xsm[j][0] + K3, Xsm[j][1] + K3) * dt
    mo[j+1] = mo[j] + 1.0 / 6.0 * (K1 + 2 * K2 + 2 * K3 + K4)
t2 = np.linspace(0.0, T, N)

# Plot results
fig = plt.figure(figsize=(12, 10))

# Plot u as a function of time
ax1 = fig.add_subplot(4, 2, 1)
ax1.plot(t1, det[:, 0], 'b')
ax1.set_xlabel('time')
ax1.set_ylabel('OCs')
plt.title('Deterministic Model')

# Plot v as a function of time
ax2 = fig.add_subplot(4, 2, 3)
ax2.plot(t1, det[:, 1], 'r')
ax2.set_xlabel('time')
ax2.set_ylabel('OBs')

# Plot w as a function of time
ax3 = fig.add_subplot(4, 2, 5)
ax3.plot(t1, det[:, 2], 'g')
ax3.set_xlabel('time')
ax3.set_ylabel('CCs')

# Plot z as a function of time
ax4 = fig.add_subplot(4, 2, 7)
ax4.plot(t1, det[:, 3], 'm')
ax4.set_xlabel('time')
ax4.set_ylabel('MO')

# Plot u as a function of time
ax5 = fig.add_subplot(4, 2, 2)
ax5.plot(t2, Xsm[:, 0], 'b')
ax5.set_xlabel('time')
ax5.set_ylabel('OCs')
plt.title('Modelo Estocastico')

# Plot v as a function of time
ax5 = fig.add_subplot(4, 2, 4)
ax5.plot(t2, Xsm[:, 1], 'r')
ax5.set_xlabel('time')
ax5.set_ylabel('OBs')

# Plot w as a function of time
ax6 = fig.add_subplot(4, 2, 6)
ax6.plot(t2, Xsm[:, 2], 'g')
ax6.set_xlabel('time')
ax6.set_ylabel('CCs')

# Plot z as a function of time
ax8 = fig.add_subplot(4, 2, 8)
ax8.plot(t2, mo, 'm')
ax8.set_xlabel('time')
ax8.set_ylabel('MO')

plt.tight_layout()
file_name = time.clock()
file_name = str(file_name)
file_name += '.png'
plt.savefig(file_name)
