import numpy as np
import matplotlib.pyplot as plt
from StoPLBRJerezChen import StoPLBRM
# Model parameters
#  k1 = 0.093,  k2 = .00032, r=290.625
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
a1 = 0.15    # Model's parameters ai = \alpha_i
a2 = 0.16
b1 = 0.2
b2 = 0.02
ns = 0.1    # Noise intensity
gamma1 = -0.837
gamma2 = 0.88
sigma = np.array([ns * b1, ns * b2])
U0 = [10, 1.0]    # Initial Conditions
# Stencil Parameters
k = 18
p = 0
r = p
T0 = 0.0
T = 650 * 48  # almost 48 periods
Ji1 = ((b1 + 0.5 * (sigma[0] ** 2)) / a1) ** (1.0 / gamma1)
Ji2 = ((b2 + 0.5 * (sigma[1] ** 2)) / a2) ** (1.0 / gamma2)
print '\n\n'
condition_H1 = (np.sign(gamma1) < 0.0) and (np.sign(gamma2) > 0.0)
condition_H2 = (np.abs(gamma1) <= np.sign(gamma2) )
condition_H3 = a1 * gamma2 < a2 * np.abs(gamma1)
print'\t\t Xi1, Xi2:= '+'['+str(Ji1)+', '+str(Ji2)+']'
if (condition_H1 and (condition_H2 and condition_H3)):
    print '\n\t\t\tH3 holds =) \n\n'
else:
    print '\n\t\t\tH3 not holds =(\n\n'
    print 'p1:=', a1 * gamma2, 'p2:=', a2 * np.abs(gamma1)

StoPlbrmJC = StoPLBRM()
StoPlbrmJC.initialize_mesh(k, p, r, T0, T)
k1 = 0.093
k2 = 0.00032
ussls_det = StoPlbrmJC.ssls_det()
ussls = StoPlbrmJC.ssls()
t = StoPlbrmJC.t_k
ones = np.ones(t.shape[0])
j = 0
StoPlbrmJC.set_parameters_sto_plbrm(a1, b1, a2, b2, 1.0, gamma2, gamma1, 1.0, k1, k2, sigma, U0)
z = StoPlbrmJC.bone_mass()
np.save('BoneMass.npy', np.transpose(np.array([t[0:-1:100], StoPlbrmJC.z[0:-1:100]])))
np.save("OneLongPathSolutionSto.npy", np.transpose(np.array([t, StoPlbrmJC.u_ssls[0, :],
                                                             StoPlbrmJC.u_ssls[:, 1]])))