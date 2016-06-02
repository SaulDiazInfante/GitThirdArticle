import numpy as np
import matplotlib.pyplot as plt
from StoPLBRJerezChen import StoPLBRM
# Model parameters
#  k1 = 0.093,  k2 = .00032, r=290.625
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
# Stencil Parameters
k = 20
p = 0
r = p
T0 = 0.0
T =  2**11# 650 * 4  # almost 48 periods


StoPlbrmJC = StoPLBRM()
StoPlbrmJC.load_parameters()
StoPlbrmJC.initialize_mesh(k, p, r, T0, T)
StoPlbrmJC.noise_update()
gamma1 = StoPlbrmJC.g21
gamma2 = StoPlbrmJC.g12
Ji1 = StoPlbrmJC.ji_1
Ji2 = StoPlbrmJC.ji_2
print '\n\n'
condition_H1 = (np.sign(gamma1) < 0.0) and (np.sign(gamma2) > 0.0)
condition_H2 = (np.abs(gamma1) <= np.sign(gamma2) )
condition_H3 = StoPlbrmJC.a1 * gamma2 < StoPlbrmJC.a2 * np.abs(gamma1)
print'\t\t Xi1, Xi2:= '+'['+str(Ji1)+', '+str(Ji2)+']'
if condition_H1 and (condition_H2 and condition_H3):
    print '\n\t\t\tH3 holds =) \n\n'
else:
    print '\n\t\t\tH3 not holds =(\n\n'
    print 'p1:=', StoPlbrmJC.a1 * gamma2, 'p2:=', StoPlbrmJC.a2 * np.abs(gamma1)

t = StoPlbrmJC.t_k

u_det = StoPlbrmJC.rk()
np.save("OneLongPathSolutionDet.npy", np.transpose(np.array([t, u_det[:, 0], u_det[:, 1]])))
StoPlbrmJC.noise_update()
# ussls = StoPlbrmJC.ssls()
np.save("OneLongPathSolutionSto.npy", np.transpose(np.array([t, StoPlbrmJC.u_ssls[:, 0],
                                                             StoPlbrmJC.u_ssls[:, 1]])))
ones = np.ones(t.shape[0])
z = StoPlbrmJC.bone_mass()
# np.save('BoneMass.npy', np.transpose(np.array([t[0:-1:100], StoPlbrmJC.z[0:-1:100]])))
plt.plot(t, z)
plt.show()
