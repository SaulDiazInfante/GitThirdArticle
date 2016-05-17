import numpy as np
import matplotlib.pyplot as plt
from StoPLBRJerezChen import StoPLBRM
# Model parameters
#  k1 = 0.0029055,  k2 = 3.9e-6, r=603
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
a1 = 0.3    # Model's parameters ai = \alpha_i
a2 = 0.18
b1 = 0.2
b2 = 0.02
ns = 0.1    # Noise intensity
gamma1 = -0.9
gamma2 = 0.5
sigma = np.array([ns*b1, ns*b2])
U0 = [10, 0.7]    # Initial Conditions
# Stencil Parameters
k = 16
p = 0
r = p
T0 = 0.0
T = 650 * 8 # almost 48 periods
Ji1 = ((b1 + 0.5 * (sigma[0] ** 2)) / a1) ** (1.0 / gamma1)
Ji2 = ((b2 + 0.5 * (sigma[1] ** 2)) / a2) ** (1.0 / gamma2)
print '\n\n'
print'\t\t Xi1, Xi2:= '+'['+str(Ji1)+', '+str(Ji2)+']'
if a1 * gamma2 < a2 * np.abs(gamma1):
    print '\n\t\t\tH3 holds =) \n\n'
else:
    print '\n\t\t\tH3 not holds =(\n\n'
bunche = 1
k2_list = []
r_list = []
r_det = 10  # 603
StoPlbrmJC = StoPLBRM()
StoPlbrmJC.initialize_mesh(k, p, r, T0, T)
i = 0
k2_mean = 39e-5
StoPlbrmJC.print_progress(i, bunche, 'Progress:', 'Complete',
                          bar_length=50, ratio=True)
ussls_det = StoPlbrmJC.ssls_det()
ones = np.ones(ussls_det.shape[0])
t = StoPlbrmJC.t
j = 0
k2 = 3.9e-6
r=603

for i in np.arange(bunche):
    """
    r_seek = r_det + r_det * 10 * np.random.rand()
    k2 = k2_mean + k2_mean * 10 * np.random.rand()
    """
    # r_seek = np.random.uniform(79,79.5)
    # k2 = np.random.uniform(69e-6, 7e-5)
    r_seek = r
    #k2= k2_mean[i]
    k1 = r_seek * k2
    StoPlbrmJC.set_parameters_sto_plbrm(a1, b1, a2, b2, 1.0,
                                     gamma2, gamma1, 1.0, k1, k2,
                                    sigma, U0)
    # seed = np.random.random_integers(1, 123456789)
    # np.random.seed(123456789)
    # StoPlbrmJC.NoiseUpdate(seed, flag=1)
    #
    # ussls = StoPlbrmJC.SSLS(seed, U0, fn=1.0)
    # StoPlbrmJC.SaveData()
    StoPlbrmJC.bone_mass()
    condition_1 = not (StoPlbrmJC.z > 1.04).any()
    condition_2 = not (StoPlbrmJC.z < .73).any()
    if condition_1 and condition_2:
        print '\n'
        print u"\tGood candidates on run {0:d}:".format(i)
        print '\tr:='+str(r_seek) + '\t k2:='+str(k2)
    k2_list.append(k2)
    r_list.append(r_seek)
    np.save('BoneMass.npy',
            np.transpose(
                np.array([t[0:-1:100], StoPlbrmJC.z[0:-1:100]])))
    tag = 'r:='+str(r_seek) + ', k2:='+str(k2)
    file_name = 'parameter_output' + str(j) + '.png'
    plt.plot(t[0: -1: 100], 100 * StoPlbrmJC.z[0:-1:100], label=tag)
    plt.plot(t[0: -1: 100], ones[0: -1: 100] * 100)
    plt.plot(t[0: -1: 100], .75 * 100 * ones)
    plt.legend(loc=0)
    plt.savefig(file_name)
    plt.close()
    j += 1
    StoPlbrmJC.print_progress(i + 1, bunche, 'Progress:', 'Complete',
                              bar_length=50, ratio=True)
'''
plt.close()
plt.figure()
urk_max0 = np.max(urk[:, 0])
urk_max1 = np.max(urk[:, 1])
ussls_max0 = np.max(ussls[:, 0])
ussls_max1 = np.max(ussls[:, 1])
plt.plot(urk[:, 0]/urk_max0, urk[:, 1]/urk_max1, alpha=0.5, label='RK')
plt.plot(ussls[:, 0]/ussls_max0, ussls[:, 1]/ussls_max1, alpha=0.5,
         label='SSLS')
plt.legend(loc=0)
plt.show()
'''
a = np.array(k2_list)
a = a.reshape([a.shape[0], 1])
b = np.array(r_list)
b = b.reshape([b.shape[0], 1])
c = np.concatenate((a, b), axis=1)
np.savetxt("parametersBoneMass.txt", c)
print "\t\t\t\t {0:d} runs".format(bunche)
print "\t-------------------------------------------------------------"
print "\t i \t\tk2 \t\t\tr"
print "\t-------------------------------------------------------------"
size = len(k2_list)
for i in np.arange(size):
    print '\t{}\t {} \t {}'.format(size - 1 - i, k2_list.pop(), \
                                        r_list.pop())