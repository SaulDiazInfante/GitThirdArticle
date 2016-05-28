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
T = 650 * 4  # almost 48 periods
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
bunche = 10
k2_list = []
r_list = []
r_det = 10  # 603
StoPlbrmJC = StoPLBRM()
StoPlbrmJC.initialize_mesh(k, p, r, T0, T)
i = 0
k1 = 0.093
k2 = 0.00032
StoPlbrmJC.print_progress(i, bunche, 'Progress:', 'Complete',
                          bar_length=50, ratio=True)
ussls_det = StoPlbrmJC.ssls_det()
t = StoPlbrmJC.t_k
ones = np.ones(t.shape[0])
j = 0
for i in np.arange(bunche):
    """
    r_seek = r_det + r_det * 10 * np.random.rand()
    k2 = k2_mean + k2_mean * 10 * np.random.rand()
    """
    # r_seek = np.random.uniform(79,79.5)
    # k2 = np.random.uniform(69e-6, 7e-5)
    r_seek = r
    k1 += 1e-3
    #k2 += 1e-5
    # k2= k2_mean[i]
    # k1 = r_seek * k2
    StoPlbrmJC.set_parameters_sto_plbrm(a1, b1, a2, b2, 1.0, gamma2,
                                        gamma1, 1.0, k1, k2, sigma, U0)
    # seed = np.random.random_integers(1, 123456789)
    # np.random.seed(123456789)
    # StoPlbrmJC.NoiseUpdate(seed, flag=1)
    #
    # ussls = StoPlbrmJC.SSLS(seed, U0, fn=1.0)
    # StoPlbrmJC.SaveData()
    z = StoPlbrmJC.bone_mass()
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
    tag = 'k1:='+str(k1) + ', k2:='+str(k2)
    file_name = 'parameter_output' + str(j) + '.png'
    plt.title(tag)
    plt.plot(t[0: -1: 100],  z[0:-1:100])
    plt.plot(t[0: -1: 100], 100 * ones[0: -1: 100])
    plt.plot(t[0: -1: 100], 75 * ones[0: -1: 100])
    #plt.legend(loc=0)
    plt.savefig(file_name)
    j += 1
    StoPlbrmJC.print_progress(i + 1, bunche, 'Progress:', 'Complete',
                              bar_length=50, ratio=True)

# plt.figure()
# plt.plot(ussls_det[:, 0], ussls_det[:, 1], alpha=0.7, label='SSLS')
# plt.legend(loc=0)
# plt.show()
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
                                        r_list.pop())'''
