import numpy as np
import matplotlib.pyplot as plt
from StoPLBRJerezChen import StoPLBRM
# Model parameters
#  k1 = 0.0029055,  k2 =3.9e-6, r=603
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
a1 = 0.3    # Model's parameters ai = \alpha_i
a2 = 0.18
b1 = 0.2
b2 = 0.02
ns = 0.1    # Noise intensity
gamma1 = -0.9
gamma2 = 0.5
#
sigma = np.array([ns*b1, ns*b2])
U0 = [10, 0.7]    # Condicion inicial
# Stencil Parameters
k = 5
p = 0
r = p
T0 = 0.0
T = 48 * 650    # almost 48 periods
Ji1 = ((b1 + 0.5 * (sigma[0] ** 2)) / a1) ** (1.0 / gamma1)
Ji2 = ((b2 + 0.5 * (sigma[1] ** 2)) / a2) ** (1.0 / gamma2)
print '\n\n'
print'\t\t Xi1, Xi2:= '+'['+str(Ji1)+', '+str(Ji2)+']'
if a1 * gamma2 < a2 * np.abs(gamma1):
    print '\n\t\t\tH3 holds =) \n\n'
else:
    print '\n\t\t\tH3 not holds =(\n\n'
bunche = 5
k2_list = []
r_list = []
r_det = 79.5
StoPlbrmJC = StoPLBRM()
StoPlbrmJC.initialize_mesh(k, p, r, T0, T)
i = 0
k2_mean = 0.000004
StoPlbrmJC.print_progress(i, bunche, 'Progress:', 'Complete',
                          bar_length=50, ratio=True)
for i in np.arange(bunche):
    # StoPlbrmJC.NoiseUpdate(123456789)
    r_seek = r_det + r_det * 0.1 * np.random.rand()
    # r_sto = rr  # + np.random.randint(1000)
    k2 = k2_mean + k2_mean * .01 * np.random.rand()
    k1 = r_seek * k2
    StoPlbrmJC.set_parameters_sto_plbrm(a1, b1, a2, b2, 1.0,
                                     gamma2, gamma1, 1.0, k1, k2,
                                    sigma, U0)
    t = StoPlbrmJC.t[0: -1: np.int64(StoPlbrmJC.R)]
    t = t.reshape([StoPlbrmJC.t[0: -1: np.int64(StoPlbrmJC.R)].shape[0], 1])
    urk = StoPlbrmJC.rk()
    # seed = np.random.random_integers(1, 123456789)
    # np.random.seed(123456789)
    # StoPlbrmJC.NoiseUpdate(seed, flag=1)
    # ustk = StoPlbrmJC.SSLS(123456789, U0, fn=0.0)
    # ussls = StoPlbrmJC.SSLS(seed, U0, fn=1.0)
    # StoPlbrmJC.SaveData()
    StoPlbrmJC.bone_mass()
    condition_1 = np.abs(StoPlbrmJC.z.max() - 1) < 0.05
    condition_2 = np.abs(StoPlbrmJC.z.min() - .75) < 0.05
    StoPlbrmJC.print_progress(i + 1, bunche, 'Progress:', 'Complete',
                              bar_length=50, ratio=True)
    if condition_1 and condition_2:
        print '\n'
        print u"\tGood candidates on run {0:d}:".format(i)
        print '\tr:='+str(r_seek) + '\t k2:='+str(k2)
    k2_list.append(k2)
    r_list.append(r_seek)
    plt.plot(t, 100 * StoPlbrmJC.z, label=str(i))
    '''
    np.save('BoneMass.npy',
            np.transpose(np.array(
                [t[:, 0], StoPlbrmJC.z, StoPlbrmJC.zSto])
                        )
            )
    plt.close()
    plt.figure()
    plt.plot(t, ussls[:,0])
    plt.plot(t, ustk[:, 0])
    plt.ylabel('OC')
    #
    plt.figure()
    plt.plot(t, ussls[:, 1])
    plt.plot(t, Ustk[:, 1])
    plt.ylabel('OB')
    '''
    # plt.figure()
    # plt.plot(t, 100 * StoPlbrmJC.zSto, label=str(i))

plt.plot(t, 100 * np.ones(StoPlbrmJC.z.shape[0]))
plt.legend(loc=0)
plt.savefig('parameter_output.png')
plt.show()
a = np.array(k2_list)
a = a.reshape([a.shape[0], 1])
b = np.array(r_list)
b = b.reshape([b.shape[0], 1])
c = np.concatenate((a, b), axis=1)
np.savetxt("parametersBoneMass.txt", c)
print "----------------------------------------------------------------"
print "\t\t {0:d} runs".format(bunche)
print "----------------------------------------------------------------"
for i in np.arange(len(k2_list)):
    print '{}\t {} \t {}'.format(i, k2_list.pop(), r_list.pop())
