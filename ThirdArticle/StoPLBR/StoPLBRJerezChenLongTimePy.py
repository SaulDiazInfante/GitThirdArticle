import numpy as np
import matplotlib.pyplot as plt
from StoPLBRJerezChen import StoPLBRM
import time
#     Model parameters  ------------------------------------------------
a1 = 0.15
a2 = 0.16
b1 = 0.2
b2 = 0.02
ns = 0.1    # Noise intensity
gamma1 = -0.837
gamma2 = 0.88
sigma = np.array([ns*b1, ns*b2])
k1 = 0.03
k2 = 0.0017
#
Ji1 = ((b1 + 0.5 * (sigma[0] ** 2)) / a1) ** (1.0 / gamma1)
Ji2 = ((b2 + 0.5 * (sigma[1] ** 2)) / a2) ** (1.0 / gamma2)
# Stencil Parameters
U0 = [10, 0.700]
k = 20
p = 0
r = p
T0 = 0.0
T = 650 * 48
#seeds = []
seeds = [109966818, 123225530, 114068981, 115396582, 112428882, 114793524]
for j in seeds:
    # seed = np.random.randint(109730187, 123456789)
    # seed = 118329736
    # np.random.seed(seed)
    # seeds.append(seed)
    #
    seed = j
    np.random.seed(seed)
    StoPlbrmJC = StoPLBRM()
    StoPlbrmJC.initialize_mesh(k, p, r, T0, T)
    # ----------------------------------------------------------

    StoPlbrmJC.set_parameters_sto_plbrm(a1, b1, a2, b2, 1.0, gamma2,
                                        gamma1, 1.0, k1, k2, sigma, U0)
    StoPlbrmJC.noise_update(seed)
    Ussls = StoPlbrmJC.ssls(seed, [1, 1], fn=1.0)
    stoU1 = Ussls[:, 0]
    stoU2 = Ussls[:, 1]
    # ---------------Long Path----------------------------------------
    seed_prefix = str(seed)
    file_name1 = 'U1'+ seed_prefix
    file_name2 = 'U2' + seed_prefix
    file_name1 += '.png'
    file_name2 += '.png'
    t = StoPlbrmJC.t_k
    #
    plt.figure()
    plt.title('seed: '+ str(seed))
    plt.plot(t, stoU1)
    plt.savefig(file_name1)
    plt.close()
    #
    plt.figure()
    plt.title('seed: ' + str(seed))
    plt.plot(t, stoU2)
    plt.savefig(file_name2)
    plt.close()
    #
    np.save('OneLongPathSolutionSto' + seed_prefix + '.npy',
             np.transpose(np.array([t, stoU1[:], stoU2[:]])))
#   np.save('OneLongPathSolutionDet.npy',
#         np.transpose(np.array([t[:, 0], U1[:], U2[:]])))
