import numpy as np
import matplotlib as mpl
mpl.use('PS')
from matplotlib import rcParams
import matplotlib.pyplot as plt
from StoPLBRJerezChen import StoPLBRM
import time
w = 1.0
fig_width_mm = 120 * w  # 120 mm
inches_per_mm = 1.0 / 25.4  # Convert mm to inch
golden_mean = (np.sqrt(5) - 1.0) / 2.0  # Aesthetic ratio
fig_width = fig_width_mm * inches_per_mm  # width in inches
fig_height = fig_width * golden_mean  # height in inches
fig_size = [fig_width, fig_height]
params = {
    'axes.labelsize': 10,
    'font.size': 10,
    'legend.fontsize': 8,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'text.usetex': True,
    'figure.figsize': fig_size}
rcParams.update(params)
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
U0 = [10, .7]
k = 20
p = 0
r = p
T0 = 0.0
T = 650 * 48
eps = 1e-4
color_pallet = ['#588C7E', '#F2E394', '#F2AE72', '#D96459', '#8C4646']
# seeds = []
seeds = [109966818, 115396582, 114793524]
color1 = color_pallet[0]
color2 = color_pallet[1]
color3 = color_pallet[2]
color4 = color_pallet[3]
color5 = color_pallet[4]

det_color = color1
seed = 115396582
StoPlbrmJC = StoPLBRM()
StoPlbrmJC.initialize_mesh(k, p, r, T0, T)
# u_ssls_ref = StoPlbrmJC.ssls(seed, [1, 1], fn=1.0)
for i in seeds:
    seed = i
    U0 = [10.0, .7]
    plt.close()
    eps_prefix = str(eps) + '-'
    file_name1 = 'OB-OC-' + eps_prefix + str(seed)
    file_name1 += '.eps'
    fig1, (ax1, ax2) = plt.subplots(nrows=2)
    ax1 = plt.subplot(211)
    ax2 = plt.subplot(212)
    for j in np.arange(3):
        np.random.seed(i)
        # seed = np.random.randint(109730187, 123456789)
        # np.random.seed(seed)
        # seeds.append(seed)-----------------------------------------
        U0 += j * eps
        StoPlbrmJC.u_zero = U0
        print U0
        StoPlbrmJC.set_parameters_sto_plbrm(a1, b1, a2, b2, 1.0, gamma2,
                                            gamma1, 1.0, k1, k2, sigma, U0)
        #
        # seed = j
        # ax1.set_xlabel(r'$t$ (days)')
        ax1.set_ylabel(r'Number of OCs')
        ax1.set_xlim([-200, 650 * 48])
        ax1.grid(False)
        ax2.set_xlabel(r'$t$ (days)')
        ax2.set_ylabel(r'Number of OBs')
        ax2.set_xlim([-200, 650 * 48])
        ax2.grid(False)
        # -----------------
        # StoPlbrmJC.noise_update(seed)
        u_ssls = StoPlbrmJC.ssls(seed, [1, 1], fn=1.0)
        # ---------------Long Path----------------------------------------
        t = StoPlbrmJC.t_k
        #
        sto_color = color_pallet[j]
        #
        label_legend = 'eps:=' + str(j * eps)
        ax1.plot(t, u_ssls[:, 0],
                 color=sto_color,
                 lw=.5,
                 linestyle='-',
                 label=label_legend
                 )
        ax2.plot(t, u_ssls[:, 1],
                 color=sto_color,
                 lw=.5,
                 linestyle='-',
                 label=label_legend
                 )
        np.save('OneLongPathSolutionSto' + eps_prefix + '.npy',
            np.transpose(np.array([t, u_ssls[:, 0], u_ssls[:, 1]])))
    plt.legend(loc=0)
    plt.tight_layout()
    plt.savefig(file_name1, resolution=1000)
