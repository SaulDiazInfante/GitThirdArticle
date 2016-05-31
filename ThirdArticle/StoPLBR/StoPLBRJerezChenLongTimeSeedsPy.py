import numpy as np
import matplotlib as mpl
mpl.use('PS')
from matplotlib import rcParams
import matplotlib.pyplot as plt
from StoPLBRJerezChen import StoPLBRM
import time
w = 1.0
fig_width_mm = 170 * w  # 120 mm
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
color_pallet = ['#588C7E', '#F2E394', '#F2AE72', '#D96459', '#8C4646']
# seeds = []
seeds = [109966818, 115396582, 114793524]
color1 = color_pallet[0]
color2 = color_pallet[1]
color3 = color_pallet[2]
color4 = color_pallet[3]
color5 = color_pallet[4]
det_color = color1
i = 1
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
    u_ssls = StoPlbrmJC.ssls(seed, [1, 1], fn=1.0)
    # ---------------Long Path----------------------------------------
    seed_prefix = str(seed)
    file_name1 = 'OB-OC-' + seed_prefix
    file_name1 += '.eps'
    t = StoPlbrmJC.t_k
    #
    fig1, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
    #
    sto_color = color3
    ax1 = plt.subplot(221)
    # ax1.set_xlabel(r'$t$ (days)')
    ax1.set_ylabel(r'Number of OCs')
    ax1.set_xlim([-200, 650 * 48])
    ax1.plot(t, u_ssls[:, 0],
             color=sto_color,
             lw=.5,
             linestyle='-',
             label='Deterministic'
             )
    ax1.grid(False)
    ax3 = plt.subplot(223)
    # ax1.set_xlabel(r'$t$ (days)')
    ax3.set_ylabel(r'Number of OCs')
    ax3.set_xlim([-200, 650 * 48])
    ax3.plot(t, u_ssls[:, 1],
             color=sto_color,
             lw=.5,
             linestyle='-',
             label='Deterministic'
             )
    ax3.grid(False)
    #
    StoPlbrmJC.noise_update()
    u_ssls = StoPlbrmJC.ssls(seed, [1, 1], fn=1.0)
    sto_color = color5
    ax2 = plt.subplot(222)
    ax2.set_xlabel(r'$t$ (days)')
    ax2.set_ylabel(r'Number of OBs')
    ax2.set_xlim([-200, 650 * 48])
    ax2.plot(t, u_ssls[:, 0],
             color=sto_color,
             lw=.5,
             linestyle='-',
             label='Deterministic'
             )
    ax2.grid(False)

    ax4 = plt.subplot(224)
    ax4.set_xlabel(r'$t$ (days)')
    ax4.set_ylabel(r'Number of OBs')
    ax4.set_xlim([-200, 650 * 48])
    ax4.plot(t, u_ssls[:, 1],
             color=sto_color,
             lw=.5,
             linestyle='-',
             label='Deterministic'
             )
    ax4.grid(False)
    plt.tight_layout()
    plt.savefig(file_name1, resolution=1000)
    #
    np.save('OneLongPathSolutionSto' + seed_prefix + '(' + str(i) + ')' + '.npy',
             np.transpose(np.array([t, u_ssls[:, 0], u_ssls[:, 1]])))
    i += 1
#   np.save('OneLongPathSolutionDet.npy',
#         np.transpose(np.array([t[:, 0], U1[:], U2[:]])))
