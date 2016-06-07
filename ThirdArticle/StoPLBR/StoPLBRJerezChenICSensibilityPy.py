import numpy as np
import matplotlib as mpl
mpl.use('PS')
from matplotlib import rcParams
import matplotlib.pyplot as plt
from StoPLBRJerezChen import StoPLBRM
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

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
# eps = [.0001, .001, .01]
eps =[0.01, 0.1, 1.0]
color_pallet = ['#588C7E', '#F2E394', '#F2AE72', '#D96459', '#8C4646']
# seeds = []
# seeds = [109966818, 115396582, 114793524]
seeds = [114793524]
color1 = color_pallet[0]
color2 = color_pallet[1]
color3 = color_pallet[2]
color4 = color_pallet[3]
color5 = color_pallet[4]
det_color = color1
seed = 114793524
StoPlbrmJC = StoPLBRM()
StoPlbrmJC.initialize_mesh(k, p, r, T0, T)
rg_old_state = StoPlbrmJC.old_rg_state
#
# for i in seeds:
#     seed = i
U0 = [10.0, .7]
print '\n\t\t seed: {:,}'.format(seed)
plt.close()
eps_prefix = str(eps) + '-'
file_name1 = 'OB-OC-' + eps_prefix + str(seed)
file_name1 += '.eps'
t = StoPlbrmJC.t_k
u_ssls_ref = StoPlbrmJC.ssls(seed, [1, 1], fn=1.0)
k = 0
sto_color = color_pallet[k]
fig1, (ax1, ax2) = plt.subplots(nrows=2)
label_legend = r'$\epsilon:=$' + str(0)
ax1 = plt.subplot(211)
ax2 = plt.subplot(212)
ax1.set_xlabel(r'$t$ (days)')
ax1.set_ylabel(r'Number of OCs')
ax1.set_xlim([-200, 650 * 48])
ax1.grid(False)
ax1.plot(t, u_ssls_ref[:, 0],
         color=sto_color,
         lw=.5,
         linestyle='-',
         label=label_legend
         )
ax2.set_xlabel(r'$t$ (days)')
ax2.set_ylabel(r'Number of OBs')
ax2.set_xlim([-200, 650 * 48])
ax2.grid(False)
ax2.plot(t, u_ssls_ref[:, 1],
         color=sto_color,
         lw=.5,
         linestyle='-',
         label=label_legend
         )
axins_a = inset_axes(ax2, 1, 0.5, loc=3, bbox_to_anchor=(0.5, 0.25), bbox_transform=ax2.figure.transFigure)
t1, t2 = 5000, 8000
y1, y2 = -10, 1300
# t1, t2 =
axins_a.set_xlim(t1, t2)
axins_a.set_ylim(y1, y2)

axins_a.plot(t, u_ssls_ref[:, 1],
             color=sto_color,
             lw=1,
             linestyle='-'
             )
plt.yticks([y1, y2], visible=True)
plt.xticks([t1, t2], visible=True)
# plt.show()
#    ##################################################3

for j in eps:
    U0 = np.array([10.0, .7])
    # np.random.set_state(rg_old_state)
    # StoPlbrmJC.noise_update()
    U0 += j * np.array([1.0, 1.0])
    StoPlbrmJC.u_zero = U0
    print '\n\t\tu_zero = [{:.5f}, {:.5f}]'.format(StoPlbrmJC.u_zero[0],StoPlbrmJC.u_zero[1])
    # StoPlbrmJC.set_parameters_sto_plbrm(a1, b1, a2, b2, 1.0, gamma2,
    #                                    gamma1, 1.0, k1, k2, sigma, U0)
    #
    # seed = j
    # -----------------
    # StoPlbrmJC.noise_update(seed)
    u_ssls = StoPlbrmJC.ssls(seed, [1, 1], fn=1.0)
    # ---------------Long Path----------------------------------------

    #
    k += 1
    sto_color = color_pallet[k]
    label_legend = r'$\epsilon:=$' + str(j)
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
    axins_a.plot(t, u_ssls[:, 1],
                 color=sto_color,
                 lw=1,
                 linestyle='-'
                 )
    np.save('OneLongPathSolutionSto' + eps_prefix + '.npy',
             np.transpose(np.array([t, u_ssls[:, 0], u_ssls[:, 1]])))


#
plt.yticks([y1-100, y2-100], visible=True)
plt.xticks([t1-200, t2-200], visible=True)
mark_inset(
    ax2,
    axins_a,
    loc1=3,
    loc2=4,
    fc="none",
    ec="0.5"
)


ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.315), ncol=4,
           fancybox=False, shadow=False)
plt.tight_layout()
plt.savefig(file_name1, resolution=1000)
