"""
# ===============================================================================
# Code for create the eps-figures of the Article:
#    Jerez, Diaz-Infante
# INPUT:
#   Binary data files:
#   SolutionDet.txt, SolutionSto.txt
#   SolutionStoAs.txt, SolutionDetAs.txt
#   U1PathsShortTime.npy, U2PathsShortTime.npy
#   U1PathsLongTime.npy, U2PathsLongTime.npy
# OUTPUT:
#   CellsShortTime.eps, ,
#   OsteoclastShortTime.eps, OsteoclastLongTime.eps
#   OsteoblastShortTime.eps, OsteoblastLongTime.eps
#   PhasePotrait.eps, PhasePotrait3d.eps
# ===============================================================================
"""
import numpy as np
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
from StoPLBRJerezChen import StoPLBRM
# from matplotlib.pyplot import flag
# from matplotlib import colors
from matplotlib import rcParams
from mpl_toolkits.mplot3d import Axes3D
# from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
# import matplotlib.animation as animation
# ---------------------------------------------------------------------


def fig_cells():
    urk = np.load('SolutionDet.npy')
    ussls = np.load('SolutionSto.npy')
    # uas = np.loadtxt('SolutionStoAs.txt')
    # t = ussls[:, 0]
    # tas = uas[:, 0]
    # Pots Vs Time
    plt.close()
    w = 1.0
    fig_width_pt = 538.58 * w    # \showthe\columnwidth(190mm)
    inches_per_pt = 1.0 / 72.27    # Convert pt to inch
    golden_mean = (np.sqrt(5) - 1.0) / 2.0    # Aesthetic ratio
    fig_width = fig_width_pt * inches_per_pt    # width in inches
    fig_height = fig_width * golden_mean	   # height in inches
    fig_size = [fig_width, fig_height]
    params = {
        'backend': 'ps',
        'axes.labelsize': 10,
        'text.fontsize': 10,
        'legend.fontsize': 8,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'text.usetex': True,
        'figure.figsize': fig_size
    }
    rcParams.update(params)
#
    fig1, ax1 = plt.subplots()
    # plt.title('Jerez-Chen Model Experimental Parameters')
    color1 = color_pallet[0]
    color2 = color_pallet[1]
    color3 = color_pallet[2]
    color4 = color_pallet[3]
    oc_line = mlines.Line2D(
        [], [],
        color=color1,
        linestyle='-',
        lw=2, marker='',
        markersize=1
    )
    ob_line = mlines.Line2D(
        [], [],
        color=color2,
        linestyle='-',
        lw=2, marker='',
        markersize=1
    )
    soc_line = mlines.Line2D(
        [], [],
        color=color3,
        linestyle='-.',
        lw=4,
        marker='',
        markersize=1
    )
    sob_line = mlines.Line2D(
        [], [],
        color=color4,
        linestyle='-.',
        lw=4,
        marker='',
        markersize=1
    )
    #
    y1max = np.max([ussls[:, 1].max(), urk[:, 1].max()])
    y2max = np.max([ussls[:, 2].max(), urk[:, 2].max()])
    ax1.set_xlabel('time (days)')
    ax1.set_ylabel('OCs (u)', color='k')
    ax1.set_xlim([-200, 5400])
    ax1.set_ylim([-1, 1.05 * y1max])
    ax2 = ax1.twinx()
    ax2.set_ylim([-267, 1.05 * y2max])
    ax2.set_xlim([-200, 5400])
    ax1.plot(urk[:, 0], urk[:, 1],
             color=color1,
             linestyle='-',
             lw=2,
             marker='',
             label='osteoclast'
             )
    ax1.plot(ussls[:, 0], ussls[:, 1],
             color=color3,
             linestyle='-.',
             lw=3,
             marker='',
             label='Sto osteoclast'
             )
    #
    ax2.plot(urk[:, 0], urk[:, 2],
             color=color2,
             lw=2,
             linestyle='-',
             marker='',
             label='osteoblast'
             )
    ax2.plot(ussls[:, 0], ussls[:, 2],
             color=color4,
             linestyle='-.',
             lw=3,
             marker='',
             label='Sto osteoblast'
             )
    ax2.set_xlabel('time (days)')
    ax2.set_ylabel(r'OBs (v)')
    ax2.grid(True)
    plt.legend(
        [oc_line, ob_line, soc_line, sob_line],
        ["OCs", "OBs", "stoOCs", "stoOBs"],
        bbox_to_anchor=(0., 1.02, 1., .5),
        loc=3,
        ncol=4,
        mode="expand",
        borderaxespad=0
    )
    plt.savefig("CellsTime.eps")
    plt.show()
#


def fig_oc():
    urk = np.load('SolutionDet.npy')
    ussls = np.load('SolutionSto.npy')
    detuas = np.loadtxt('SolutionDetAs.txt')
    uas = np.loadtxt('SolutionStoAs.txt')
    # t = ussls[:, 0]
    tas = uas[:, 0]
    # Pots Vs Time
    plt.close()
    w = 1.0
    fig_width_pt = 538.58 * w    # \columnwidth(190mm)
    inches_per_pt = 1.0 / 72.27					# Convert pt to inch
    golden_mean = (np.sqrt(5) - 1.0) / 2.0   	# Aesthetic ratio
    fig_width = fig_width_pt * inches_per_pt 	# width in inches
    fig_height = fig_width * golden_mean	   		# height in inches
    fig_size = [fig_width, fig_height]
    params = {
        'backend': 'ps',
        'axes.labelsize': 10,
        'text.fontsize': 10,
        'legend.fontsize': 8,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'text.usetex': True,
        'figure.figsize': fig_size
    }
    rcParams.update(params)
#
##
    color1 = color_pallet[0]
    # color2 = colorPallete[1]
    color3 = color_pallet[2]
    # color4 = colorPallete[3]
    fig1, ax1 = plt.subplots()
    oc_line = mlines.Line2D([], [],
                            color=color1,
                            linestyle='-',
                            lw=2,
                            marker='',
                            markersize=1
                            )

    soc_line = mlines.Line2D([], [],
                             color=color3,
                             linestyle='-.',
                             lw=4, marker='',
                             markersize=1
                             )
    #
    ax1.set_xlabel('time (days)')
    ax1.set_ylabel('OCs (u)', color='k')
    y1max = np.max([ussls[:, 1].max(), urk[:, 1].max()])
    ax1.set_xlim([-200, 5400])
    ax1.set_ylim([-2, 1.05 * y1max])
    ax1.grid(True)
    ax1.plot(urk[:, 0], urk[:, 1],
             color=color1,
             linestyle='-',
             lw=2,
             marker='',
             label='osteoclast'
             )
    ax1.plot(ussls[:, 0], ussls[:, 1],
             color=color3,
             linestyle='-.',
             lw=3,
             marker='',
             label='Sto Osteoclast'
             )
    plt.legend(
        [oc_line, soc_line],
        ["OCs", "stoOCs"],
        bbox_to_anchor=(0., 1.02, 1., .5),
        loc=3,
        ncol=4,
        mode="expand",
        borderaxespad=0.
    )
    plt.savefig("OsteoclastShortTime.eps")
    # Long Time
    fig2, ax2 = plt.subplots()
    ax2.set_xlabel('time (days)')
    ax2.set_ylabel('OCs (u)', color='k')
    ax2.set_xlim([tas[0] - 200, tas[-1] + 200])
    ax2.set_ylim([-2, 1.05 * y1max])
    ax2.grid(True)
    ax2.plot(tas, detuas[:, 1],
             color=color1,
             linestyle='-',
             lw=2,
             marker='',
             label='osteoclast'
             )
    ax2.plot(tas, uas[:, 1],
             color=color3,
             linestyle='-.',
             lw=3,
             marker='',
             label='Sto osteoclast'
             )
    plt.legend(
        [oc_line, soc_line],
        ["OCs", "stoOCs"],
        bbox_to_anchor=(0., 1.02, 1., .5),
        loc=3,
        ncol=4,
        mode="expand",
        borderaxespad=0.
    )
    plt.savefig("OsteoclastLongTime.eps")
    plt.show()
#


def fig_ob():
    urk = np.load('SolutionDet.npy')
    ussls = np.load('SolutionSto.npy')
    detuas = np.loadtxt('SolutionDetAs.txt')
    uas = np.loadtxt('SolutionStoAs.txt')
    t = ussls[:, 0]
    tas = uas[:, 0]
    y1max = np.max([ussls[:, 2].max(), urk[:, 2].max()])
    # y2max = np.max([uas[:, 2].max(), detuas[:, 2].max()])
    # Pots Vs Time
    # color1 = colorPallete[0]
    color2 = color_pallet[1]
    # color3 = colorPallete[2]
    color4 = color_pallet[3]

    plt.close()
    w = 1.0
    fig_width_pt = 538.58 * w    # columnwidth(190mm)
    inches_per_pt = 1.0 / 72.2    # Convert pt to inch
    golden_mean = (np.sqrt(5) - 1.0) / 2.0    # Aesthetic ratio
    fig_width = fig_width_pt * inches_per_pt    # width in inches
    fig_height = fig_width * golden_mean	   # height in inches
    fig_size = [fig_width, fig_height]
    params = {
        'backend': 'ps',
        'axes.labelsize': 10,
        'text.fontsize': 10,
        'legend.fontsize': 8,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'text.usetex': True,
        'figure.figsize': fig_size
    }
    rcParams.update(params)
# ---------------------------------------------------------------
    fig1, ax1 = plt.subplots()
    # plt.title('Jerez-Chen Model Experimental Parameters')

    ob_line = mlines.Line2D([], [],
                            color=color2,
                            linestyle='-',
                            lw=2,
                            marker='',
                            markersize=1
                            )
    sob_line = mlines.Line2D([], [],
                             color=color4,
                             linestyle='-.',
                             lw=4,
                             marker='',
                             markersize=1
                             )
    ax1.set_xlabel('time (days)')
    ax1.set_ylabel('OCs(u)', color='k')
    ax1.set_ylim([-267, 1.05 * y1max])
    ax1.set_xlim([-200, 5400])
    ax1.plot(urk[:, 0], urk[:, 2],
             color=color2,
             lw=2,
             linestyle='-',
             marker='',
             label='osteoblast'
             )
    ax1.plot(t, ussls[:, 2],
             color=color2,
             linestyle='-.',
             lw=3,
             marker='',
             label='Sto osteoblast'
             )
    ax1.set_xlabel('time (days)')
    ax1.set_ylabel(r'OBs (v)')
    ax1.grid(True)
    plt.legend(
        [ob_line, sob_line],
        ["OBs", "stoOBs"],
        bbox_to_anchor=(0., 1.02, 1., .5),
        loc=3,
        ncol=4,
        mode="expand",
        borderaxespad=0.
    )
    plt.savefig("OsteoblasShortTime.eps")
#
    fig2, ax2 = plt.subplots()
    ax2.set_xlabel('time (days)')
    ax2.set_ylabel('OCs (u)', color='k')
    ax2.set_ylim([-267, 1.05 * y1max])
    ax2.set_xlim([tas[0] - 200, tas[-1] + 200])
    ax2.plot(tas, detuas[:, 2],
             color=color2,
             lw=2,
             linestyle='-',
             marker='',
             label='osteoblast'
             )
    ax2.plot(tas, uas[:, 2],
             color=color4,
             linestyle='-.',
             lw=3,
             marker='',
             label='Sto osteoblast'
             )
    ax2.set_xlabel('time (days)')
    ax2.set_ylabel(r'OBs (v)')
    ax2.grid(True)
    plt.legend([
        ob_line, sob_line],
        ["OBs", "stoOBs"],
        bbox_to_anchor=(0., 1.02, 1., .5),
        loc=3,
        ncol=4,
        mode="expand",
        borderaxespad=0.
    )
    plt.savefig("OsteoblasLongTime.eps")
    plt.show()
# #####################################################################


def fig_phase_potrait():
    urk = np.load('SolutionDet.npy')
    ussls = np.load('SolutionSto.npy')
    uas = np.loadtxt('SolutionStoAs.txt')
#    t = ussls[:, 0]
#    tas = uas[:, 0]
    # Pots Vs Time
    color1 = color_pallet[0]
    color2 = color_pallet[1]
    # color3 = colorPallete[2]
    color4 = color_pallet[3]
    plt.close()
    w = 1.0
    fig_width_pt = 538.58 * w    # \columnwidth(190mm)
    inches_per_pt = 1.0 / 72.27	    # Convert pt to inch
    golden_mean = (np.sqrt(5) - 1.0) / 2.0    # Aesthetic ratio
    fig_width = fig_width_pt * inches_per_pt    # width in inches
    fig_height = fig_width * golden_mean	     # height in inches
    fig_size = [fig_width, fig_height]
    params = {
        'backend': 'ps',
        'axes.labelsize': 10,
        'text.fontsize': 10,
        'legend.fontsize': 8,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'text.usetex': True,
        'figure.figsize': fig_size
    }
    rcParams.update(params)
    # #################################################################
    # Phase Potrait
    fig1, ax3 = plt.subplots()
    # plt.title('Stochastic Power Law Bone Remodeling')
    ax3.set_xlim(-1, 30)
    ax3.set_ylim(-200, 9500)
    ax3.plot(urk[:, 1], urk[:, 2],
             color=color2,
             linestyle='-',
             lw=3,
             label='Deterministic'
             )
    ax3.plot(StoPlbrmJC.Ubar[0], StoPlbrmJC.Ubar[1],
             'o',
             mfc='none',
             ms=8,
             # label=r'$\bar{u}$'
             )
    ax3.plot(StoPlbrmJC.Uzero[0], StoPlbrmJC.Uzero[1],
             's',
             mfc=color1,
             ms=8
             )
    ax3.plot(ussls[:, 1], ussls[:, 2],
             color=color4,
             ms=1,
             label='Stochastic Short Time'
             )
    ax3.plot(uas[:, 1], uas[:, 2],
             color=color4,
             lw=1,
             linestyle='-.',
             label='Long Time Stochastic'
             )
    # ax3.plot(uas[0,0], uas[0,1], 'x', mfc='green', ms=12)
    ax3.set_xlabel(r'$u_1$')
    ax3.set_ylabel(r'$u_2$')
    ax3.grid(True)
    ax3.legend(loc=0)
    plt.savefig("PhasePotrait.eps")
    plt.show()
# #####################################################################################################


def fig_phase_potrait3d(file_name1="PhasePotrait3d(a).eps",
                        file_name2="PhasePotrait3d(b).eps"):
    urk = np.load('SolutionDet.npy')
    ussls = np.load('OneLongPathSolutionSto.npy')
    detuas = np.load('SolutionDetAs.npy')
#    uas = np.loadtxt('SolutionStoAs.txt')
    short = 16000 * 6
    t = ussls[0:short, 0]
    tas = ussls[-short:-1, 0]
    # Pots Vs Time
    color1 = color_pallet[0]
    # color2 = colorPallete[1]
    color3 = color_pallet[2]
    det_color = color1
    sto_color = color3

    plt.close()
    w = 1.0
    fig_width_mm = 63.5 * w                      # 63.5 mm
    inches_per_mm = 1.0 / 25.4                   # Convert mm to inch
    golden_mean = (np.sqrt(5) - 1.0) / 2.0       # Aesthetic ratio
    fig_width = fig_width_mm * inches_per_mm     # width in inches
    fig_height = fig_width * golden_mean         # height in inches
    fig_size = [fig_width, fig_height]
    params = {
        'backend': 'ps',
        'axes.labelsize': 10,
        'text.fontsize': 10,
        'legend.fontsize': 8,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'text.usetex': True,
        'figure.figsize': fig_size
    }
    rcParams.update(params)
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111, projection='3d')
    ax1.frameon = False
    ax1.zaxis.set_rotate_label(False)
    ax1.set_xlabel(r'$OC(u)$')
    ax1.set_zlabel(r'$OB(v)$', rotation=90)
    #
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(True)
    ax1.spines['bottom'].set_visible(True)
    ax1.spines['left'].set_visible(True)
    yzmax = np.max(
        [ussls[0:short, 2].max(),
         urk[0:short, 2].max(),
         detuas[0:short, 2].max()]
    )
    ax1.set_zlim([-200, yzmax])
    ax1.view_init(elev=20., azim=-143)
    #
    ax1.plot(urk[0:short, 1], t, urk[0:short, 2],
             color=det_color,
             lw=1,
             linestyle='-',
             label='Det. Short Time'
             )
    ax1.plot(ussls[0:short, 1], t, ussls[0:short, 2],
             color=sto_color,
             lw=1,
             linestyle='-',
             # marker='.',
             ms=1,
             mfc='none',
             label='Sto. Short Time'
             )
    ax1.zaxis._axinfo['label']['space_factor'] = 1.75
    ax1.yaxis._axinfo['label']['space_factor'] = 2.5
    ax1.xaxis._axinfo['label']['space_factor'] = 2.75
#
    ax1.yaxis.set_ticks(np.arange(t[0], t[-1], 1500))
    ax1.set_ylabel(r'$t$')
    plt.tight_layout()
    axbox = ax1.get_position()
    x_value = .49    # Offset by eye
    y_value = .7
    ax1.legend(loc=(axbox.x0 + x_value, axbox.y0 + y_value),
               numpoints=1)
    plt.savefig(file_name1, bbox_inches='tight')
    # plt.show()
    plt.close()
    rcParams.update(params)
    #
    offset = 10
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111, projection='3d')
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(True)
    ax2.spines['bottom'].set_visible(True)
    ax2.spines['left'].set_visible(True)
    ax2.plot(detuas[0:short, 1], ussls[-1 - short:-1, 0],
             detuas[0:short, 2],
             color=det_color,
             lw=1,
             linestyle='-',
             label='Det Long Time'
             )
    ax2.plot(ussls[-short:-1:offset, 1],
             ussls[-short:-1:offset, 0],
             ussls[-short:-1:offset, 2],
             color=sto_color,
             lw=1,
             linestyle='-',
             # marker='.',
             ms=1,
             label="Sto. Long Time"
             )
    ax2.yaxis.set_ticks(np.arange(tas[0], tas[-1], 1500))
    ax2.zaxis._axinfo['label']['space_factor'] = 1.75
    ax2.yaxis._axinfo['label']['space_factor'] = 2.75
    ax2.xaxis._axinfo['label']['space_factor'] = 2.75
    ax2.set_xlabel(r'$OC(u)$')
    ax2.set_ylabel(r'$t$ (days)')
    ax2.zaxis.set_rotate_label(False)
    ax2.set_zlabel(r'$OB(v)$', rotation=90)
    ax2.zaxis.set_rotate_label(False)
    ax2.set_zlim([-200, yzmax])
    ax2.set_ylim([tas[0] - 200, tas[-1] + 200])
    ax2.view_init(elev=20., azim=-145)
    plt.tight_layout()    # pad=0.4, w_pad=0.5, h_pad=1.0)
    axbox = ax1.get_position()
    x_value = .49    # Offset by eye
    y_value = .52
    ax2.legend(loc=(axbox.x0 + x_value, axbox.y0 + y_value),
               numpoints=1)
    plt.savefig(file_name2, bbox_inches='tight')
# ---------------------------------------------------------------------


def fig_long_path(file_name1="LongPath.eps",
                  file_name2="LongPathXi.eps",
                  file_name3="PathHalfTime.eps"):
    ussls = np.load("OneLongPathSolutionSto.npy")
    urk = np.load("OneLongPathSolutionDet.npy")
    t = ussls[:, 0]
    # \Xi axe
    n = ussls[:, 0].shape[0]
    j1 = np.ones(n) * Ji1
    j2 = np.ones(n) * Ji2
    color1 = color_pallet[0]
    color2 = color_pallet[1]
    color3 = color_pallet[2]
    color4 = color_pallet[3]
    color5 = color_pallet[4]
    sto_color = color3
    det_color = color1
    ji_color = color5
    plt.close()
# EPS setup
    w = 1.0
    fig_width_pt = 120 * w                      # 120 mm
    inches_per_mm = 1.0 / 25.4                  # Convert mm to inch
    golden_mean = (np.sqrt(5) - 1.0) / 2.0      # Aesthetic ratio
    fig_width = fig_width_pt * inches_per_mm    # width in inches
    fig_height = fig_width * golden_mean        # height in inches
    fig_size = [fig_width, fig_height]
    params = {
        'backend': 'ps',
        'axes.labelsize': 10,
        'text.fontsize': 10,
        'legend.fontsize': 8,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'text.usetex': True,
        'figure.figsize': fig_size}
    rcParams.update(params)
    fig1, (ax1, ax2) = plt.subplots(nrows=2)
#
    ax1 = plt.subplot(211)
    # ax1.set_xlabel(r'$t$ (days)')
    ax1.set_ylabel(r'$OC(u)$')
    ax1.set_xlim([-200, 650 * 50])
    ax1.plot(t, ussls[:, 1],
             color=sto_color,
             lw=.5,
             linestyle='-',
             label='Deterministic'
             )
    ax1.grid(False)
    ax2 = plt.subplot(212)
    ax2.set_xlabel(r'$t$ (days)')
    ax2.set_ylabel(r'$OB(v)$')
    ax2.set_xlim([-200, 650 * 50])
    ax2.plot(t, ussls[:, 2],
             color=sto_color,
             lw=.5,
             linestyle='-',
             label='Deterministic'
             )
    ax2.grid(False)
    plt.tight_layout()
    plt.savefig(file_name1)
# ---------------------------------------------------------------
#
#  Figure likening at medium time
# ---------------------------------------------------------------
    fig2, (ax3, ax4) = plt.subplots(nrows=2)
    tm1 = 300000
    tm2 = 500000
    tml1 = t[tm1]
    tml2 = t[tm2]
    ax3 = plt.subplot(211)
    # ax3.set_xlabel(r'$t$ (days)')
    ax3.set_ylabel(r'$OC(u)$')
    ax3.set_xlim([tml1, tml2])
    ax3.plot(t[tm1:tm2], urk[tm1:tm2, 1],
             color=det_color,
             lw=.5,
             linestyle='-',
             label='Deterministic'
             )
    ax3.plot(t[tm1:tm2], ussls[tm1:tm2, 1],
             color=sto_color,
             lw=.5,
             linestyle='-',
             label='Sto. OC'
             )
    ax3.grid(False)
    #
    #
    #
    ax4 = plt.subplot(212)
    ax4.set_xlabel(r'$t$ (days)')
    ax4.set_ylabel(r'$OB(v)$')
    ax4.set_xlim([tml1, tml2])
    ax4.plot(t[tm1:tm2], urk[tm1:tm2, 2],
             color=det_color,
             lw=.5,
             linestyle='-',
             label='Deterministic'
             )
    ax4.plot(t[tm1:tm2], ussls[tm1:tm2, 2],
             color=sto_color,
             lw=.5,
             linestyle='-',
             label='Sto. OB'
             )
    ax4.grid(False)
    plt.tight_layout()
    plt.savefig(file_name3)
# ===================================================================
# Figure Xi indicator 120mm format
# ====================================================================
    fig3, (ax5, ax6) = plt.subplots(nrows=2)
    ax5 = plt.subplot(211)
    ax5.set_xlabel(r'$t$')
    ax5.set_ylabel(r'$OC(u)$')
    ax5.set_xlim([-200, 650 * 50])
    ax5.plot(
        t, ussls[:, 1],
        color=sto_color,
        lw=.5,
        linestyle='-',
        label='Sto. OCs'
    )
    ax5.plot(t, j2,
             color=ji_color,
             lw=2,
             linestyle='-',
             label=r'$\xi_2$'
             )
#
    ax5.legend(
        bbox_to_anchor=(0., 1.0, 1., .5),
        loc=3,
        ncol=2,
        mode="expand",
        borderaxespad=0.0
    )
    axins_a = inset_axes(ax5, 1, 0.5,
                         loc=2,
                         bbox_to_anchor=(0.515, 0.927),
                         bbox_transform=ax4.figure.transFigure
                         )
    y1 = Ji2 - 1.1 * Ji2
    y2 = Ji2 + Ji2 * 1.1
    t1, t2 = 25000, 30000
    axins_a.set_xlim(t1, t2)
    axins_a.set_ylim(y1, y2)
    axins_a.plot(t, j2,
                 color=ji_color,
                 lw=2,
                 linestyle='-'
                 )
    axins_a.plot(t, ussls[:, 1],
                 color=sto_color,
                 lw=1,
                 linestyle='-'
                 )
#
    plt.xticks(visible=False)
    mark_inset(
        ax5,
        axins_a,
        loc1=3,
        loc2=1,
        fc="none",
        ec="0.5"
    )
    ax6 = plt.subplot(212)
    ax6.set_xlabel(r'$t$')
    ax6.set_ylabel(r'$OB(v)$')
    ax6.set_xlim([-200, 650 * 50])
    ax6.plot(
        t, ussls[:, 2],
        color=sto_color,
        lw=.5,
        linestyle='-',
        label='Sto. OBs'
    )
    ax6.plot(t, j1,
             color=ji_color,
             lw=2,
             linestyle='-',
             label=r'$\xi_1$'
             )
    ax6.legend(
        bbox_to_anchor=(0., 1.0, 1., .5),
        loc=3,
        ncol=2,
        mode="expand",
        borderaxespad=0.0
    )
    axins_b = inset_axes(
        ax6, 1, 0.5,
        loc=3,
        bbox_to_anchor=(0.515, 0.25),
        bbox_transform=ax6.figure.transFigure
    )
    t1, t2, y1, y2 = 25000, 30000, (Ji1 - 1.1 * Ji1), (Ji1 + 1.1 * Ji1)
    axins_b.set_xlim(t1, t2)
    axins_b.set_ylim(y1, y2)
    axins_b.plot(t, j1,
                 color=ji_color,
                 lw=2,
                 linestyle='-'
                 )
    axins_b.plot(t, ussls[:, 2],
                 color=sto_color,
                 lw=1,
                 linestyle='-'
                 )
    plt.xticks(visible=False)
    mark_inset(ax6, axins_b, loc1=3, loc2=1, fc="none", ec="0.5")
    ax6.grid(False)
    plt.tight_layout()
    plt.savefig(file_name2)
    plt.close()
# ====================================================================
#
#       Second Format of Xi figures
#
# =====================================================================
    plt.close()
    rcParams.update(params)
# ======================================================================
    fig4, ax7 = plt.subplots()
    ax7.set_xlabel(r'$t$')
    ax7.set_ylabel(r'$OC(u)$')
    ax7.set_xlim([-200, 650 * 50])
    ax7.plot(
        t, ussls[:, 1],
        color=sto_color,
        lw=.5,
        linestyle='-',
        label='Sto. OCs'
    )
    ax7.plot(t, j2,
             color=ji_color,
             lw=2,
             linestyle='-',
             label=r'$\xi_2$'
             )
    ax7.legend(
        bbox_to_anchor=(0., 1.0, 1., .5),
        loc=3,
        ncol=2,
        mode="expand",
        borderaxespad=0.0
    )
    axins = inset_axes(ax7, 1, .5,
                       loc=2,
                       bbox_to_anchor=(0.35, 0.55),
                       bbox_transform=ax4.figure.transFigure
                       )
    y1 = Ji2 - 1.1 * Ji2
    y2 = Ji2 + Ji2 * 1.1
    t1, t2 = 25000, 30000
    axins.set_xlim(t1, t2)
    axins.set_ylim(y1, y2)
    axins.plot(t, j2,
               color=ji_color,
               lw=2,
               linestyle='-'
               )
    axins.plot(t, ussls[:, 1],
               color=sto_color,
               lw=1,
               linestyle='-'
               )
    plt.xticks(visible=False)
    mark_inset(
        ax7,
        axins,
        loc1=3,
        loc2=1,
        # fc="none",
        ec="0.5"
    )
    plt.tight_layout()
    plt.savefig("LongPath(OCs).eps")
# ==============================================================
    plt.close()
    rcParams.update(params)
    fig5, ax8 = plt.subplots()
    ax8.set_xlabel(r'$t$')
    ax8.set_ylabel(r'$OB(v)$')
    ax8.set_xlim([-200, 650 * 50])
    ax8.plot(t, j1,
             color=ji_color,
             lw=2,
             linestyle='-',
             label=r'$\xi_1$'
             )
    ax8.plot(
        t, ussls[:, 2],
        color=color1,
        lw=.5,
        linestyle='-',
        label='Sto. OBs'
    )
    ax8.legend(
        bbox_to_anchor=(0., 1.0, 1., .5),
        loc=3,
        ncol=2,
        mode="expand",
        borderaxespad=0.0
    )
    axins = inset_axes(
        ax8, 1, .5,
        loc=2,
        bbox_to_anchor=(0.35, .55),
        bbox_transform=ax6.figure.transFigure
    )
    # axins = zoomed_inset_axes(ax5, 2.5, loc=2)
    t1, t2, y1, y2 = 25000, 30000, (Ji1 - 1.1 * Ji1), (Ji1 + 1.1 * Ji1)
    axins.set_xlim(t1, t2)
    axins.set_ylim(y1, y2)
    axins.plot(t, j1,
               color=ji_color,
               lw=2,
               linestyle='-'
               )
    axins.plot(t, ussls[:, 2],
               color=sto_color,
               lw=1,
               linestyle='-'
               )
    plt.xticks(visible=False)
    mark_inset(ax8, axins, loc1=3, loc2=1, fc="none", ec="0.5")
    ax7.legend(loc='upper center',
               bbox_to_anchor=(0.5, 1.05),
               ncol=3,
               fancybox=True,
               shadow=True)
    plt.tight_layout()
    plt.savefig("LongPath(OBs).eps")
#    ================================================================


def fig_bone_mass():
    z = np.load('BoneMass.npy')
    t = z[:, 0]
    z_det = z[:, 1]
    z_sto = z[:, 2]
    fig_width_pt = 256    # 538.58*w	 \columnwidth(140mm)
    inches_per_pt = 1.0 / 72.27    # Convert pt to inch
    golden_mean = (np.sqrt(5) - 1.0) / 2.0    # Aesthetic ratio
    fig_width = fig_width_pt * inches_per_pt    # width in inches
    fig_height = fig_width * golden_mean    # height in inches
    fig_size = [fig_width, fig_height]
    params = {
        'backend': 'ps',
        'axes.labelsize': 10,
        'text.fontsize': 10,
        'legend.fontsize': 8,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'text.usetex': True,
        'figure.figsize': fig_size
    }
    rcParams.update(params)
    color1 = color_pallet[0]
#    color2 = colorPallete[1]
    color3 = color_pallet[2]
    color4 = color_pallet[3]

    plt.close()
    rcParams.update(params)
    plt.figure()
    plt.xlabel(r'$t$ (days)')
    plt.ylabel(r'Bone mass (z)')
    plt.xlim([0, 5000])
    plt.plot(t, 100 * np.ones(t.shape[0]),
             color=color4,
             lw=.5,
             linestyle=':'
             )
    plt.plot(t, 100 * z_det,
             color=color1,
             lw=1,
             linestyle='-',
             label='Deterministic'
             )
    plt.plot(t, 100 * z_sto,
             color=color3,
             lw=1,
             linestyle='-',
             label='Stochastic'
             )
    plt.legend(loc=0)
    plt.tight_layout()
    plt.savefig("BoneMass.eps")
#


def fig_short_time_moments():
    mean_short = np.load('MeanUShortTime.npy')
    ms_short = np.load('MSUShortTime.npy')
    u1_mean = mean_short[:, 0]
    u2_mean = mean_short[:, 1]
    u1_ms = ms_short[:, 0]
    u2_ms = ms_short[:, 1]
    n = u1_mean.shape[0]
    t = StoPlbrmJC.Dt * 100 * np.arange(n)
    #  Pots Vs Time
    color1 = color_pallet[0]
    color2 = color_pallet[1]
    color3 = color_pallet[2]
    color4 = color_pallet[3]
    plt.close()
    w = 1.0
    fig_width_pt = 120 * w    # \columnwidth(190mm)
    inches_per_pt = 1.0 / 25.4    # Convert pt to inch
    golden_mean = (np.sqrt(5) - 1.0) / 2.0    # Aesthetic ratio
    fig_width = fig_width_pt * inches_per_pt    # width in inches
    fig_height = fig_width * golden_mean    # eight in inches
    fig_size = [fig_width, fig_height]
    params = {'backend': 'ps',
              'axes.labelsize': 10,
              'text.fontsize': 10,
              'legend.fontsize': 8,
              'xtick.labelsize': 8,
              'ytick.labelsize': 8,
              'text.usetex': True,
              'figure.figsize': fig_size}
    rcParams.update(params)
    fig = plt.figure()
    ax1 = plt.subplot(221)
    ax1.set_xlabel(r'$t$')
    ax1.set_ylabel(r'$OC(u)$')
    ax1.set_xlim([0, 5200])
    ax1.plot(t, u1_mean,
             color=color2,
             lw=1.0,
             linestyle='-',
             label=r'Mean Short'
             )
    ax1.plot(t, u1_ms,
             color=color4,
             lw=1.0,
             linestyle='-',
             label='MS-Short'
             )
    ax1.legend(loc=0)
#
    ax2 = plt.subplot(223)
    ax2.set_xlabel(r'$t$')
    ax2.set_ylabel(r'$OB(v)$')
    ax2.set_xlim([0, 5200])
    ax2.plot(t, u2_mean,
             color=color2,
             lw=1.0,
             linestyle='-',
             label=r'Mean Short'
             )
#
    ax2.plot(t, u2_ms,
             color=color4,
             lw=1,
             linestyle='-',
             label='MS-Short'
             )
    # ax2.legend(loc=0)
#
    ax3 = plt.subplot(122)
    ax3.set_xlabel(r'$OC(u)$')
    ax3.set_ylabel(r'$OB(v)$')
    ax3.plot(u1_mean, u2_mean,
             color=color2,
             lw=1.0,
             linestyle='-',
             label='Mean Short'
             )
    ax3.plot(u1_ms, u2_ms,
             color=color4,
             lw=1.0,
             linestyle='-',
             marker='.',
             ms=1,
             # mfc="none",
             label='MS-Short'
             )
    ax3.legend(loc=0)
    ax3.grid(True)
#
    plt.tight_layout()
    plt.savefig("ShortPathMoments.eps")
    plt.show()
#


def fig_long_time_moments():
    mean_long = np.load('MeanULongTime.npy')
    ms_long = np.load('MSULongTime.npy')
    u1_mean = mean_long[:, 0]
    u2_mean = mean_long[:, 1]
    u1_ms = ms_long[:, 0]
    u2_ms = ms_long[:, 1]
    n = u1_mean.shape[0]
    t = 31200 + StoPlbrmJC.Dt * 100 * np.arange(n)
    # Pots Vs Time
    plt.close()
    # color1 = color_pallet[0]
    color2 = color_pallet[1]
    # color3 = color_pallet[2]
    color4 = color_pallet[3]

    w = 1.0
    fig_width_pt = 120 * w    # \columnwidth(190mm)
    inches_per_pt = 1.0 / 25.4    # Convert pt to inch
    golden_mean = (np.sqrt(5) - 1.0) / 2.0     # Aesthetic ratio
    fig_width = fig_width_pt * inches_per_pt     # width in inches
    fig_height = fig_width * golden_mean    # height in inches
    fig_size = [fig_width, fig_height]
    params = {'backend': 'ps',
              'axes.labelsize': 10,
              'text.fontsize': 10,
              'legend.fontsize': 8,
              'xtick.labelsize': 8,
              'ytick.labelsize': 8,
              'text.usetex': True,
              'figure.figsize': fig_size}
    rcParams.update(params)

#    fig = plt.figure()
    ax1 = plt.subplot(221)
    ax1.set_xlabel(r'$t$')
    ax1.set_ylabel(r'$OC(u)$')
    ax1.set_xlim([31200, 36400])
    ax1.plot(t, u1_mean,
             color=color2,
             lw=1.0,
             linestyle='-',
             label='Mean-Long'
             )
    ax1.plot(t, u1_ms,
             color=color4,
             lw=1.0,
             linestyle='-',
             label='MS-Long'
             )
    ax1.legend(loc=0)
#
    ax2 = plt.subplot(223)
    ax2.set_xlabel(r'$t$')
    ax2.set_ylabel(r'$OB(v)$')
    ax2.set_xlim([31200, 36400])
    ax2.plot(t, u2_mean,
             color=color2,
             lw=1.0,
             linestyle='-',
             label='Mean-Long'
             )
    ax2.plot(t, u2_ms,
             color=color4,
             lw=1.0,
             linestyle='-',
             label='MS-Long'
             )
    # ax2.legend(loc=0)
#
    ax3 = plt.subplot(122)
    ax3.set_xlabel(r'$OC(u)$')
    ax3.set_ylabel(r'$OB(v)$')
    # ax3.set_xlim([-5,40])
    # ax3.set_ylim([-200, 12700])
    ax3.plot(u1_mean, u2_mean,
             color=color2,
             lw=1.0,
             linestyle='-',
             label='Mean-Long'
             )
    ax3.plot(u1_ms, u2_ms,
             color=color4,
             lw=1.0,
             linestyle='-',
             marker='.',
             ms=1,
             # mfc="none",
             label='MS-Long'
             )
    ax3.legend(loc=0)
    ax3.grid(True)
#
    plt.tight_layout()
    plt.savefig("LongPathMoments.eps")

# Model parameters ---------------------------------------------------
a1 = 0.3        # Model parameters  ai = \alpha_i
a2 = 0.18
b1 = 0.2
b2 = 0.02
ns = 0.1        # Noise intensity
gamma1 = -0.8
gamma2 = 0.5
sigma = np.array([ns * b1, ns * b2])
k1 = 0.03
k2 = 0.0017
#
Ji1 = ((b1 + 0.5 * (sigma[0] ** 2)) / a1) ** (1.0 / gamma1)
Ji2 = ((b2 + 0.5 * (sigma[1] ** 2)) / a2) ** (1.0 / gamma2)
# Stencil Parameters
U0 = [10, .7]
k = 5
p = 0
r = p
T0 = 0.0
# T = 1300
T = 650 * 8
LTM = 5
#
StoPlbrmJC = StoPLBRM()
StoPlbrmJC.InitializeMesh(k, p, r, T0, T)
StoPlbrmJC.NoiseUpdate(123456789)
StoPlbrmJC.SetParametersStoPLBRM(a1, b1, a2, b2, 1.0,
                                 gamma2, gamma1, 1.0, k1, k2, sigma, U0)
# color Pallet
# colorPallete = ['#f0f9e8', '#bae4bc', '#7bccc4', '#2b8cbe']
# colorPallete = ['#F6372D', '#202529', '#F6A178', '#C6E0F1']
color_pallet = ['#588C7E', '#F2E394', '#F2AE72', '#D96459', '#8C4646']
fig_long_path(file_name1="Fig2.eps", file_name2="Fig3.eps",
              file_name3="Fig4.eps")
fig_phase_potrait3d(file_name1="Fig5.eps", file_name2="Fig6.eps")

