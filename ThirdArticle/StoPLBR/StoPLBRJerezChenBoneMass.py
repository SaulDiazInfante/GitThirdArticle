import numpy as np
import matplotlib.pyplot as plt
from cPickle import dump
from cPickle import load
from StoPLBRJerezChen import StoPLBRM

# Model parameters
#  k1 = 0.093,  k2 = .00032, r=290.625
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
# Stencil Parameters
k = 5
p = 0
r = p
T0 = 0.0
offset = 48
T =650 * offset  # almost 48 periods
"""
with open('state.obj', 'wb') as f:
    dump(np.random.get_state(), f)


"""
with open('state.obj', 'rb') as f:
    np.random.set_state(load(f))


bm = StoPLBRM()
bm.__init__(k=5, p=0, r=0, t_0=0, t_f=T)
# bm.load_parameters()
bm.initialize_mesh(k, p, r, T0, T)
# bm.noise_update()
gamma1 = bm.g21
gamma2 = bm.g12
Ji1 = bm.ji_1
Ji2 = bm.ji_2
print '\n\n'
condition_H1 = (np.sign(gamma1) < 0.0) and (np.sign(gamma2) > 0.0)
condition_H2 = (np.abs(gamma1) <= np.sign(gamma2))
condition_H3 = bm.a1 * gamma2 < bm.a2 * np.abs(gamma1)
print'\t\t Xi1, Xi2:= '+'['+str(Ji1)+', '+str(Ji2)+']'
if condition_H1 and (condition_H2 and condition_H3):
    print '\n\t\t\tH3 holds =) \n\n'
else:
    print '\n\t\t\tH3 not holds =(\n\n'
    print 'p1:=', bm.a1 * gamma2, \
          'p2:=', bm.a2 * np.abs(gamma1)



# u_det, z = bm.ssls_det()
# u_det, z = bm.rk()
u_det = np.load("bone_mass_xppaut_graph.npy")
t = u_det[:, 0]
z = u_det[:, 3]
iter_max = 1000
bone_condition = True
j = 1
eps = 0.01
k1 = bm.k1
k2 = bm.k2
while bone_condition and j < iter_max:
    bm.k1_sto = k1 + .25 * k1 * np.random.randn()
    bm.k2_sto = k2 + .25 * k2 * np.random.randn()
    print '{:d}  k1_sto:={:1.8f} k2_sto={:1.8f}'.\
        format(j, bm.k1_sto, bm.k2_sto)
    print '\n'
    # print str(j) + '\t' + str(bm.k1_sto) + '\t' + str(bm.k2_sto)
    #  print '\n'
    ussls, z_sto = bm.ssls()
    bone_condition = not bm.bone_flag
    plt.close()
    plt.plot(t, z)
    plt.plot(bm.t_k, z_sto)
    plt.savefig("sto_bonemass.png")

    plt.close()
    plt.plot(bm.t_k, ussls[:, 1])
    plt.savefig("ob_path_guide.png")
    j += 1
    if bone_condition:
        del ussls
        del z_sto
    j += 1

# np.save("OneLongPathSolutionSto.npy", np.transpose(np.array([t,
# bm.u_ssls[:, 0], bm.u_ssls[:, 1]])))
# ones = np.ones(t.shape[0])
# z = bm.bone_mass()
# np.save('BoneMass.npy',
#  np.transpose(np.array([t[0:-1:100], bm.z[0:-1:100]])))
