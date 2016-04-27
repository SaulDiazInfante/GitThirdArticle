import numpy as np
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
from DeterministicModelJerezChen import PLBRMJerezChen
# Model parameters -----------------------------------------------------------------------------------------------------
# a1 = 3.0
# a2 = 4.0
# b1 = 0.2
# b2 = 0.02
# g11 = 1.1
# g12 = 0.5
# g21 = -0.5
# g22 = 0.0
# k1 = 0.03
# k2 = 0.0017
# #U0=[10.0+2.3094010767585034,3.0517578125]
# U0=[10.0, 1.0]
# #Stencil Parameters
# k=5.0
# T0=0.0
# T=1000
# plbrm=PLBRMJerezChen()
# plbrm.InitializeMesh(k,T0,T)
# plbrm.SetParametersKomarova(a1,b1,a2,b2,g11,g12,g21,g22,k1,k2,U0)
# plbrm.OdeSolve()
# plt.close()
# #
# fig1=plt.figure()
# plt.title('Komarova model Ayati parameters')
# ax1 = fig1.add_subplot(121)
# ax1.plot(plbrm.t,plbrm.U1,'k-',label='osteoclast')
# ax1.plot(plbrm.t,plbrm.Ubar[0]*np.ones(plbrm.U1.shape[0]),'r--',label='Steady State')
# ax1.set_xlabel(r't')
# ax1.set_ylabel(r'Osteoclast')
# ax1.grid(True)
# ax1.legend(loc=0)
# #
# ax2 = fig1.add_subplot(122)
# ax2.plot(plbrm.t,plbrm.U2,'k-',label='osteoblast')
# ax2.plot(plbrm.t,plbrm.Ubar[1]*np.ones(plbrm.U1.shape[0]),'r--',label='Steady State')
# ax2.set_xlabel(r't')
# ax2.set_ylabel(r'Osteoblast')
# ax2.grid(True)
# ax2.legend(loc=0)
# plt.show()
# #
# fig2=plt.figure()
# plt.title('Komarova model Ayati parameters')
# ax1 = fig2.add_subplot(111)
# ax1.plot(plbrm.Ubar[0], plbrm.Ubar[1],'o', mfc='red',ms=8)
# ax1.plot(plbrm.Uzero[0], plbrm.Uzero[1],'s', mfc='green',ms=8)
# ax1.plot(plbrm.U1, plbrm.U2,'k-')
# ax1.set_xlabel(r'$u_1$')
# ax1.set_ylabel(r'$u_2$')
# ax1.grid(True)
# ax1.legend(loc=0)
# plt.show()
# #
# fig3, ax1 = plt.subplots()
# plt.title('Komarova model Ayati parameters')
# OC_line = mlines.Line2D([], [], color='black', linestyle='-', marker='', markersize=1)
# OB_line = mlines.Line2D([], [], color='black', linestyle=':', marker='', markersize=1)
# ax1.plot(plbrm.t,plbrm.U1,'k-',label='osteoclast')
# ax1.set_xlabel('time (days)')
# ax1.set_ylabel('OCs (u1)', color='b')
# ax2 = ax1.twinx()
# ax2.plot(plbrm.t,plbrm.U2,'k:',label='osteoblast')
# ax2.set_xlabel('time (days)')
# ax2.set_ylabel(r'OBs (u2)')
# ax2.grid(True)
# plt.legend([OC_line, OB_line], ["OCs", "OBs"], loc=0)
# plt.show()
#=======================================================================================================================
#=======================================================================================================================
# Stability analysis of a Komarova type model for the interactions of osteo blast and 
# osteo clast cells during bone remodeling
# Jerez-Chen 
#=======================================================================================================================
# Model parameters -----------------------------------------------------------------------------------------------------
a1 = 0.3
a2 = 0.1
b1 = 0.2
b2 = 0.02
gamma1 = -0.3
gamma2 = 0.5
k1 = 0.03
k2 = 0.0017
#Stencil Parameters
U0=[10, 1.0]
k = 6.0
T0 = 0.0
T = 2000
plbrmJC = PLBRMJerezChen()
plbrmJC.InitializeMesh(k,T0,T)
plbrmJC.SetParametersJerezChen(a1, b1, a2, b2, gamma1, gamma2, k1, k2, U0)
plbrmJC.OdeSolve()
Ussls = plbrmJC.SSLS()
#
fig4 = plt.figure()
plt.title('Jerez-Chen parameters')
ax1 = fig4.add_subplot(121)
ax1.plot(plbrmJC.t, plbrmJC.U1,'k-', label='osteoclast')
ax1.plot(plbrmJC.t, plbrmJC.Ubar[0]*np.ones(plbrmJC.U1.shape[0]), 'r--', label='Steady State')
ax1.set_xlabel(r't')
ax1.set_ylabel(r'Osteoclast')
ax1.grid(True)
ax1.legend(loc=0)
#
ax2 = fig4.add_subplot(122)
ax2.plot(plbrmJC.t, plbrmJC.U2,'k-',label='osteoblast')
ax2.plot(plbrmJC.t, plbrmJC.Ubar[1]*np.ones(plbrmJC.U1.shape[0]),'r--',label='Steady State')
ax2.set_xlabel(r't')
ax2.set_ylabel(r'Osteoblast')
ax2.grid(True)
ax2.legend(loc=0)
plt.show()
#
fig5 = plt.figure()
plt.title('Jerez-Chen parameters')
ax1 = fig5.add_subplot(111)
ax1.plot(plbrmJC.Ubar[0], plbrmJC.Ubar[1],'o', mfc='red',ms=8)
ax1.plot(plbrmJC.Uzero[0], plbrmJC.Uzero[1],'s', mfc='green',ms=8)
ax1.plot(plbrmJC.U1, plbrmJC.U2,'k-')
ax1.plot(Ussls[:,0], Ussls[:,1],'r--')
ax1.set_xlabel(r'u_1')
ax1.set_ylabel(r'u_2')
ax1.grid(True)
ax1.legend(loc=0)
plt.show()
#
fig6, ax1 = plt.subplots()
plt.title('Jerez-Chen parameters')
OC_line = mlines.Line2D([], [], color='black', linestyle='-', marker='', markersize=1)
OB_line = mlines.Line2D([], [], color='black', linestyle=':', marker='', markersize=1)
ax1.plot(plbrmJC.t, plbrmJC.U1,'k-',label='osteoclast')
ax1.set_xlabel('time (days)')
ax1.set_ylabel('OCs (u1)', color='k')
ax2 = ax1.twinx()
ax2.plot(plbrmJC.t, plbrmJC.U2,'k:',label='osteoblast')
ax2.set_xlabel('time (days)')
ax2.set_ylabel(r'OBs (u2)')
ax2.grid(True)
plt.legend([OC_line, OB_line], ["OCs", "OBs"], loc=0)
plt.show()
