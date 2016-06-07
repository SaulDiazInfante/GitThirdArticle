import numpy as np
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
from DeterministicModelKomarova import PLBRM
a1=3.0
a2=4.0
b1=0.2
b2=0.02
g11=1.1
g12=1.0
g21=-0.5
g22=0.0
k1=0.03
k2=0.0017
#U0=[10.0+2.3094010767585034,3.0517578125]
U0=[1.0, 0.25]
#Stencil Parameters
k=5.0
T0=0.0
T=1000
plbrm=PLBRM()
plbrm.InitializeMesh(k,T0,T)
plbrm.SetParametersKomarova(a1,b1,a2,b2,g11,g12,g21,g22,k1,k2,U0)
plbrm.OdeSolve()
#
fig1=plt.figure()
ax1 = fig1.add_subplot(121)
ax1.plot(plbrm.t,plbrm.U1,'k-',label='osteoclast')
ax1.plot(plbrm.t,plbrm.Ubar[0]*np.ones(plbrm.U1.shape[0]),'r--',label='Steady State')
ax1.set_xlabel(r't')
ax1.set_ylabel(r'Osteoclast')
ax1.grid(True)
ax1.legend(loc=0)
#
ax2 = fig1.add_subplot(122)
ax2.plot(plbrm.t,plbrm.U2,'k-',label='osteoblast')
ax2.plot(plbrm.t,plbrm.Ubar[1]*np.ones(plbrm.U1.shape[0]),'r--',label='Steady State')
ax2.set_xlabel(r't')
ax2.set_ylabel(r'Osteoblast')
ax2.grid(True)
ax2.legend(loc=0)
plt.show()
#
fig2=plt.figure()
ax1 = fig2.add_subplot(111)
ax1.plot(plbrm.Ubar[0],plbrm.Ubar[1],'o', mfc='red',ms=8)
ax1.plot(plbrm.Uzero[0],plbrm.Uzero[1],'s', mfc='green',ms=8)
ax1.plot(plbrm.U1,plbrm.U2,'k-')
ax1.set_xlabel(r'u_1')
ax1.set_ylabel(r'u_2')
ax1.grid(True)
ax1.legend(loc=0)
plt.show()
#
fig3, ax1 = plt.subplots()
OC_line = mlines.Line2D([], [], color='black', linestyle='-', marker='', markersize=1)
OB_line = mlines.Line2D([], [], color='black', linestyle=':', marker='', markersize=1)
ax1.plot(plbrm.t,plbrm.U1,'k-',label='osteoclast')
ax1.set_xlabel('time (days)')
ax1.set_ylabel('OCs (u1)', color='b')
ax2 = ax1.twinx()
ax2.plot(plbrm.t,plbrm.U2,'k:',label='osteoblast')
ax2.set_xlabel('time (days)')
ax2.set_ylabel(r'OBs (u2)')
ax2.grid(True)
plt.legend([OC_line, OB_line], ["OCs", "OBs"], loc=0)
plt.show()