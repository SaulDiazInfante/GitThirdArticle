import numpy as np
import matplotlib.pyplot as plt
from StoPLBRJerezChen import StoPLBRM

# Model parameters ----------------------------------------------------
a1 = 0  # Model's parameters
a2 = 0.18
b1 = 0.2
b2 = 0.02

ns = 0.1 # Noise Amplitude
gamma1 = -0.9
gamma2 = 0.5
#
sigma = np.array([ns*b1, ns*b2])
r=603
rSto = 745
k2 = 0.0000038
k1 =  r*k2
U0=[10, 0.7]    # Initial Condition
# Stencil Parameters
k = 5
p = 0
r = p
T0 = 0.0
T = 8*650
Ji1 = (( b1 + 0.5 * (sigma[0] ** 2) ) / a1) ** (1.0 / gamma1)
Ji2 = (( b2 + 0.5 * (sigma[1] ** 2) ) / a2) ** (1.0 / gamma2)
print 'Xi1, Xi2:= '+'['+str(Ji1)+', '+str(Ji2)+']'
if a1 * gamma2 < a2 * np.abs(gamma1): 
	print 'H3 holds =)' 
else:
	print 'H3 not holds =(' 
#
StoPlbrmJC = StoPLBRM()
StoPlbrmJC.InitializeMesh(k, p, r, T0, T)
#StoPlbrmJC.NoiseUpdate(123456789)
StoPlbrmJC.SetParametersStoPLBRM(a1, b1, a2, b2, 1.0, gamma2, gamma1, 1.0, k1, k2, sigma, U0)
t = StoPlbrmJC.t[0 : -1 : StoPlbrmJC.R].reshape([StoPlbrmJC.t[0 : - 1 : StoPlbrmJC.R].shape[0],1])
#Urk=StoPlbrmJC.RK()
seed = np.random.random_integers(1, 123456789)
np.random.seed(123456789)
StoPlbrmJC.NoiseUpdate(seed, flag=1)
Ustk=StoPlbrmJC.SSLS(123456789, U0, fn=0.0)

ussls=StoPlbrmJC.SSLS(seed, U0, fn=1.0)
StoPlbrmJC.SaveData()

StoPlbrmJC.BoneMass(r, rSto,.0000039)
np.save('BoneMass.npy', np.transpose(np.array([t[:,0], StoPlbrmJC.z , StoPlbrmJC.zSto])))
#'''
plt.close()
plt.figure()
plt.plot(t, ussls[:,0])
plt.plot(t, Ustk[:, 0])
plt.ylabel('OC')
#
plt.figure()
plt.plot(t, ussls[:,1])
plt.plot(t, Ustk[:, 1])
plt.ylabel('OB')
#'''
plt.figure()
plt.plot(t, 100 * np.ones(StoPlbrmJC.z.shape[0]))
plt.plot(t, 100 * StoPlbrmJC.z)
plt.plot(t, 100 * StoPlbrmJC.zSto)
plt.show()