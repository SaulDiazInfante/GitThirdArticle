import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from expStoPLBR2Noises import StochasticPLBRM
from DeterministicModelKomarova import PLBRM
a1 = .30
a2 = .10
b1 = 0.2
b2 = 0.02
g11 = 1.0
g12 = 0.5
g21 = -0.3
g22 = 1.0
k1 = 0.03
k2 = 0.0017
#U0=[10.0+2.3094010767585034,3.0517578125]
U0=[.4, 3.85]
sigma=np.array([b1/500.0,b2/5000.0])
#sigma=np.array([0.0, 0.0 ])
#Stencil parameters
k=8.0
r=1.0
T0=0.0
T=2000.0
#Deterministic Part
plbrm=PLBRM()
plbrm.InitializeMesh(k-3,T0,T)
plbrm.SetParametersKomarova(a1, b1, a2, b2, g11, g12, g21, g22, k1, k2, U0)
plbrm.OdeSolve()
plbrm.SaveData()
#Stochastic Part
StoPlbrm=StochasticPLBRM()
StoPlbrm.InitializeMesh(k, r, T0, T)
StoPlbrm.SetParametersStoPLBRM(a1, b1, a2, b2, g11, g12, g21, g22, k1, k2, sigma, U0)
StoPlbrm.SodeSolver()
StoPlbrm.SaveData()
