import numpy as np
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from DeterministicModelKomarova import PLBRM
class PLBRMJerezChen(PLBRM):
	def SetParametersJerezChen(self, a1, b1, a2, b2, gamma1, gamma2, k1, k2, U0):
		self.a1 = a1
		self.b1 = b1
		self.a2 = a2
		self.b2 = b2
		self.g11 = 1.0
		self.g12 = gamma2
		self.g21 = gamma1
		self.g22 = 1.0
		self.k1 = k1
		self.k2 = k2
		#steady states
		g12 = self.g12
		g21 = self.g21
		g11 = self.g11
		g22 = self.g22
		self.gamma = g12*g21 - (1.0 - g11) * (1.0 - g22)
		u1bar = ((b1 / a1) ** ((1.0 - g22) / self.gamma)) * (b2 / a2) ** ( g21 / self.gamma)
		u2bar = ((b1 / a1) ** (g12 / self.gamma)) * ((b2 / a2) **( (1.0 - g11)/self.gamma))
		self.Ubar=np.array([u1bar,u2bar])
		self.Uzero=np.array([U0[0],U0[1]])
#	
	def SSLS(self):
		alpha1 = self.a1
		alpha2 = self.a2
		beta1 = self.b1
		beta2 = self.b2
		gamma1 = self.g21
		gamma2 = self.g12
		h = self.dt
		size = self.N
		self.Ussls = np.zeros([size+1, 2])
		self.Ussls[0, :] = self.Uzero
		for j in np.arange(size):
			uj1 = self.Ussls[j, 0]
			uj2 = self.Ussls[j, 1]
			a11 = (alpha1 * np.exp(gamma1* np.log(uj2)) - beta1) 
			Uj1 = uj1 * np.exp(h * a11)
			a12 = (alpha2 * np.exp(gamma2* np.log(Uj1)) - beta2)
			Uj2 = uj2 * np.exp(h * a12)
			self.Ussls[j+1, 0] = Uj1
			self.Ussls[j+1, 1] = Uj2
		Usol = self.Ussls
		return Usol