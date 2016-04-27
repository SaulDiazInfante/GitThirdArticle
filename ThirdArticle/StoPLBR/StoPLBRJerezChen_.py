import numpy as np
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
from DeterministicModelJerezChen import PLBRMJerezChen
#from sphinx.ext.todo import Todo
from scipy.optimize.minpack import fsolve
class StoPLBRM(PLBRMJerezChen):
	def InitializeMesh(self, k, p, r, T0, T):
		self.k = k
		self.p = p
		self.r = r
		self.T0 = T0
		self.N = 10.0**k
		self.P = 10.0**p
		self.R = 10.0**r
		#self.N = 2.0 ** k
		#self.P = 2.0 ** p
		#self.R = 2.0 ** r
		self.T = T
		self.dt = self.T / np.float(self.N)
		self.IndexN = np.arange(self.N + 1)
		# set of index to Ito integral ------------------------------------------------------
		self.tau = self.IndexN[0:self. N + 1:self.P]
		self.t = np.linspace(0, self.T, self.N + 1)
		self.Dt = np.float(self.R) * self.dt
		self.L = self.N / self.R
#
	def NoiseUpdate(self, flag=1):
		self.dWj = np.random.randn(2, np.int(self.N + 1))
		self.dWj = np.sqrt(self.dt) * self.dWj
		#self.dWj = np.cumsum(self.dWj,axis=1)
		if flag==1:
			self.dWj[:,0] = 0.0		
#
	def SetParametersStoPLBRM(self, a1, b1, a2, b2, g11, g12, g21, g22, k1, k2, sigma, U0):
		self.a1=a1
		self.b1=b1
		self.a2=a2
		self.b2=b2
		self.g11=g11
		self.g12=g12
		self.g21=g21
		self.g22=g22
		self.k1=k1
		self.k2=k2
		self.sigma=sigma
		#steady states
		self.gamma=g12*g21-(1.0-g11)*(1.0-g22)
		u1bar=(b1/a1)**((1.0-g22)/self.gamma)*(b2/a2)**(g21/self.gamma)
		u2bar=((b1/a1)**(g12/self.gamma))*((b2/a2)**((1.0-g11)/self.gamma))
		self.Ubar=np.array([u1bar,u2bar])
		self.Uzero=np.array([U0[0], U0[1]])
#
	def a(self, u):
		#===================================================================================
		# Diffusion
		#===================================================================================
		alpha1 = self.a1
		alpha2 = self.a2
		beta1 = self.b1
		beta2 = self.b2
		gamma1 = self.g21
		gamma2 = self.g12
		u1 = u[0]
		u2 = u[1]
		U1 = u1 * (alpha1 * np.exp(gamma1 * np.log(u2)) - beta1)
		U2 = u2 * (alpha2 * np.exp(gamma2 * np.log(u1)) - beta2)
		return np.array([U1, U2]).reshape([2,1])
	def b(self, U):
		sigma1 = self.sigma[0]
		sigma2 = self.sigma[1]
		u1 = U[0,0]
		u2 = U[1,0] 
		return np.diag([sigma1 * u1, sigma2 * u2])
#
	def ai(self, X, j):
		h = self.Dt
		DW = self.Winc
		X0 = self.Ubem[j].reshape([2, 1])
		x0  = X0+ np.dot(self.b(X0),DW)
		fx = self.a(X)
		return (h * fx.ravel() + x0.ravel() - X.ravel())
#
	def EM(self, flag=0):
		h = self.Dt
		L = self.L
		R = self.R
		if flag == 1:
			Dt = self.dt
			L = self. N
			R =1
		self.Uem = np.zeros([L, 2])
		self.Uem[0] = self.Uzero
		for j in np.arange(L - 1):
			self.Winc = np.sum(self.dWj[:,R * (j):R * (j + 1)],axis=1)
			self.Winc = self.Winc.reshape([2,1])
			uj = self.Uem[j,:].reshape([2,1])
			aj = self.a(uj)
			increment = uj + h *aj + np.dot(self.b(uj),self.Winc)
			self.Uem[j+1, :] = increment[:, 0]
		uem = self.Uem
		return uem
#
	def TamedEuler(self, flag=0):
		h = self.Dt
		L = self.L
		R = self.R
		if flag == 1:
			Dt = self.dt
			L = self. N
			R =1
		self.Utem = np.zeros([L, 2])
		self.Utem[0] = self.Uzero
		for j in np.arange(L - 1):
			self.Winc = np.sum(self.dWj[:,R * (j):R * (j + 1)],axis=1)
			self.Winc = self.Winc.reshape([2,1])
			uj = self.Utem[j,:].reshape([2,1])
			aj = h*self.a(uj)
			increment = uj + aj/(1+np.linalg.norm(aj)) + np.dot(self.b(uj),self.Winc)
			self.Utem[j+1, :] = increment[:, 0]
			utem = self.Utem
		return utem
#
	def RK(self, Uzero = [1, 1]):
		h = self.Dt
		Dt =self.Dt
		L = self.L
		R = self.R
		self.Urk = np.zeros([L, 2])
		self.Urk[0] = self.Uzero
		if Uzero[0]!=1:	
			self.Urk[0] = Uzero
#
		for j in np.arange(self.L - 1):
			uj = self.Urk[j].reshape([2,1])
			K1 = h * self.a(uj)
			K2 = h * self.a(uj + 0.5 *K1)
			K3 = h * self.a(uj + K2/2.0)
			K4 = h * self.a(uj + K3)
			increment = 1/6.0 * (K1 + 2.0 * K2 + 2.0 * K3 + K4)
			self.Urk[j+1, :] = uj[:, 0] + increment[:, 0]
		urk = self.Urk
		return urk
#
	def BIM(self):
		h = self.Dt
		Dt =self.Dt
		L = self.L
		R = self.R
		b = np.diag([self.sigma[0], self.sigma[1]])
		self.Ubim = np.zeros([L, 2])
		self.Ubim[0] = self.Uzero
		# Find parameters  of implicitness ---------------------------------
#
	def BEM(self):
		L = self.L
		R = self.R
		self.Ubem = np.zeros([L, 2])
		self.Ubem[0] = self.Uzero
		for j in np.arange(L-1):
			self.Winc = np.sum(self.dWj[:, R * (j) : R * (j + 1)], axis=1)
			self.Winc = self.Winc.reshape([2, 1])
			uj = self.Ubem[j,:].reshape([2, 1])
			increment = fsolve(self.ai, uj, args=(j)).reshape([2, 1])
			self.Ubem[j+1, :] = increment[:, 0]
		ubem = self.Ubem
		return ubem
#
	def SSLS(self,Uzero=[1,1]):
		h = self.Dt
		L = self.L
		R = self.R
		alpha1 = self.a1
		alpha2 = self.a2
		beta1 = self.b1
		beta2 = self.b2
		gamma1 = self.g21
		gamma2 = self.g12
		self.Ussls = np.zeros([L, 2])
		self.Ussls[0] = self.Uzero
		if Uzero[0]!=1:	
			self.Ussls[0] = Uzero
		
		for j in np.arange(L-1):
			self.Winc = np.sum(self.dWj[:,R * (j):R * (j + 1)], axis=1)
			self.Winc = self.Winc.reshape([2, 1])
			uj = self.Ussls[j,:].reshape([2, 1])
			a11 = alpha1 * np.exp(gamma1 * np.log(uj[1,0])) - beta1
			Uj1 = uj[0, 0] * np.exp(h * a11) 
			a12 = alpha2 * np.exp(gamma2 * np.log(Uj1)) - beta2
			Uj2 = uj[1, 0] * np.exp(h * a12)
			Ustar = np.array([Uj1, Uj2]).reshape([2, 1])
			increment = Ustar + np.dot(self.b(Ustar), self.Winc)
			self.Ussls[j+1, :] = increment[:, 0]
		ussls = self.Ussls
		return ussls
	def LongTimeBehavior(self, K=1):
		self.LongTimeM = K
		self.Ust = self.Ussls
		L=self.L
		self.Uas = np.zeros([L,2])
		self.detUas = np.zeros([L,2])
		for j in np.arange(K):
			self.NoiseUpdate(flag=0)
			Uzero1=self.Ussls[-1]
			Uzero2=self.Urk[-1]
			self.Uas = self.SSLS(Uzero1)
			self.detUas = self.RK(Uzero2)
		Uas = self.Uas
		detUas = self.detUas
		return detUas, Uas
#	
	def Moments(self, M,LTM=10):
		self.LTM = LTM
		self.M=M
		L=self.L
		self.M_uls = np.zeros([L,2])
		self.MS_uls = np.zeros([L,2])
		S_uls=np.zeros([M,L,2])
		print '=================================================================='
		for j in np.arange(M):
			if np.mod(j,10)==0:
				print '\t Sample: '+str(j)+' of '+str(M)
			S_uls[j]=self.SSLS()
			self.NoiseUpdate(flag=1)
		self.M_uls = S_uls.mean(axis=0)
		self.MS_uls = (S_uls**2).mean(axis=0)
#
	def SaveData(self):
#
		U1=self.Urk[:, 0]
		U2=self.Urk[:, 1]
		tagPar = np.array([\
			'a1=',\
			'b1=',\
			'a2=',\
			'b2=',\
			'g11=',\
			'g12=',\
			'g21=',\
			'g22=',\
			'sigma1',\
			'sigma2',\
			'k1=',\
			'k2=',\
			'Ubar1=',\
			'Ubar2=',\
			'Uzero1=',\
			'Uzero2=',\
			'k=',\
			'T0=',\
			'N=',\
			'T=',\
			'dt='\
			])
		ParValues=np.array([\
			self.a1,\
			self.b1,\
			self.a2,\
			self.b2,\
			self.g11,\
			self.g12,\
			self.g21,\
			self.g22,\
			self.sigma[0],\
			self.sigma[1],\
			self.k1,\
			self.k2,\
			self.Ubar[0],\
			self.Ubar[1],\
			self.Uzero[0],\
			self.Uzero[1],\
			self.k,\
			self.T0,\
			self.N,\
			self.T,\
			self.dt\
			])
		PARAMETERS=np.column_stack((tagPar,ParValues))
		np.savetxt('ParametersDet.txt', PARAMETERS, delimiter=" ", fmt="%s")
	#	np.savetxt('SolutionDet.txt',\
	#			np.transpose(np.array([t[:,0], U1[:], U2[:]])),\
	#			fmt=['%1.8f','%1.8f','%1.8f'],\
	#			delimiter='\t'\
	#	)
	#
	#Mu1 = self.M_uls[:,0]
	#Mu2 = self.M_uls[:,1]
	#MSu1 = self.MS_uls[:,0]
	#MSu2 = self.MS_uls[:,1]
	#
	#np.savetxt('SolutionMean.txt',\
	#				np.transpose(np.array([t[:,0], Mu1[:] , Mu2[:]])),\
	#			fmt=['%1.8f','%1.8f','%1.8f'],\
	#			delimiter='\t'\
	#	)
	#np.savetxt('SolutionMeanSquare.txt',\
	#			np.transpose(np.array([t[:,0], MSu1[:] , MSu2[:]])),\
	#			fmt=['%1.8f','%1.8f','%1.8f'],\
	#			delimiter='\t'\
	#	)
	#'''