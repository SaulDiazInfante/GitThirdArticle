import numpy as np
from DeterministicModelJerezChen import PLBRMJerezChen
import multiprocessing as mp
from scipy.optimize.minpack import fsolve
import multiprocessing as mp
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
	def NoiseUpdate(self, seed, flag=1):
		np.random.seed(seed)
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
	def TamedEuler(self, Uzero, seed, flag=0):
		h = self.Dt
		L = self.L
		R = self.R
		ss= 100
		if flag == 1:
			Dt = self.dt
			L = self. N
			R =1
		self.Utem = np.zeros([L, 2])
		self.Utem[0] = self.Uzero
		if Uzero.any != self.Uzero.any:
			self.Utem[0] = Uzero
		self.NoiseUpdate(seed)
		for j in np.arange(L - 1):
			self.Winc = np.sum(self.dWj[:,R * (j):R * (j + 1)],axis=1)
			self.Winc = self.Winc.reshape([2,1])
			uj = self.Utem[j,:].reshape([2,1])
			aj = h*self.a(uj)
			increment = uj + aj/(1+np.linalg.norm(aj)) + np.dot(self.b(uj),self.Winc)
			self.Utem[j+1, :] = increment[:, 0]
		utem = self.Utem[0:L:ss,:]
		return utem
#
	def RK(self):
		h = self.Dt
		Dt =self.Dt
		L = self.L
		R = self.R
		self.Urk = np.zeros([L, 2])
		self.Urk[0] = self.Uzero
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
#
	def SSLS(self, seed, Uzero=[1,1]):
		h = self.Dt
		L = self.L
		ss = 100
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
		self.NoiseUpdate(seed)
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
		ussls = self.Ussls[0:L:ss,:]
		return ussls
#
	def LongTimeBehavior(self, seed, K=1):
		self.LongTimeM = K
		L = self.L
		self.Ust = self.SSLS(seed, self.Uzero)
		#self.Ust=self.TamedEuler(seed, self.Uzero, flag=0)
		for j in np.arange(K):
			self.NoiseUpdate(seed, flag=0)
			self.Uas=np.zeros([L,2])
			self.Uas[0]=self.Ussls[-1]
			#self.Uas[0] = self.Utem[-1]
			self.Uas = self.SSLS(seed, self.Uas[0])
			#self.Uas = self.TamedEuler(self.Uas[0], seed, flag=0)
		Uas = self.Uas
		return Uas
#	
	def Moments(self, M, LTM=10):
		self.LTM = LTM
		self.M=M
		L=self.L
		Uzero = self.Uzero
		print '======================Calculating Short Time Moments============================================'
		pool = mp.Pool(processes=8)
		UpathsSymb = [pool.apply_async(self.SSLS, args=(selfUzero, i)) for i in np.random.random_integers(1, 123456789, M)]
		Upaths=[p.get() for p in UpathsSymb]
		Upaths = np.array(Upaths)
		self.meanU = Upaths.mean(axis=0)
		self.meanSU = (Upaths**2).mean(axis=0)
		np.save("MeanShortTime.npy", self.meanU)
		np.save("MeanSquareShortTime.npy", self.meanSU)
		np.save("ShortTimePaths.npy")
		
		print '======================Calculating Long Time Moments============================================'
		UpathsSymb = [pool.apply_async(self.LongTimeBehavior, args=(i,LMT)) for i in np.random.random_integers(1, 123456789, M)]
		Upaths = [p.get() for p in LongUpathsSymb]
		Upaths = np.array(LongUpaths)
		meanU = LongUpaths.mean(axis=0)
		meanSU = (LongUpaths**2).mean(axis=0)
		self.meanU = Upaths.mean(axis=0)
		self.meanSU = (Upaths**2).mean(axis=0)
		np.save("MeanLongTime.npy", self.meanU)
		np.save("MeanLongSquareShortTime.npy", self.meanSU)
		np.save("LongTimePaths.npy")

#
	def SaveData(self):
		t=self.t[0 : -1 : self.R].reshape([self.t[0 : - 1 : self.R].shape[0],1])
		tas = self.LongTimeM * t
		U1=self.Urk[:, 0]
		U2=self.Urk[:, 1]
		stoU1 = self.Ust[:, 0]
		stoU2 = self.Ust[:, 1]
		stoU1as = self.Uas[:, 0]
		stoU2as = self.Uas[:, 1]
		Mu1 = self.M_uls[:,0]
		Mu2 = self.M_uls[:,1]
		MSu1 = self.MS_uls[:,0]
		MSu2 = self.MS_uls[:,1]
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
		np.savetxt('SolutionDet.txt',\
				np.transpose(np.array([t[:,0], U1[:], U2[:]])),\
				fmt=['%1.8f','%1.8f','%1.8f'],\
				delimiter='\t'\
		)
		np.savetxt('SolutionSto.txt',\
				np.transpose(np.array([t[:,0], stoU1[:] , stoU2[:]])),\
				fmt=['%1.8f','%1.8f','%1.8f'],\
				delimiter='\t'\
		)
		np.savetxt('SolutionStoAs.txt',\
				np.transpose(np.array([tas[:,0], stoU1as[:] , stoU2as[:]])),\
				fmt=['%1.8f','%1.8f','%1.8f'],\
				delimiter='\t'\
		)
