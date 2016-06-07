#======================================================================================================================
#This code calculate the solution of Deterministic Bone Remodeling Model
#proposal by Komarova in the paper
#Mathematical model predicts a critical role for osteoclast autocrine regulation in the control of bone remodeling.
#Here we reproduce the simulation presented by Ayati using odeint from scipy and the Linear Steklov Method. The 
#Komarova model reads
# du_1(t) = \alpha_1 u_1^{\gama_{11}} u_2^{\gamma_{21}}-\beta_1 u_1 + \sigma_1 u_1
# du_2(t) = \alpha_2 u_1^{\gama_{12}} u_2^{\gamma_{22}}-\beta_2 u_2 + \sigma_2 u_2
#======================================================================================================================
import numpy as np
from scipy.integrate import odeint
class PLBRM:
	def InitializeMesh(self,k,T0,T):
	#Stensil of the mesh
		self.k=k
		self.T0=T0
		self.N=10.0**k
		self.T=T
		self.dt = self.T/np.float(self.N)
		self.t=np.linspace(0,self.T,self.N+1)	
	def SetParametersKomarova(self,a1,b1,a2,b2,g11,g12,g21,g22,k1,k2,U0):
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
		#
		#steady states
		self.gamma = g12*g21 - (1.0 - g11) * (1.0 - g22)
		u1bar = ((b1 / a1) ** ((1.0 - g22) / self.gamma)) * (b2 / a2) ** ( g21 / self.gamma)
		u2bar = ((b1 / a1) ** (g12 / self.gamma)) * ((b2 / a2) **( (1.0 - g11)/self.gamma))
		self.Ubar=np.array([u1bar,u2bar])
		self.Uzero=[U0[0],U0[1]]
		#self.Uzero=np.array([11.0,5.0,Zzero])
	def OdeSolve(self):
		def f(u,t):
			u1=self.a1*(u[0]**self.g11)*(u[1]**self.g21)-self.b1*u[0]
			u2=self.a2*(u[0]**self.g12)*(u[1]**self.g22)-self.b2*u[1]
			return np.array([u1,u2])
		self.soln = odeint(f, self.Uzero, self.t)
		self.U1 = self.soln[:, 0]
		self.U2 = self.soln[:, 1]
		# Solve using Runge Kutta. -------------------------------------------------------------------------------------
		self.Xrk = np.zeros([self.N+1,2])
		self.Xrk[0]=self.Uzero
		for j in np.arange(self.N):
			Y1=self.Xrk[j]
			Y2=self.Xrk[j]+0.5*self.dt*f(Y1,self.t[j])
			Y3=self.Xrk[j]+0.5*self.dt*f(Y2,self.t[j])
			Y4=self.Xrk[j]+self.dt*f(Y3,self.t[j])
			self.Xrk[j+1]=self.Xrk[j]\
			+self.dt*(\
			1.0/6.0*f(Y1,self.t[j])+\
			1.0/3.0*f(Y2,self.t[j])+\
			1.0/3.0*f(Y3,self.t[j])+\
			1.0/6.0*f(Y4,self.t[j])\
			)
	def SaveData(self):
		t=self.t
		U1=self.U1.reshape([self.U1.shape[0],1])
		U2=self.U2.reshape([self.U2.shape[0],1])
		tagPar = np.array([\
			'a1=',\
			'b1=',\
			'a2=',\
			'b2=',\
			'g11=',\
			'g12=',\
			'g21=',\
			'g22=',\
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
		np.savetxt('SolutionDet.txt', np.array(np.transpose([t, U1[:, 0] , U2[:, 0]])),fmt=['%1.8f','%1.8f','%1.8f'])
#plt.ion()
