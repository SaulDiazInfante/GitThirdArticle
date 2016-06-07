"""
#==============================================================================
#This code calculate the solution of Stochastic Bone remodeling model
#proposal by Chen and Jerez. This code consider the following SDE
	dx(t)=(\alpha_1\exp(y(t))-v1)dt+\sigma_1dW_1(t)
	dy(t)=(\alpha_2\exp(-x(t))-v1)dt+\sigma_2dW_2(t)
where:
	x(t)=ln(u_1(t))	v1=\beta_1+1/2\sigma_{1}^{2}
	y(t)=ln(u_2(t))	v2=\beta_2+1/2\sigma_{2}^{2}
	x_0=ln(u_1(0))	y_0=ln(u_2(0))
#==============================================================================
"""
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
from scipy.integrate import odeint
class StochasticPLBRM:
	def InitializeMesh(self,k,r,T0,T):
	#Stensil of the mesh
		self.k=k
		self.r=r
		self.T0=T0
		self.N=10.0**k
		self.R=10.0**r
		self.T=T
		self.dt = (self.T-self.T0)/np.float(self.N)
		self.IndexN=np.arange(self.N+1)
		self.tau=self.IndexN[0:self.N+1:self.R]
		self.t=np.linspace(0,self.T,self.N)
		#
		self.Dt = np.float(self.R) * self.dt;
		self.L = self.N / self.R;
		#1-dimensional Brownian Paths generation
		'''
		self.DistNormal=np.random.randn(np.int(self.N))
		self.dW = np.sqrt(self.dt)*self.DistNormal
		self.W = np.cumsum(self.dW)
		self.W = np.concatenate(([0], self.W))
		self.DeltaW=np.zeros(self.L)
		#2-dimensional Brownian paths
		'''
		self.DistNormal=np.random.randn(np.int(self.N),2)
		self.dW = np.sqrt(self.dt)*self.DistNormal
		self.W = np.cumsum(self.dW,axis=0)
		self.W = np.concatenate((np.zeros([1,2]), self.W))
		
	def SetParametersStoPLBRM(self,a1,b1,a2,b2,g11,g12,g21,g22,k1,k2,sigma,U0):
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
		self.Ubar=[u1bar,u2bar]
		self.Uzero=[np.log(U0[0]),np.log(U0[1])]
		#self.Uzero=np.array([11.0,5.0,Zzero])
	#
	def SodeSolver(self):
		def a(u):
			v1=self.b1+0.5*self.sigma[0]**2
			v2=self.b2+0.5*self.sigma[1]**2
			u1=self.a1*np.exp(u[1])-v1
			u2=self.a2*np.exp(-u[0])-v2
			return np.array([u1,u2])
		def b(u):
			u1=self.sigma[0]#*(u[0]-self.Ubar[0])
			u2=self.sigma[1]#*(u[1]-self.Ubar[1])
			return np.array([u1,u2])
			# Implementation of EulerScheme.
		self.Uem = np.zeros([self.L,2])
		self.Uem[0]=self.Uzero
		h = self.Dt
		for j in np.arange(self.L-1):
			self.Winc =np.sum(self.dW[self.R*(j):self.R*(j+1)],axis=0)
			aj = a(self.Uem[j])
			bj = b(self.Uem[j])
			eu = h * aj + np.dot(bj, self.Winc)
			norm = LA.norm(eu,2) 
			tamed = np.max([1.0, h*norm])
			increment=self.Uem[j] + eu/tamed
			self.Uem[j+1]=increment
		self.sU1=np.exp(self.Uem[:,0])
		self.sU2=np.exp(self.Uem[:,1])
	#
	def SaveData(self):
		t=self.t[0:-1:self.R].reshape([self.t[0:-1:self.R].shape[0],1])
		U1=self.sU1.reshape([self.sU1.shape[0],1])
		U2=self.sU2.reshape([self.sU2.shape[0],1])
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
			'sigma1=',\
			'sigma2=',
			'Ubar1=',\
			'Ubar2=',\
			'Uzero1=',\
			'Uzero2=',\
			'k=',\
			'r=',\
			'T0=',\
			'N=',\
			'R=',\
			'T=',\
			'dt=',\
			'Dt=',\
			'L='\
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
			self.sigma[0],\
			self.sigma[1],
			self.Ubar[0],\
			self.Ubar[1],\
			np.exp(self.Uzero[0]),\
			np.exp(self.Uzero[1]),\
			self.k,\
			self.r,\
			self.T0,\
			self.N,\
			self.R,\
			self.T,\
			self.dt,\
			self.Dt,\
			self.L\
			])
		PARAMETERS=np.column_stack((tagPar,ParValues))
		np.savetxt('ParametersSto.txt', PARAMETERS, delimiter=" ", fmt="%s")
		np.savetxt('SolutionSto.txt', np.hstack([t,U1,U2]),fmt=['%1.8f','%1.8f','%1.8f'])
	#
	def Plots(self):
		#Load Data.
		A=np.loadtxt('SolutionDet.txt')
		B=np.loadtxt('SolutionSto.txt')
		t=A[:,0]
		U1=A[:,1]
		U2=A[:,2]
		st=B[:,0]
		sU1=B[:,1]
		sU2=B[:,2]
		#
		fig1=plt.figure()
		ax1 = fig1.add_subplot(121)
		ax1.plot(t,U1,'k-',label='osteoclast')
		ax1.plot(t,self.Ubar[0]*np.ones(U1.shape[0]),'r--',label='Steady State')
		ax1.plot(st,sU1,'g-',label='noise')
		ax1.set_xlabel(r't')
		ax1.set_ylabel(r'Osteoclast')
		ax1.grid(True)
		ax1.legend(loc=0)
		ax1.set_title(r'$dx(t)=(\alpha_1\exp(y(t))-v1)dt+\sigma_1dW(t)$')
		#
		ax2 = fig1.add_subplot(122)
		ax2.plot(t,U2,'k-',label='osteoblast')
		ax2.plot(t,self.Ubar[1]*np.ones(U1.shape[0]),'r--',label='Deter Steady State')
		#
		ax2.plot(st,self.sU2,'g-',label='noise')
		ax2.set_xlabel(r't')
		ax2.set_ylabel(r'Osteoblast')
		ax2.grid(True)
		ax2.legend(loc=0)
		plt.savefig('figure1.png')
		#
		fig2=plt.figure()
		ax1 = fig2.add_subplot(111)
		ax1.plot(U1,U2,'k-',label='Deterministic')
		ax1.plot(sU1,sU2,'r-',alpha=0.3,label='noise perturbation')
		ax1.plot(self.Ubar[0],self.Ubar[1],'go',ms=7)
		ax1.plot(sU1[0],sU2[1],'rs',ms=7)
		ax1.set_xlabel(r'u_1')
		ax1.set_ylabel(r'u_2')
		ax1.grid(True)
		ax1.legend(loc=0)
		plt.savefig('figure2.png')