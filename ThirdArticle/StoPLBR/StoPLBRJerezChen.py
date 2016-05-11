import numpy as np
from scipy.optimize.minpack import fsolve
from scipy.integrate import ode
from DeterministicModelJerezChen import PLBRMJerezChen
import sys as Sys

class StoPLBRM(PLBRMJerezChen):
    def __init__(self, k=5, p=0, r=0, t_0=0, t_f=36250, flag=1,
                 a1=0.3, b1=0.2, a2=0.18, b2=0.02, g11=1.0, g12=0.5,
                 g21=-0.9, g22=1.0, k1=0.022914, k2=0.0000038,
                 sigma_1=0.02, sigma_2=0.002, u_0=10.0, v_0=0.7):

        self.k = k
        self.p = p
        self.r = r
        self.T0 = t_0
        self.N = np.int64(10 ** k)
        self.P = np.int64(10 ** p)
        self.R = np.int64(10 ** r)
        # self.N = 2.0 ** k
        # self.P = 2.0 ** p
        # self.R = 2.0 ** r
        self.T = t_f
        self.dt = self.T / np.float(self.N)
        self.IndexN = np.arange(self.N, dtype=np.uint64)
        # set of index to Ito integral --------------------------------
        self.tau = self.IndexN[0:self.N:self.P]
        self.t = np.linspace(0, self.T, self.N , dtype=np.float32)
        self.t_k = self.t[self.tau]
        self.Dt = np.float(self.R) * self.dt
        self.L = self.N / self.R

        self.dWj = np.random.randn(2, np.int64(self.N ))
        self.dWj = np.sqrt(self.dt) * self.dWj
        if flag == 1:
            self.dWj[:, 0] = 0.0
        self.a1 = a1
        self.b1 = b1
        self.a2 = a2
        self.b2 = b2
        self.g11 = g11
        self.g12 = g12
        self.g21 = g21
        self.g22 = g22
        self.k1 = k1
        self.k2 = k2
        self.sigma = np.array([sigma_1, sigma_2])
        # steady states
        self.gamma = g12 * g21 - (1.0 - g11) * (1.0 - g22)
        #
        self.u_bar = np.array([0.0, 0.0])
        self.u_zero = np.array([u_0, v_0])
        self.ji_1 = ((b1 + 0.5 * (sigma_1 ** 2)) / a1) ** (1.0 / g21)
        self.ji_2 = ((b2 + 0.5 * (sigma_2 ** 2)) / a2) ** (1.0 / g12)
        # arrays for methods
        self.u_em = np.zeros([self.L, 2])
        self.u_tem = np.zeros([self.L, 2])
        self.u_rk = np.zeros([self.L, 2])
        self.u_bem = np.zeros([self.L, 2])
        self.u_ssls = np.zeros([self.L, 2])
        # Long time simulation
        self.u_as = np.zeros([self.L, 2])
        self.det_uas = np.zeros([self.L, 2])
        # Numerical methods parameters
        self.long_time_m = 2 * 8 * 650
        self.w_inc = np.array([0.0, 0.0])
        # Bone Mass
        self.z = np.zeros(self.L)
        self.z_sto = np.zeros(self.L)
        self.k_t = 0

    # Print iterations progress
    def print_progress(self, iteration, total, prefix='', suffix='',
                      decimals=2, bar_length=100, ratio=False):
        """
        Call in a loop to create terminal progress bar
        @params:
            iteration   - Required  : current iteration (Int)
            total       - Required  : total iterations (Int)
            prefix      - Optional  : prefix string (Str)
            suffix      - Optional  : suffix string (Str)
        """
        filled_length = int(round(bar_length * iteration / float(
                total)))
        percents = round(100.00 * (iteration / float(total)), decimals)
        bar = '#' * filled_length + '-' * (bar_length - filled_length)
        if ratio:
            Sys.stdout.write(
                    '%s [%s] %s%s%s %s\r' %
                    (prefix, bar, iteration, '/', total, suffix)),
            Sys.stdout.flush()
        else:
            Sys.stdout.write('%s [%s] %s%s %s\r' % (prefix, bar,
                                                    percents, '%',
                                                    suffix)),
            Sys.stdout.flush()
        if iteration == total:
            print("\n")

    # noinspection PyTypeChecker,PyUnresolvedReferences

    def initialize_mesh(self, k, p, r, t_0, t_f):
        self.k = k
        self.p = p
        self.r = r
        self.T0 = t_0
        self.N = np.int64(10 ** k)
        self.P = np.int64(10 ** p)
        self.R = np.int64(10 ** r)
        # self.N = 2.0 ** k
        # self.P = 2.0 ** p
        # self.R = 2.0 ** r
        self.T = t_f
        self.dt = self.T / np.float(self.N)
        self.IndexN = np.arange(self.N, dtype=np.float64)
        # set of index to Ito integral --------------------------------
        self.tau = self.IndexN[0: self.N: self.P]
        self.t = np.linspace(0, self.T, self.N )
        self.Dt = np.float(self.R) * self.dt
        self.L = self.N / self.R

    #
    def noise_update(self, flag=True):
        # flag: boolean
        self.dWj = np.random.randn(2, np.int(self.N))
        self.dWj = np.sqrt(self.dt) * self.dWj
        if flag == True:
            self.dWj[:, 0] = 0.0
        #

    def set_parameters_sto_plbrm(self, a1, b1, a2, b2, g11, g12, g21,
                                 g22, k1, k2, sigma, u0):
        self.a1 = a1
        self.b1 = b1
        self.a2 = a2
        self.b2 = b2
        self.g11 = g11
        self.g12 = g12
        self.g21 = g21
        self.g22 = g22
        self.k1 = k1
        self.k2 = k2
        self.sigma = sigma
        # steady states
        self.gamma = g12 * g21 - (1.0 - g11) * (1.0 - g22)
        u1bar = (b1 / a1) ** ((1.0 - g22) / self.gamma) * \
                 (b2 / a2) ** (g21 / self.gamma)
        u2bar = ((b1 / a1) ** (g12 / self.gamma)) * \
                 ((b2 / a2) ** ((1.0 - g11) / self.gamma))
        self.u_bar = np.array([u1bar, u2bar])
        self.u_zero = np.array([u0[0], u0[1]])

    #
    def a(self, u):
        # type: (object) -> np.array
        # ==============================================================
        #         Diffusion
        # ==============================================================
        alpha1 = self.a1
        alpha2 = self.a2
        beta1 = self.b1
        beta2 = self.b2
        gamma1 = self.g21
        gamma2 = self.g12
        u1 = u[0]
        u2 = u[1]
        u_1 = u1 * (alpha1 * np.exp(gamma1 * np.log(u2)) - beta1)
        u_2 = u2 * (alpha2 * np.exp(gamma2 * np.log(u1)) - beta2)
        return np.array([u_1, u_2]).reshape([2, 1])

    #
    def b(self, u):
        sigma1 = self.sigma[0]
        sigma2 = self.sigma[1]
        u1 = u[0, 0]
        u2 = u[1, 0]
        return np.diag([sigma1 * u1, sigma2 * u2])

    #
    def ai(self, x, j):
        h = self.Dt
        d_w = self.w_inc
        x0 = self.u_bem[j].reshape([2, 1])
        x_0 = x0 + np.dot(self.b(x0), d_w)
        fx = self.a(x)
        return h * fx.ravel() + x_0.ravel() - x.ravel()

    def euler(self, flag=0):
        h = self.Dt
        l = self.L
        r = self.R
        if flag == 1:
            dt = self.dt
            l = self.N
            r = 1
        self.u_em[0] = self.Uzero
        for j in np.arange(l - 1):
            self.w_inc = np.sum(self.dWj[:, r * j:r * (j + 1)], axis=1)
            self.w_inc = self.w_inc.reshape([2, 1])
            uj = self.u_em[j, :].reshape([2, 1])
            aj = self.a(uj)
            increment = uj + h * aj + np.dot(self.b(uj), self.w_inc)
            self.u_em[j + 1, :] = increment[:, 0]
        uem = self.u_em
        return uem

    # noinspection PyTypeChecker
    def tamed_euler(self, uzero, seed, flag=0):
        h = self.Dt
        l = self.L
        r = self.R
        ss = 100
        if flag == 1:
            d_t = self.dt
            l = self.N
            r = 1

        self.u_tem[0] = self.u_zero
        if uzero.any != self.u_zero.any:
            self.u_tem[0] = uzero
        self.noise_update(seed)
        for j in np.arange(l - 1):
            self.w_inc = np.sum(self.dWj[:, r * j:r * (j + 1)], axis=1)
            self.w_inc = self.w_inc.reshape([2, 1])
            uj = self.u_tem[j, :].reshape([2, 1])
            aj = h * self.a(uj)
            increment = uj + aj / (1 + np.linalg.norm(aj)) + np.dot(
                self.b(uj), self.w_inc)
            self.u_tem[j + 1, :] = increment[:, 0]
        u_tem = self.u_tem[0:l:ss, :]
        return u_tem

    def rk(self, uzero=np.array([1.0, 1.0])):
        ss = 1
        h = self.Dt
        l = self.L
        self.u_rk[0] = self.u_zero
        if np.int64(uzero[0]) != np.int64(1.0):
            self.u_rk[0] = uzero
        for j in np.arange(l - 1):
            uj = self.u_rk[j].reshape([2, 1])
            k1 = h * self.a(uj)
            k2 = h * self.a(uj + 0.5 * k1)
            k3 = h * self.a(uj + k2 / 2.0)
            k4 = h * self.a(uj + k3)
            increment = 1 / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
            self.u_rk[j + 1, :] = uj[:, 0] + increment[:, 0]
        urk = self.u_rk[0:l:ss, :]
        return urk

    def bem(self):
        l = self.L
        r = self.R

        self.u_bem[0] = self.u_zero
        for j in np.arange(l - 1):
            self.w_inc = np.sum(self.dWj[:, r * (j): r * (j + 1)],
                                axis=1)
            self.w_inc = self.w_inc.reshape([2, 1])
            uj = self.u_bem[j, :].reshape([2, 1])
            increment = fsolve(self.ai, uj, args=j).reshape([2, 1])
            self.u_bem[j + 1, :] = increment[:, 0]
        u_bem = self.u_bem
        return u_bem

    #
    def ssls(self, seed=123456789, uzero=np.array([1, 1]), fn=1.0):
        h = self.Dt
        l = self.L
        ss = 1
        r = self.R
        alpha1 = self.a1
        alpha2 = self.a2
        beta1 = self.b1
        beta2 = self.b2
        gamma1 = self.g21
        gamma2 = self.g12

        self.u_ssls[0] = self.u_zero
        if uzero[0] != 1:
            self.u_ssls[0] = uzero
        for j in np.arange(l - 1):
            self.w_inc = np.sum(self.dWj[:, r * j:r * (j + 1)], axis=1)
            self.w_inc = self.w_inc.reshape([2, 1])
            uj = self.u_ssls[j, :].reshape([2, 1])
            a11 = alpha1 * np.exp(gamma1 * np.log(uj[1, 0])) - beta1
            uj1 = uj[0, 0] * np.exp(h * a11)
            a12 = alpha2 * np.exp(gamma2 * np.log(uj1)) - beta2
            uj2 = uj[1, 0] * np.exp(h * a12)
            ustar = np.array([uj1, uj2]).reshape([2, 1])
            increment = ustar + fn * np.dot(self.b(ustar), self.w_inc)
            self.u_ssls[j + 1, :] = increment[:, 0]
        u_ssls = self.u_ssls
        return u_ssls

    #
    def long_time_behavior(self, seed, k_times=1):
        """

        :type k_times: np.int64
        :type seed: np.float32
        """
        self.long_time_m = np.int64(k_times)
        for j in np.arange(k_times):
            self.u_as[0] = self.Ussls[-1]
            self.det_uas[0] = self.u_rk[-1]
            self.noise_update(flag=False)
            self.u_as = self.ssls(seed, self.u_as[0])
            self.det_uas = self.rk(self.det_uas[0])
        u_as = self.u_as
        det_uas = self.det_uas
        return u_as, det_uas

    #
    def f(self, t, x):
        k1 = self.k1
        k2 = self.k2
        ji_1 = self.ji_1
        ji_2 = self.ji_2
        u0 = self.u_bar[0]
        v0 = self.u_bar[1]
        k = self.k_t
        t_k = self.t_k[k]
        while  t_k < t:
            k += 1
            t_k = self.t_k[k]
        u = self.u_rk[k, 0]
        v = self.u_rk[k, 1]
        self.k_t = k

        z = - k1 * np.sqrt(np.abs(u - ji_2)) \
            + k2 * np.sqrt(np.abs(v - ji_1))
        """
        z = - k1 * np.max([0.0, u - u0]) \
            + k2 * np.max([0.0, v - v0])

        z = - k1 * np.sqrt(np.max([0.0, u - ji_2])) \
            + k2 * np.sqrt(np.max([0.0, v - ji_1]))
        """
        return z

    def bone_mass(self, k1_sto=0.00004):
        h = self.Dt
        l = self.L
        r = self.R
        k1 = self.k1
        k2 = self.k2
        self.z[0] = .96
        self.z_sto[0] = .96

        integrator = ode(self.f)
        integrator.set_integrator('dopri5')
        integrator.set_initial_value(self.z[0], 0)
        j = 0
        for t in self.t_k[1:]:
            ans = integrator.integrate(t)
            if not integrator.successful():
                print "fuck the integrator does not work"
                self.z[j+1] = self.z[j]
            else:
                self.z[j+1] = ans
            j += 1
        """"
        for j in np.arange(l - 1):
            uj = self.u_rk[j]
            k1 = self.f(uj)
            k2 = self.f(uj + 0.5 * k1)
            k3 = self.f(uj + k2 / 2.0)
            k4 = self.f(uj + k3)
            increment = 1.0 / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
            self.z[j + 1] = self.z[j] + h * increment
            # self.z_sto[j + 1] = self.z_sto[j] + h * z_sto
        """
    #
    def save_data(self):
        t = self.t[0: -1: self.R].reshape(
                [self.t[0: - 1: self.R].shape[0], 1])
        # tas = self.LongTimeM * t
        # self.RK()
        u_0 = self.u_zero
        # U1 = self.u_rk[:, 0]
        # U2 = self.u_rk[:, 1]
        u_stk = self.ssls(123456789, u_0, fn=0.0)
        u_1 = u_stk[:, 0]
        u_2 = u_stk[:, 1]
        #
        sto_u1 = self.u_ssls[:, 0]
        sto_u2 = self.u_ssls[:, 1]
        # stoU1as = self.Uas[:, 0]
        # stoU2as = self.Uas[:, 1]
        # Mu1 = self.M_uls[:,0]
        # Mu2 = self.M_uls[:,1]
        # MSu1 = self.MS_uls[:,0]
        # MSu2 = self.MS_uls[:,1]
        tag_par = np.array([
            'a1=',
            'b1=',
            'a2=',
            'b2=',
            'g11=',
            'g12=',
            'g21=',
            'g22=',
            'sigma1',
            'sigma2',
            'k1=',
            'k2=',
            'Ubar1=',
            'Ubar2=',
            'Uzero1=',
            'Uzero2=',
            'k=',
            'T0=',
            'N=',
            'T=',
            'dt='
            ])
        par_values = np.array([
            self.a1,
            self.b1,
            self.a2,
            self.b2,
            self.g11,
            self.g12,
            self.g21,
            self.g22,
            self.sigma[0],
            self.sigma[1],
            self.k1,
            self.k2,
            self.u_bar[0],
            self.u_bar[1],
            self.u_zero[0],
            self.u_zero[1],
            self.k,
            self.T0,
            self.N,
            self.T,
            self.dt
            ])
        parameters = np.column_stack((tag_par, par_values))
        np.savetxt('ParametersDet.txt', parameters, delimiter=" ",
                   fmt="%s")
        np.save('SolutionDet.npy',
                np.transpose(np.array([t[:, 0], u_1[:], u_2[:]]))
                )
        np.save('SolutionSto.npy',
                np.transpose(np.array([t[:, 0], sto_u1[:], sto_u2[:]]))
                )
        # np.savetxt('SolutionStoAs.txt',\
        # np.transpose(np.array([tas[:,0], stoU1as[:] , stoU2as[:]])),\
        # fmt=['%1.8f','%1.8f','%1.8f'],\
        # delimiter='\t'\
        # )
