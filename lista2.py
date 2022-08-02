from cmath import log
from re import X
import numpy as np
import scipy.integrate
import math
c = 3*10**5


class cosmo:
    def __init__(self, omega_m,omega_r,omega_L,omega_k,H0):
        self.omega_m = omega_m
        self.omega_r = omega_r
        self.omega_L = omega_L
        self.omega_k = omega_k
        self.H0      = H0 
    
    

    def E(self,z):
        HH0 = (self.omega_L + self.omega_k * (1+z)**2 + self.omega_m * (1+z)**3 + self.omega_r * (1+z)**4)**(1/2)   
        return HH0
    


    def Ea(self, param, z):
        if param == 'H0':
            return 0
        
        if param == 'omega_m':
            return (1+z)**3/self.E(z)
        
        if param == 'omega_r':
            return (1+z)**4/self.E(z)

        if param == 'omega_L':
            return 1/self.E(z)
        
        if param == 'omega_k':
            return (1+z)**2/self.E(z)
        

         


    def Dc(self,z):
        cosmo_model = cosmo(self.omega_m,self.omega_r,self.omega_L,self.omega_k,self.H0)
        def Dc_integrand(z_prime):
            return 1/cosmo_model.E(z_prime)
        D_c = scipy.integrate.quad(Dc_integrand,0,z)[0]
        return D_c

    
    def Dca(self,param,z):
        cosmo_model = cosmo(self.omega_m,self.omega_r,self.omega_L,self.omega_k,self.H0)
        def Dca_integrand(z_prime):
            return cosmo_model.Ea(param,z_prime)/cosmo_model.E(z_prime)**2
        D_ca = -scipy.integrate.quad(Dca_integrand,0,z)[0]
        return D_ca    


    
    def Dt(self, z):
        cosmo_model = cosmo(self.omega_m,self.omega_r,self.omega_L,self.omega_k,self.H0)
        if self.omega_k < 0:
            D_t = np.sinh(np.sqrt(-self.omega_k) * cosmo_model.Dc(z))/(np.sqrt(-self.omega_k))

        if self.omega_k > 0:
            D_t = np.sin(np.sqrt(self.omega_k) * cosmo_model.Dc(z))/(np.sqrt(self.omega_k))

        if self.omega_k == 0:
            D_t = cosmo_model.Dc(z)

        return D_t


    def Dta(self, param ,z):
        cosmo_model = cosmo(self.omega_m,self.omega_r,self.omega_L,self.omega_k,self.H0)
        if self.omega_k < 0:
            D_ta = np.cosh(np.sqrt(-self.omega_k) * cosmo_model.Dc(z))*cosmo_model.Dca(param,z)

        if self.omega_k > 0:
            D_ta = np.cos(np.sqrt(self.omega_k) * cosmo_model.Dc(z))**cosmo_model.Dca(param,z)

        if self.omega_k == 0:
            D_ta = cosmo_model.Dca(param,z)

        return D_ta

    def Dl(self,z):
        cosmo_model = cosmo(self.omega_m,self.omega_r,self.omega_L,self.omega_k,self.H0)
        return cosmo_model.Dt(z) * (1+z)

    def mua(self, param, z):
        cosmo_model = cosmo(self.omega_m,self.omega_r,self.omega_L,self.omega_k,self.H0)
        if param == 'H0':
            return -5/(self.H0*np.log(10))
        
        else:
            mu_a = 5/(cosmo_model.Dt(z)*np.log(10))* cosmo_model.Dta(param,z)
            return mu_a

    


    
def super_gauss(params, mu, sigma, z, arg):
    cosmo_params = {'H0': params[0], 'omega_m': params[1], 'omega_r': params[2], 'omega_k': params[3], 'omega_L': params[4]}
    cosmo_model = cosmo(cosmo_params['omega_m'], cosmo_params['omega_r'], cosmo_params['omega_L'], cosmo_params['omega_k'], cosmo_params['H0'])
        
    if arg == 'func':
        mu_z = []
        x = []
        mu = mu 
        sigma = sigma
        for zi in z:
            mu_z.append(5 * np.log10(cosmo_model.Dl(zi)) + 25 + 5 * np.log10(c/cosmo_params['H0']))    #H0 em m/sMpc
        
        for i in range(len(z)):
            x.append((mu_z[i] - mu[i])**2/sigma[i]**2)
        
        L = np.exp(-1/2*(np.sum(x)))

        return -2*np.log(L) 
        
    if arg == 'params':
        return cosmo_params
        
    else:
        print('Escolha arg como func ou params')

        
def super_gauss_flatr0(params, mu, sigma, z, arg):
    cosmo_params = {'H0': params[0], 'omega_m': params[1], 'omega_L': params[2]}
    cosmo_model = cosmo(cosmo_params['omega_m'], 0, cosmo_params['omega_L'], 0, cosmo_params['H0'])
        

    if arg == 'func':
        mu_z = []
        x = []
        mu = mu 
        sigma = sigma
        for zi in z:
            mu_z.append(5 * np.log10(cosmo_model.Dl(zi)) + 25 + 5 * np.log10(c/cosmo_params['H0']))    #H0 em m/sMpc
        
        for i in range(len(z)):
            x.append((mu_z[i] - mu[i])**2/sigma[i]**2)
        
        L = np.exp(-1/2*(np.sum(x)))

        return -2*np.log(L) 
        
    if arg == 'params':
        return cosmo_params
        
    else:
        print('Escolha arg como func ou params')

    


