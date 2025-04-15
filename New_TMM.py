#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 14 15:14:20 2025

@author: sidharthgoyal
"""

import numpy as np
from math import sin, pi
import matplotlib.pyplot as plt

class TMM_TE:
    def __init__(self, eps, mu, d, wavelength, incidence_angle):
        self.c = 3e8
        self.wave_length = wavelength
        self.omega = (2 * pi * self.c) / self.wave_length
        self.k0 = self.omega / self.c

        self.eps = eps
        self.mu = mu
        self.d = d
        self.inc = incidence_angle

        self.eps0 = 8.85e-12
        self.mu0 = pi * 4e-7
    
    
    def kz(self):
        
        kz = []
        k0 = self.k0
        inc = self.inc
        eps = self.eps
        mu = self.mu
        
        kx = np.sqrt(eps[0]*mu[0]*(k0)**2)*sin(inc)
        
        
        for i in range(0, len(eps)):
            k_mag_squared = eps[i]*mu[i]*(k0)**2 +0j
            
            kz.append(np.sqrt(k_mag_squared - kx**2 +0j))
            '''print('kx', kx**2)
            print('kmag', k_mag_squared)
            print('kz = ', np.sqrt(k_mag_squared - kx**2 +0j))'''
                
                
            
        return kz
    
    def F_matrix(self, kz, d):
        return np.array([
            [np.exp(1j * kz * d), 0],
            [0, np.exp(-1j * kz * d)]
        ], dtype=complex)

    def N_matrix(self, kz, mu):
        return np.array([
            [1, 1],
            [-kz / (self.omega * mu * self.mu0), kz / (self.omega * mu * self.mu0)]
        ], dtype=complex)
    
    
    def Transfer_Matrix(self):
        kz = self.kz()
        mu = self.mu
        d = self.d
        
        T = np.identity(2, dtype = complex)
        for which_corrd in range(0, len(d)):
            #print('This is: ', which_corrd)
            
            F2 = self.F_matrix(kz[which_corrd + 1], d[which_corrd])
            F1 = self.F_matrix(kz[which_corrd], d[which_corrd])
            N2 = self.N_matrix(kz[which_corrd + 1], mu[which_corrd])
            N1 = self.N_matrix(kz[which_corrd], mu[which_corrd])
            
           # print(-kz[which_corrd] / (self.omega * mu[which_corrd] * self.mu0))


            F2 = np.linalg.pinv(F2)
            N2 = np.linalg.pinv(N2)
            
            T = F2@N2@N1@F1@T
            
            
        #print('T', T)  
        
        
        return T
            
    
    def R_T_coeff(self):
        T = self.Transfer_Matrix()
        
        
        t11 = T[0,0]
        t12 = T[0,1]
        t21 = T[1,0]
        t22 = T[1,1]
        
        R = abs(t21/t22)**2
        T = 1 - R
            
            
        
        #print(R)
        
        return R,T
    
    def plot_R_T(self, angle_range):
        R_values = []
        T_values = []
        angles = []
        

        init_angle = 0
        fin_angle = angle_range*pi/180
        
        while init_angle <= fin_angle:
            self.inc = init_angle
            coeff = self.R_T_coeff()
            R_values.append(coeff[0])
            T_values.append(coeff[1])
            angles.append(init_angle*180/pi)
            
            init_angle+=0.001
            
        return R_values,T_values,angles
    
class TMM_TM:
    def __init__(self, eps, mu, d, wavelength, incidence_angle):
        self.c = 3e8
        self.wave_length = wavelength
        self.omega = (2 * pi * self.c) / self.wave_length
        self.k0 = self.omega / self.c

        self.eps = eps
        self.mu = mu
        self.d = d
        self.inc = incidence_angle

        self.eps0 = 8.85e-12
        self.mu0 = pi * 4e-7
    
    
    def kz(self):
        
        kz = []
        k0 = self.k0
        inc = self.inc
        eps = self.eps
        mu = self.mu
        
        kx = np.sqrt(eps[0]*mu[0]*(k0)**2)*sin(inc)
        
        
        for i in range(0, len(eps)):
            k_mag_squared = eps[i]*mu[i]*(k0)**2 +0j
            
            kz.append(np.sqrt(k_mag_squared - kx**2 +0j))
            '''print('kx', kx**2)
            print('kmag', k_mag_squared)
            print('kz = ', np.sqrt(k_mag_squared - kx**2 +0j))'''
                
                
            
        return kz
    
    def F_matrix(self, kz, d):
        return np.array([
            [np.exp(1j * kz * d), 0],
            [0, np.exp(-1j * kz * d)]
        ], dtype=complex)

    def N_matrix(self, kz, eps):
        return np.array([
            [1, 1],
            [-kz / (self.omega * eps * self.eps0), kz / (self.omega * eps * self.eps0)]
        ], dtype=complex)
    
    
    def Transfer_Matrix(self):
        kz = self.kz()
        eps = self.eps
        d = self.d
        
        T = np.identity(2, dtype = complex)
        for which_corrd in range(0, len(d)):
            #print('This is: ', which_corrd)
            
            F2 = self.F_matrix(kz[which_corrd + 1], d[which_corrd])
            F1 = self.F_matrix(kz[which_corrd], d[which_corrd])
            N2 = self.N_matrix(kz[which_corrd + 1], eps[which_corrd + 1])
            N1 = self.N_matrix(kz[which_corrd], eps[which_corrd])
            #print(eps[which_corrd])
            
           # print(-kz[which_corrd] / (self.omega * mu[which_corrd] * self.mu0))


            F2 = np.linalg.pinv(F2)
            N2 = np.linalg.pinv(N2)
            
            T = F2@N2@N1@F1@T
            
            
        #print('T', T)  
        
        
        return T
            
    
    def R_T_coeff(self):
        T = self.Transfer_Matrix()
        
        
        t11 = T[0,0]
        t12 = T[0,1]
        t21 = T[1,0]
        t22 = T[1,1]
        
        R = abs(t21/t22)**2
        T = 1 - R
            
            
        
        #print(R)
        
        return R,T
    
    def plot_R_T(self, angle_range):
        R_values = []
        T_values = []
        angles = []
        

        init_angle = 0
        fin_angle = angle_range*pi/180
        
        while init_angle <= fin_angle:
            self.inc = init_angle
            coeff = self.R_T_coeff()
            R_values.append(coeff[0])
            T_values.append(coeff[1])
            angles.append(init_angle*180/pi)
            
            init_angle+=0.001
            
        return R_values,T_values,angles
            

## Testing
eps = [2.25,-10+1j,1]
#eps = [1.3,1,1.3]
mu = [1,1,1]
d = [0,25e-9]
wavelength = 1e-6
incidence_angle = 89


TMM = TMM_TE(eps, mu, d, wavelength, incidence_angle)   
#print(TMM.kz())
TMM2 = TMM_TM(eps, mu, d, wavelength, incidence_angle)
#print(TMM.kz(eps, mu))
#print(np.linalg.inv(TMM.F_matrix(4360919.820072095j, 25e-9)))

R_values, T_values, angles = TMM.plot_R_T(incidence_angle)

R_val_TM = TMM2.plot_R_T(incidence_angle)[0]
T_val_TM = TMM2.plot_R_T(incidence_angle)[1]


#print(angles)
#print(R_values)

plt.figure(figsize=(8, 6))
plt.plot(angles, R_values, label="Reflection_TE (R)", color='red')
plt.plot(angles, R_val_TM, label="Reflection_TM (R)", color='blue')
#plt.plot(angles, T_values, label="Transmission_TE (T)", color='blue')
plt.xlabel("Incidence Angle (degrees)")
plt.ylabel("Coefficient Value")
plt.title("Reflection Coefficients vs. Incidence Angle")
plt.legend()
plt.grid(True)
plt.show() 

plt.figure(figsize=(8, 6))
#plt.plot(angles, R_val_TM, label="Reflection_TM (R)", color='red')
plt.plot(angles, T_values, label="Transmission_TE (T)", color='red')
plt.plot(angles, T_val_TM, label="Transmission_TM (T)", color='blue')
plt.xlabel("Incidence Angle (degrees)")
plt.ylabel("Coefficient Value")
plt.title("Transmission Coefficients vs. Incidence Angle")
plt.legend()
plt.grid(True)
plt.show()
            