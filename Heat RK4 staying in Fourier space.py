# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 16:23:34 2024

@author: baile
"""

import numpy as np 
import matplotlib.pyplot as plt

def F(uhat, wavenumbers):
    return uhat * -np.power(wavenumbers,2)

def analytic(t, x):   
    u_actual = -np.exp(-t) * np.cos(x)
    plt.plot(x, u_actual, 'r--', linewidth="8")   
    
#Heat Equation conditions
L=2*np.pi
N=105
xpre=np.linspace(0, L, N+1)
x=xpre[0:N]
t_iteration = 850
t=0
dt = 0.001
u=np.zeros([t_iteration,len(x)],dtype=complex)
u0=-np.cos(x)
u[:][0]=np.fft.fft(u0)
wavenumbers = np.concatenate((np.arange(0, N//2), [0], np.arange(-N//2+1, 0)))



for i in range(t_iteration-1): 
    k1= F(u[:][i],wavenumbers)
    k2= F(u[:][i]+dt*k1/2,wavenumbers)
    k3= F(u[:][i]+dt*k2/2,wavenumbers)
    k4= F(u[:][i]+dt*k3,wavenumbers)
    u[:][i+1]=u[:][i] + dt/6 * (k1+2*k2+2*k3+k4)
    t=t+dt
    if i==800:
        uplot=np.fft.ifft(u[:][i]).real
        plt.xlabel("x", fontsize=35)
        plt.ylabel("u(x,t)", fontsize=35)
        plt.xticks(fontsize=20) 
        plt.yticks(fontsize=20) 
        plt.axis([0,L,-1,1])
        plt.tight_layout(rect=[0, 0, 1, 1])
        plt.plot(x, uplot, 'b',linewidth="8")
        analytic(t, x)
        plt.style.use("seaborn")
        plt.show()  