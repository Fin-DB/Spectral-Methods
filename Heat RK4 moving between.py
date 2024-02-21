# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 16:31:06 2024

@author: baile
"""

import numpy as np 
import matplotlib.pyplot as plt

  
def analytic(t, x):   
    u_actual = -np.exp(-t) * np.cos(x)
    plt.plot(x, u_actual, 'r--', linewidth="8") 
    
def diff(u, wavenumbers):
    uhat=np.fft.fft(u)
    return uhat * -np.power(wavenumbers,2)

#Heat Equation conditions
L=2*np.pi
a=1
N=105
xpre=np.linspace(0, L, N+1)
x=xpre[0:N]
t_iteration = 850
t=0
dt = 0.001
u=np.zeros([t_iteration,len(x)])
u[:][0]=-np.cos(x)
wavenumbers = np.concatenate((np.arange(0, N//2), [0], np.arange(-N//2+1, 0)))


for i in range(t_iteration-1): 
    k1= np.fft.ifft(diff(u[:][i],wavenumbers)).real
    k2= np.fft.ifft(diff(u[:][i]+dt*k1/2,wavenumbers)).real
    k3= np.fft.ifft(diff(u[:][i]+dt*k2/2,wavenumbers)).real
    k4= np.fft.ifft(diff(u[:][i]+dt*k3,wavenumbers)).real
    u[:][i+1]=u[:][i] + dt/6 * (k1+2*k2+2*k3+k4)
    t=t+dt
    if i==800:
        plt.xlabel("x", fontsize=35)
        plt.ylabel("u(x,t)", fontsize=35)
        plt.xticks(fontsize=20) 
        plt.yticks(fontsize=20) 
        plt.axis([0,L,-1,1])
        plt.tight_layout(rect=[0, 0, 1, 1])
        plt.plot(x, u[:][i+1], 'b',linewidth="8")
        analytic(t, x)
        plt.style.use("seaborn")
        plt.show()  