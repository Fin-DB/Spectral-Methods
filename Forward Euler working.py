# -*- coding: utf-8 -*-
"""
Created on Tue Jan  2 02:27:15 2024

@author: baile
"""

import numpy as np
import matplotlib.pyplot as plt

def analytic(t, x):   
    u_actual = np.exp(-4*t) * np.cos(2*x)
    plt.plot(x, u_actual, 'r')    
    
def heat(u, wavenumbers,a):
    #Using FFT for derivatives and forward Euler for time stepping
    u_hat=np.fft.fft(u)
    v_hat= u_hat * -(wavenumbers)**2
    u_deriv=np.fft.ifft(v_hat).real
    u_new = u + dt*a*u_deriv
    return u_new

def accuracy(t, x, u):
    u_actual = np.exp(-4*t) * np.cos(2*x)
    if max(abs(u_actual-u))>1:
        exit()
        
#Heat Equation conditions
L=2*np.pi
a=1
N=64
xpre=np.linspace(0, L, N+1)
x=xpre[0:N]
u=np.cos(2*x)
t_iteration = 122
t=0
dt = 0.0025
wavenumbers = np.concatenate((np.arange(0, N//2), [0], np.arange(-N//2+1, 0)))
    
for i in range(t_iteration): 
    if i==20:
        #Plotting our solution for u and the analytic solution
        # plt.title("Solution at t={0:.4f}".format(t),fontsize=30)
        # plt.axis([0,L,-1,1])
        # plt.xticks(fontsize=16) 
        # plt.yticks(fontsize=16) 
        # plt.tight_layout(rect=[0, 0, 1, 1])
        # plt.style.use("seaborn")
        # plt.plot(x, u, 'b')
        # analytic(t, x)
        # plt.savefig("Euler error")
        # plt.show()
        print(t)
        print(max(abs(u-np.exp(-4*t) * np.cos(2*x))))
    u=heat(u,wavenumbers,a)
    t=t+dt
    
    