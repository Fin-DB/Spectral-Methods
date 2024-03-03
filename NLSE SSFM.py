# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 14:20:14 2024

@author: baile
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D 
    
def nlin(u,dt):
    modu=np.abs(u)
    exp=np.exp(1j*dt*modu**2)
    return exp*u

def F(nonlinear,k,dt):
    temp=np.fft.fft(nonlinear)  
    sol=np.exp(-1j*dt*k**2)
    return np.fft.ifft(sol*temp)

def analytic(x,t):
    plt.xlabel("x", fontsize=35)
    plt.ylabel("u(x,t)", fontsize=35)
    plt.xticks(fontsize=20) 
    plt.yticks(fontsize=20) 
    plt.tight_layout(rect=[0, 0, 1, 1])
    plt.plot(x, np.abs(np.sqrt(2)*np.e**(0.5j*(x+0.75*t))*1/np.cosh(x-t)), 'r', linewidth="6")
    #plt.plot(x, np.real(1-4*(1+2j*t)/(1+4*x**2+4*t**2)*np.exp(1j*t)), 'r')
    plt.style.use("seaborn")
    
#Heat Equation conditions
L=20
N=100
dx=L/N
xpre=np.linspace(-L/2,L/2, N+1)
x=xpre[0:N]
t_iteration = 1001
t=0
dt = 0.001
#tfinal=t_iteration*dt-1
u=np.zeros([t_iteration,len(x)],dtype=complex)

u[:][0]=np.sqrt(2)*np.e**(0.5j*x)*1/np.cosh(x)
wavenumbers = (2*np.pi/L)*np.concatenate((np.arange(0, N//2), [0], np.arange(-N//2+1, 0)))
#x=np.arange(0,L,dx)
#t=np.arange(0,tfinal,dt)
#xs, ys = np.meshgrid(x, t)

for i in range(t_iteration-1): 
    if i%100==0:
        plt.xlabel("x", fontsize=35)
        plt.ylabel("u(x,t)", fontsize=35)
        plt.xticks(fontsize=20) 
        plt.yticks(fontsize=20) 
        #plt.axis([-L/2,L/2,0,1])
        plt.tight_layout(rect=[0, 0, 1, 1])
        analytic(x,t)
        probdis=np.abs(u[:][i])
        plt.plot(x, np.abs(u[:][i]), 'b--', linewidth="6")
        plt.style.use("seaborn")
        plt.savefig("{} burger video.png".format(i))
        plt.show()        
    nonlinear=nlin(u[:][i],dt)
    #print(nonlinear)
    u[:][i+1]=F(nonlinear, wavenumbers,dt)
    t=t+dt
    
#fig = plt.figure()
#ax = Axes3D(fig)
#plt.tight_layout(rect=[0, 0, 1, 1])
#ax.plot_surface(xs, ys, u, rstride=1, cstride=1, cmap='jet')
#plt.savefig("2D heat equation {}.png".format(t))
#plt.show()   