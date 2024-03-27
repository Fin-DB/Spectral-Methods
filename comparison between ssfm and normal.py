# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 18:05:17 2024

@author: baile
"""

import numpy as np
import matplotlib.pyplot as plt
import time


def ssfm():
    def nlin(u,dt):
        modu=np.abs(u)
        modusquare=np.fft.ifft(dealiaseusquare(modu)).real
        exp=np.exp(1j*dt*modusquare)
        return exp*u

    def F(nonlinear,k,dt):
        temp=np.fft.fft(nonlinear)  
        sol=np.exp(-1j*dt*k**2)
        return np.fft.ifft(sol*temp)


    def dealiaseusquare(u): 
        uhat=np.fft.fft(u)
        N=len(uhat)
        K = int(3/2 *N)
        uhatpad = np.zeros(K,dtype=complex)
        indvpad = np.concatenate([np.arange(1, N//2 + 1), np.arange(K - N//2 + 1, K + 1)])
        uhatpad[indvpad-1]=uhat
        temp = np.fft.fft((np.fft.ifft(uhatpad)**2).real)
        temp = K/N * temp
        Final = temp[indvpad-1]
        Final[N//2] = 0
        return Final


    def analytic(x,t):
        plt.xlabel("x", fontsize=35)
        plt.ylabel("u(x,t)", fontsize=35)
        plt.xticks(fontsize=20) 
        plt.yticks(fontsize=20) 
        plt.tight_layout(rect=[0, 0, 1, 1])
        plt.plot(x, np.abs(np.sqrt(2)*np.e**(0.5j*(x+0.75*t))*1/np.cosh(x-t)), 'r', linewidth="6")
        plt.style.use("seaborn")
        
    #Initial conditions
    L=20
    N=128
    dx=L/N
    xpre=np.linspace(-L/2,L/2, N+1)
    x=xpre[0:N]
    t_iteration = 3002
    t=0
    dt = 0.001
    u=np.zeros([t_iteration,len(x)],dtype=complex)
    u[:][0]=np.sqrt(2)*np.e**(0.5j*x)*1/np.cosh(x)
    wavenumbers = (2*np.pi/L)*np.concatenate((np.arange(0, N//2), [0], np.arange(-N//2+1, 0)))


    for i in range(t_iteration-1):        
        nonlinear=nlin(u[:][i],dt)
        u[:][i+1]=F(nonlinear, wavenumbers,dt)
        t=t+dt  

        
def normal():
    def F(u, wavenumbers):
        pre=-wavenumbers**2*np.fft.fft(u)
        return 1j*np.fft.ifft(pre) + 1j*np.fft.ifft(dealiaseusquare(np.abs(u))).real*u

    def dealiaseusquare(u): 
        uhat=np.fft.fft(u)
        N=len(uhat)
        K = int(3/2 *N)
        uhatpad = np.zeros(K,dtype=complex)
        indvpad = np.concatenate([np.arange(1, N//2 + 1), np.arange(K - N//2 + 1, K + 1)])
        uhatpad[indvpad-1]=uhat
        temp = np.fft.fft((np.fft.ifft(uhatpad)**2).real)
        temp = K/N * temp
        Final = temp[indvpad-1]
        Final[N//2] = 0
        return Final

    def analytic(x,t):
        plt.xlabel("x", fontsize=35)
        plt.ylabel("|A|", fontsize=35)
        plt.xticks(fontsize=20) 
        plt.yticks(fontsize=20) 
        plt.tight_layout(rect=[0, 0, 1, 1])
        plt.plot(x, np.abs(np.sqrt(2)*np.e**(0.5j*(x+0.75*t))*1/np.cosh(x-t)), 'r--', linewidth="6")
        #plt.plot(x, np.real(1-4*(1+2j*t)/(1+4*x**2+4*t**2)*np.exp(1j*t)), 'r')
        plt.style.use("seaborn")
        

    #Initial conditions
    L=20
    N=128
    dx=L/N
    xpre=np.linspace(-L/2,L/2, N+1)
    x=xpre[0:N]
    t_iteration = 3002
    t=0
    dt = 0.001
    u=np.zeros([t_iteration,len(x)],dtype=complex)
    u[:][0]=np.sqrt(2)*np.e**(0.5j*x)*1/np.cosh(x)
    wavenumbers = (2*np.pi/L)*np.concatenate((np.arange(0, N//2), [0], np.arange(-N//2+1, 0)))



    for i in range(t_iteration-1): 
        k1= F(u[:][i],wavenumbers)
        k2= F(u[:][i]+dt*k1/2,wavenumbers)
        k3= F(u[:][i]+dt*k2/2,wavenumbers)
        k4= F(u[:][i]+dt*k3,wavenumbers)
        u[:][i+1]=u[:][i] + dt/6 * (k1+2*k2+2*k3+k4)
        t=t+dt
        
        
def timer():
    #Prints the time taken between two methods 
    loop=10
    start_move = time.time()
    for i in range (loop):
        normal()   
    normaltime=time.time() - start_move
    
    start_stay = time.time()
    for i in range (loop):
        ssfm()  
    ssfmtime=time.time() - start_stay
    print("normal ", normaltime)
    print("ssfm ", ssfmtime)
    
timer()  