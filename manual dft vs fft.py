# -*- coding: utf-8 -*-
"""
Created on Tue Dec 26 00:18:05 2023

@author: baile
"""

import numpy as np
import matplotlib.pyplot as plt
import time


def f(x):
    return np.sin(2*x)

def dft(a): 
    N = len(a)
    a_tilde = np.zeros_like(a, dtype=complex)
    #Compute using standard DFT formula
    for k in range(len(a)):
        m = np.arange(N)
        twiddle_factors = np.exp(-2j*np.pi*m*k/N)
        a_tilde[k] = np.sum(a*twiddle_factors)
    return a_tilde


def fft(x):
    N = len(x)  
    #Return when down to a single data point
    if N == 1:
        return x
    
    #Split the sequence into even and odd parts
    even = fft(x[0::2])
    odd = fft(x[1::2])

    #Compute the twiddle factors
    T = [np.exp(-2j*np.pi*k/N) * odd[k] for k in range(N//2)]
    
    #Butterfly operation
    return [even[k] + T[k] for k in range(N//2)] + [even[k] - T[k] for k in range(N//2)]

def timer(N):
    #Plots time taken for FFT vs DFT
    loop=50
    start_dft = time.time()
    for i in range (loop):
        dft(a)   
    dfttime=time.time() - start_dft
    
    start_fft = time.time()
    for i in range (loop):
        fft(a)  
    ffttime=time.time() - start_fft
    plt.title("DFT vs FFT Convergence", fontsize=22)
    plt.xlabel("N", fontsize=22)
    plt.ylabel("Time", fontsize=22)
    plt.xticks(fontsize=18) 
    plt.yticks(fontsize=18)
    plt.plot(N, dfttime, 'b*', markersize="10")
    plt.plot(N, ffttime, 'r*', markersize="10")
    plt.style.use("seaborn-dark")
    

for k in range(13):
    N=2**k
    L=2*np.pi
    dx=L/N
    x_pre= np.linspace(0, L, N+1)
    x=x_pre[0:N]    
    a=f(x)
    timer(N)


