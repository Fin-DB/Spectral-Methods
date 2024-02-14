# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 12:19:50 2023

@author: baile
"""

import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt

def diff_matrix(N,h):
    M1=np.zeros([N,N])
    for i in range(N-1):
        M1[i][i+1]=2/3
    M1[N-1][0]=2/3
    M2=np.zeros([N,N])
    for j in range(N-2):
        M1[j][j+2]=-1/12
    M2[N-2][0]=1/12
    M2[N-1][1]=1/12
    D=M1-M2
    return (D-np.transpose(D))/h

def fin():
    #Convergance of periodic finite difference methods 
    N = 2 ** np.arange(3, 13)
    for index in N:
        h = 2*np.pi/index
        x = np.arange(1, index+1) * h - np.pi
        u = np.exp(np.cos(x))
        uprime= -np.sin(x)*u   
        
        #Construct Finite Difference differentiation matrix 
        M_diff= diff_matrix(index,h)
        toetimesu = np.matmul(M_diff, u)
        Dif = toetimesu - uprime
        error = linalg.norm(Dif, np.inf)
        
        #Plotting the maximum error
        plt.xlabel("N", fontsize=22)
        plt.ylabel("Error ", fontsize=22)
        plt.xticks(fontsize=20) # Makes x labels digits larger
        plt.yticks(fontsize=20) # Makes y labels digits larger
        plt.title("Spectral vs Finite Differentiation",fontsize=22)
        plt.tight_layout(rect=[0, 0, 1, 1]) # Ensure the full figure is visible
        plt.loglog(index, error, "b*",markersize=10) # Plot data
        plt.style.use("seaborn-dark")
           
        
#Convergance of periodic spectral methods 
for index in range(4, 41, 2):
    h = 2*np.pi/index
    x = np.arange(1, index+1) * h - np.pi
    u = np.exp(np.cos(x))
    uprime= -np.sin(x)*u     
    
    #Construct Spectral differnetiation matrix 
    old_column = (0.5*(-1)**np.arange(1, index)*
    (np.cos(np.arange(1, index)*h/2)/np.sin(np.arange(1, index) * h/2)))
    column=np.insert(old_column, 0, 0)
    M_toeplitz = linalg.toeplitz(column, -column)
    toetimesu = np.matmul(M_toeplitz, u)
    Dif = toetimesu - uprime
    error = linalg.norm(Dif, np.inf)
    
    #Plotting the maximum error
    plt.loglog(index, error, "r*",markersize=10) # Plot data
    
fin()
#plt.savefig("Convergence of spectral differentiation vs finite.png",bbox_inches='tight')
plt.show()