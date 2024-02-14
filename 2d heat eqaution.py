# -*- coding: utf-8 -*-
"""
Created on Thu Jan  4 20:30:54 2024

@author: baile
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D

def dim2heat(u, wavenumbersx, wavenumbersy):
    #Using FFT for derivatives and forward Euler for time stepping
    u_hat=np.fft.fft2(u)
    v_hatx= u_hat * -(wavenumbersx)**2
    v_haty= u_hat * -(wavenumbersy[:, np.newaxis])**2 #Unsure
    u_derivx=np.fft.ifft2(v_hatx).real
    u_derivy=np.fft.ifft2(v_haty).real
    u_new = u + dt*(u_derivx+u_derivy)
    return u_new

        
#Heat Equation conditions
Lx=10
Ly=10
Nx=128
Ny=128
dx=Lx/Nx
dy=Ly/Ny
x=np.arange(0,Lx,dx)
y=np.arange(0,Ly,dy)
xs, ys = np.meshgrid(x, y)
u=np.zeros([Nx,Ny])
u[int(Nx/4):int(3*Nx/4),int(Ny/4):int(3*Ny/4)]=1

t_iteration = 2000
t=0
dt = 0.0001
wavenumbersx = np.concatenate((np.arange(0, Nx//2), [0], np.arange(-Nx//2+1, 0)))
wavenumbersy = np.concatenate((np.arange(0, Ny//2), [0], np.arange(-Ny//2+1, 0)))

for i in range(t_iteration): 
    #Plotting our solution for u and the analytic solution
    if i%20==0:
        plt.figure()
        plt.imshow(np.flipud(u),aspect=1)
        plt.axis('off')
        plt.set_cmap('jet') #flag cool jet normal
        plt.savefig("{} 2D heat equation .png".format(i))
        plt.show()
        # fig = plt.figure()
        # ax = Axes3D(fig)
        # plt.tight_layout(rect=[0, 0, 1, 1])
        # ax.plot_surface(xs, ys, u, rstride=1, cstride=1, cmap='jet')
        plt.savefig("2D heat equation {}.png".format(t))
        # plt.show()
    u=dim2heat(u,wavenumbersx,wavenumbersy)
    t=t+dt