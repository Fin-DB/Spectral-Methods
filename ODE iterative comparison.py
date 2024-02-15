# -*- coding: utf-8 -*-
"""
Created on Sat Jan  6 17:03:50 2024

@author: baile
"""

import numpy as np
import matplotlib.pyplot as plt

def f(y,t):
    return np.exp(t)

def euler(y,t,h):
    g=f(y,t)
    return y + h*g

def leapfrog(y,y_old,t,h):
    g=f(y,t)
    return 2*h*g + y_old

def RK4(y,t,h):
    k1=f(y,t)
    k2=f(y+h*k1/2,t+h/2)
    k3=f(y+h*k2/2,t+h/2)
    k4=f(y+h*k3,t+h)
    return y + h/6 * (k1+2*k2+2*k3+k4)

final_time=5
Nt=50
h=final_time/Nt
t=np.linspace(0, final_time,Nt)
y=np.zeros([3,Nt])
h=final_time/Nt
y[0][0]=1
y[1][0]=1
y[2][0]=1
y[0][1]=np.e
y[1][1]=np.e
y[2][1]=np.e

x=np.linspace(0,final_time,64)
for i in range(1,len(t)-1):
    y[0][i+1]=euler(y[0][i],t[i],h)
    print(y[1][i-1])
    y[1][i+1]=leapfrog(y[1][i],y[1][i-1],t[i],h) 
    y[2][i+1]=RK4(y[2][i],t[i],h)

plt.xlabel("Time (t)", fontsize=22)
plt.ylabel("y(t)", fontsize=22)
plt.xticks(fontsize=20) # Makes x labels digits larger
plt.yticks(fontsize=20) # Makes y labels digits larger
plt.title("Accuracy of Different Iterative Methods h={}".format(h),fontsize=22)
plt.style.use("seaborn")
plt.tight_layout(rect=[0, 0, 1, 1])
plt.plot(t,y[0], "b",label='Euler')
plt.plot(t,y[1], "r",label='Leapfrog')
plt.plot(t,y[2], "g",label='RK4')
plt.plot(x,np.exp(x), "m",label='Analytical')
plt.legend(loc='upper left', prop={'size': 14})
plt.savefig("Accuracy of Different Iterative Methods h={}.png".format(h))
plt.show()