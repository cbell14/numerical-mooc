#Christopher Bell
#1D Burgers' Equation
#02-04

###Animation issue

import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from matplotlib import rcParams, animation
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16
plt.ion()
from sympy import init_session
init_session()
from sympy.utilities.lambdify import lambdify


x, nu, t = sp.symbols('x nu t')
phi = sp.exp((-(x-4*t)**2)/(4*nu*(t+1))) + \
sp.exp(-((x-4*t-2*np.pi)**2)/(4*nu*(t+1)))
phi
sp.pretty_print (phi)

phiprime = phi.diff(x)
phiprime
sp.pretty_print (phiprime)

u = -2*nu*(phiprime/phi)+4
print(u)

ufunc = lambdify((t,x,nu),u)
print('The value of u at t=1, x=4, nu=3 is {}.'.format(ufunc(1,4,3)))

#Variable Declarations
nx = 101
nt = 100
dx = 2*np.pi/(nx-1)
nu = 0.07
dt = dx*nu

x = np.linspace(0,2*np.pi,nx)
un = np.empty(nx)
t = 0

u = np.asarray([ufunc(t,x0,nu) for x0 in x])
u

plt.figure(figsize=(8,5),dpi=100)
plt.plot(x,u,color='#003366',ls='--',lw=3)
plt.xlim([0,2*np.pi])
plt.ylim([0,10])
plt.show()

for n in range(nt):
    un = u.copy()
    
    u[1:-1] = un[1:-1] - un[1:-1] * dt/dx * (un[1:-1]-un[:-2]) + nu*dt/dx**2*\
    (un[2:] - 2*un[1:-1] + un[:-2])
    
    u[0] = un[0] - un[0] * dt/dx * (un[0] - un[-1]) + nu*dt/dx**2*\
    (un[1] - 2*un[0] + un[-1])
    
    u[-1] = un[-1] - un[-1]*dt/dx * (un[-1] - un[-2]) + nu*dt/dx**2*\
    (un[0] - 2*un[-1] + un[-2])
    
    u_analytical = np.asarray([ufunc(nt*dt,xi,nu) for xi in x])

    #plt.figure(figsize=(8,5), dpi=100)
    plt.plot(x,u, color='#003366',ls='--',lw=3, label='Computational')
    plt.plot(x, u_analytical, label='Analytical')
    plt.xlim([0,2*np.pi])
    plt.ylim([0,10])
    plt.legend()

u = np.asarray([ufunc(t,x0,nu) for x0 in x])

fig = plt.figure(figsize=(8,6))
ax = plt.axes(xlim=(0,2*np.pi),ylim=(0,10))
line = ax.plot([],[],color='#003366',ls='--',lw=3)[0]
line2 = ax.plot([],[],'k-',lw=2)[0]
ax.legend(['Computed','Analytical'])

def burgers(n):
    un = u.copy()
    
    u[1:-1] = un[1:-1] - un[1:-1] * dt/dx * (un[1:-1] - un[:-2]) + nu*dt/dx**2*\
                    (un[2:] - 2*un[1:-1] + un[:-2])

    u[0] = un[0] - un[0] * dt/dx * (un[0] - un[-1]) + nu*dt/dx**2*\
                (un[1] - 2*un[0] + un[-1])
    u[-1] = un[-1] - un[-1] * dt/dx * (un[-1] - un[-2]) + nu*dt/dx**2*\
                (un[0]- 2*un[-1] + un[-2])
        
    u_analytical = np.asarray([ufunc(n*dt, xi, nu) for xi in x])
    line.set_data(x,u)
    line2.set_data(x, u_analytical)
    

animation.FuncAnimation(fig, burgers,
                        frames=nt, interval=100)