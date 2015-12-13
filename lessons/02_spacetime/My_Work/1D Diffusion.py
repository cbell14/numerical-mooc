#Christopher Bell
#1-D Diffusion
#Notebook 02-03

###Cant import animation

import numpy as np
import matplotlib.pyplot as plt
plt.ion()
#from JSAnimation import display_animation
from matplotlib import rcParams, animation
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16

nx =41
dx = 2.0/(nx-1)
nt = 50
nu = 0.3 #viscosity
sigma = 0.2
dt = (sigma*(dx**2))/nu

x = np.linspace(0,2,nx)

u = np.ones(nx)
u[.5/dx : 1/dx+1] = 2

un = np.ones(nx)

for n in range(nt):
    un=u.copy()
    u[1:-1] = un[1:-1] + ((nu*dt)/(dx**2))*(un[2:] - 2*un[1:-1] + un[0:-2])    
    plt.plot(np.linspace(0,2,nx), u, color='#003366', ls='--',lw=3)
    plt.ylim(0,2.5)
    


#fig = plt.figure(figsize=(8,5))
#ax = plt.axes(xlim=(0,2), ylim=(1,2.5))
#line = ax.plot([],[],color='#003366',ls='--',lw=3)[0]

#def diffusion(i):
#    line.set_data(x,u)
    
#    un = u.copy()
#    u[1:-1] = un[1:-1] + ((nu*dt)/(dx**2))*(un[2:] - 2*un[1:-1] + un[0:-2])