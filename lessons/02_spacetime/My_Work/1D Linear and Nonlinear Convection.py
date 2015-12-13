#Christopher Bell
#Module 2
#1D Linear and nonlinear convection

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16

###Linear

nx = 1000
sigma = 0.5
dx = 2.0/(nx-1)
nt = 125
c = 1.0
dt = (dx*sigma)/c


u = np.ones(nx)
u[.5/dx : 1/dx+1] = 2 #setting u = 2 between 0.5 and 1 as per IC
#print(u)

plt.plot(np.linspace(0,2,nx), u, color='#003366', ls='--',lw=3)
plt.title('Linear Convection')
plt.ylim(0,2.5)
plt.show()

un = np.ones(nx)

#inefficient can make x positions an array avoids double loop
for n in range(nt):
    un = u.copy()
    for i in range(1,nx):
        u[i] = un[i] - c*dt/dx*(un[i]-un[i-1])
        
plt.plot(np.linspace(0,2,nx), u, color='#003366', ls='--',lw=3)
plt.ylim(0,2.5)
print nt*dt
plt.show()

####Non-Linear
nx_nl = 1000
dx_nl = 2.0/(nx_nl-1)
nt_nl = 100
sigma = 0.5

u_nl = np.ones(nx_nl)
u_nl[0.5/dx_nl : 1/dx_nl+1] = 2

un_nl = np.ones(nx_nl)

x_nl = np.linspace(0,2,nx_nl)

plt.figure()
plt.plot(np.linspace(0,2,nx_nl), u_nl, color='#003366', ls='--', lw=3)
plt.ylim(0,2.5)
plt.show()

#for q in range(nt_nl):
#    un_nl = u_nl.copy()
#    for k in range(1,nx_nl):
#        u_nl[k] = un_nl[k] - un_nl[k]*dt_nl/dx_nl*(un_nl[k]-un_nl[k-1])

#faster way

time = 0.0
count = 0

while time<=0.5:
    count+=1
    dt_nl = (sigma*dx_nl)/max(u_nl) #makes CFL condition dependent on maximum of u
    print dt_nl
    time = time + dt_nl
    un_nl = u_nl.copy()
    u_nl[1:] = un_nl[1:]-un_nl[1:]*dt_nl/dx_nl*(un_nl[1:]-un_nl[0:-1]) #backward difference
    #u_nl[1:] = un_nl[1:]-un_nl[1:]*dt_nl/dx_nl*(-un_nl[0:-1]+un_nl[1:]) #forward
    u_nl[0] = 1.0 
            
print count
plt.plot(np.linspace(0,2,nx_nl), u_nl, color='#003366', ls='--', lw=3)
plt.ylim(0,2.5)
plt.title('Non-Linear Convection')
plt.show()