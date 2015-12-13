#Christopher Bell
#Leap Frog Method

from math import sin, cos, log, ceil, pi
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.family']= 'serif'
rcParams['font.size']= 16

def centraldiff (unm1, u, f, dt):
    """Returns the solution at n+1 using a central difference Eularian Method
    
    Parameters
    ----------
    unm1 : array of float
        solution at time step n-1
    u : array of float
        solution at time step n
    f : function
        function to compute RHS of SOE's
    dt : float
        time-increment
        
    Returns
    -------
    u_n_plus_1 : array of float
        solution at time step  n+1
    """
    
    return unm1 + 2*dt*f(u)

def f(u):
    """Returns RHS of the phugoid system of equations.
    
    Params
    ------
    u: array of float
       array containing the solution at time n
       
    Returns
    ------
    dudt: array of float
        array containing RHS from given u
    """
    
    #define locations of inputs in u matrix
    v = u[0]
    theta = u[1]
    x = u[2]
    y = u[3]
    
    return np.array([-g*sin(theta) - (C_D/C_L)*(g/(v_t**2))*v**2,\
    (-g*cos(theta))/v + (g/(v_t**2))*v, v*cos(theta), v*sin(theta)])

def rk2_step(u, f, dt):
    """Returns the solution at the next time step using 2nd order Runge-Kutta
    
    Parameter:
    ----------
    u : array of float
        solution at the previous time-step.
    f : function
        function to compute the RHS of the system of equations.
    dt : float
        time-increment.
    
    Returns:
    --------
    u_n_plus_1 : array of float
        solution at the next time step.
    """
    u_star = u + 0.5*dt*f(u)
    return u + dt*f(u_star)

def get_diffgrid(u_current, u_fine, dt):
    """Returns the difference between one grid and the fine grid using L-1 norm
    
    Params:
    ------
    u_current : array of float
        solution on the current grid.
    u_finest : array of float
        solution on the fine grid.
    dt : float
        time-increment on the current grid.
        
    Returns
    -------
    diffgrid : float
        difference computed in the L-1 norm
    """
    
    N_current = len(u_current[:,0])
    N_fine = len(u_fine[:,0])
    
    grid_size_ratio = ceil(N_fine/float(N_current))
    
    diffgrid = dt * np.sum(np.abs(u_current[:,2] - u_fine[::grid_size_ratio,2]))
    
    return diffgrid

#Global Params
g= 9.8 #m/s2
v_t = 4.9 #trim velocity m/s
C_D = 1.0/5.0 #Drag Coefficient
C_L = 1.0

#ICs
v0 = 6.5
theta0 = -0.1
x0 = 0.0
y0 = 2.0

#time increment
T = 15.0
dt = 0.01
N = int(T/dt)+1

u_centdiff = np.empty((N,4))
u_centdiff[0] = np.array([v0,theta0,x0,y0])

#first step
u_centdiff[1] = rk2_step(u_centdiff[0],f,dt)

for n in range(1,N-1):
    u_centdiff[n+1] = centraldiff(u_centdiff[n-1], u_centdiff[n], f, dt)

x_centdiff = u_centdiff[:,2]
y_centdiff = u_centdiff[:,3]

idx_neg_centdiff = np.where(y_centdiff<0.0)[0]

if len(idx_neg_centdiff)==0:
    idx_ground_centdiff = N-1
    print ('The glider has not reached the ground yet!')
else:
    idx_ground_centdiff = idx_neg_centdiff[0]
    
#plotting time
plt.figure
plt.grid(True)
plt.xlabel('$x$')
plt.ylabel('$y$')
plt.plot(x_centdiff[:idx_ground_centdiff], y_centdiff[:idx_ground_centdiff], color = 'k', ls='-',lw=2)
plt.title('Distance Traveled: {:.3f}'.format(x_centdiff[idx_ground_centdiff]), fontsize=18)
plt.show()

#convergence
r = 2
h = 0.001

dt_values = np.array([h, r*h, (r**2)*h])

u_values = np.empty_like(dt_values, dtype=np.ndarray)

for i,dt in enumerate(dt_values):
    t = np.linspace(0.0,T,N)
    u=np.empty((N,4))
    u[0]=np.array([v0,theta0,x0,y0])
    u[1] = rk2_step(u[0],f,dt)
    for n in range(1,N-1):
        u[n+1] = centraldiff(u[n-1],u[n],f,dt)
    
    u_values[i] = u
    
alpha = (log(get_diffgrid(u_values[2],u_values[1],dt_values[2])) -  log(get_diffgrid(u_values[1],u_values[0],dt_values[1]))) / log(r)

print('Order of convergence is alpha = {:.3f}'.format(alpha))
