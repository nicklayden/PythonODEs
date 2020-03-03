# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 23:26:10 2020

@author: nickl
"""

# -*- coding: utf-8 -*-
"""
Basic method for solving ODEs numerically and plotting them in python.

"""

from mpl_toolkits import mplot3d
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def circle(x):
    return 1 - x**2

def sphere(x,y):
    return 1 - x**2 - y**2
###############################################################################
#
# Define the systems of ODEs to solve
#
###############################################################################
def scalarfield_cosm(vec,t,lam):
    """Single scalar field + Cosmological constant problem. 3-D ODE """
    x, y, z = vec
    # Deceleration parameter
    q = 2*x**2 - y**2 - z**2
    
    f = [x*(q-2) - (np.sqrt(6)/2.)*lam*y**2 , 
         y*((np.sqrt(6)/2.)*lam*x + q + 1),
         z*(q+1)]
    return f


###############################################################################
#
# This function is what we call on the systems to generate numerical solutions
#
###############################################################################
def solve(system, ic, params):
    # tb is our backwards integration time series
    # tf is our forwards integration time series
    # for autonomous systems, the specific time values do not matter.
    tb = np.linspace(0, -50, 5001)
    tf = np.linspace(0,  50, 5001)

    
    # For each initial condition, evolve the solution backwards to a source
    # and forwards to a sink.
    solb = odeint(system, ic, tb, args=params)
    solf = odeint(system, ic, tf, args=params)
    

    
    # 3D plot of the phase plane
    ax.plot3D(solb[:,0],solb[:,1],solb[:,2], c="b")
    ax.plot3D(solf[:,0],solf[:,1],solf[:,2], c="b")

    



###############################################################################
#
# This is just standard practice in python to run a file
#
###############################################################################
if __name__ == "__main__":
    # Define a shared figure that all solutions will be plotted on
    fig = plt.figure()
    ax = plt.axes(projection="3d")
    
    # Create initial conditions inside the sphere of radius 1. (Actually in an
    # octant of the sphere, solution is mirrored everywhere else!)
    N_ic = 25
    height = 0.1
    width = 0.1
    x_init = np.linspace(-sphere(height,width),sphere(height,width),N_ic)
    y_init = [height for i in range(N_ic)]
    z_init = [width for i in range(N_ic)]
    
    # Define extra parameters to the problem - i.e. like lambda in cosmo 
    # problems.
    # This has to be a 'tuple' in python, a tuple is a list like (a,b,c,d) with 
    # circular brackets.
    # Since I only have 1 parameter, I needed to add a comma at the end to make 
    # python think it was
    # actually a tuple.
    params = (6,)
    

    # Solve for a whole set of initial conditions
    # Note that each call to solve() will add the solution to the axes in fig
    for i in range(len(z_init)):
        z0 = [x_init[i],y_init[i],z_init[i]]
        solve(scalarfield_cosm,z0,params)























