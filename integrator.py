# -*- coding: utf-8 -*-
"""
Basic method for solving ODEs numerically and plotting them in python.

"""

import numpy as np
from scipy.integrate import odeint
import scipy.integrate as sp
import matplotlib.pyplot as plt

def circle(x):
    return 1 - x**2

###############################################################################
#
# Define the systems of ODEs to solve
#
###############################################################################
def system(z, t, a):
    """A damped harmonic oscillator."""
    return [z[1] - a*z[0], -z[0]]

def scalarfield(z,t,lam):
    """Single scalar field cosmology problem. 2-D ODE """
    x, y = z
    f = [x*(-2 + 2*x**2-y**2) + 0.5*np.sqrt(6)*y**2, 
         y*(1 + 2*x**2 - y**2 - 0.5*np.sqrt(6)*lam*x)]
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
    
    
    # Plots the time evolution of each variable
#    plt.plot(tb, solb)
#    plt.plot(tf,solf)
    
    
    # Plots the phase plane, i.e. no explicit time in these plots
    plt.plot(solb[:,0],solb[:,1],label="Backward Solution",c='b')
    plt.plot(solf[:,0],solf[:,1],label="Forward Solution",c='r')
    plt.grid(True)
    plt.xlim(-1,1)
    plt.ylim(0,1.2)
    plt.show()
    



###############################################################################
#
# This is just standard practice in python to run a file
#
###############################################################################
if __name__ == "__main__":
    # Define Initial Conditions
    z0 = [1,1]
    
    # Create initial conditions inside the circle.
    eps = 1e-3
    N_ic = 25
    height = 0.1
    z_init = np.linspace(-circle(height),circle(height),N_ic)
    y_init = [height for i in range(N_ic)]
    
    
    # Define extra parameters to the problem - i.e. like lambda in cosmo 
    # problems.
    # This has to be a 'tuple' in python, a tuple is a list like (a,b,c,d) with 
    # circular brackets.
    # Since I only have 1 parameter, I needed to add a comma at the end to make 
    # python think it was
    # actually a tuple.
    params = (np.sqrt(3)/2,)
    
    # This function takes in a system, initial conditions, and parameters
    # and solves them forwards and backwards for the time domain specified 
    # inside the function above.
#    solve(system,z0,params)
    
    # Solve for a whole set of initial conditions
    for i in range(len(z_init)):
        z0 = [z_init[i],y_init[i]]
        solve(scalarfield,z0,params)























