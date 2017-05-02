#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 10:28:32 2016

@author: andrewalferman
"""


def exactsolution(x, t, u):
    """A function that returns the exact solution for a given t and x."""
    if 0 <= (x - u*t) and (x - u*t) <= 0.2:
        temp = 1 - (10 * (x - u*t) -1)**2
    else:
        temp = 0
    return temp


def firstscheme(nextpoint, point, prevpoint, parameters):
    """A function that finds the solution at the next time step for a given
    point using the forward euler approximation for time advancement and second
    order central difference approximation for spacial derivative.  Assumes
    that the initial temperature is know throughout.  Function in 1D."""
    dx, dt, u = parameters
    return point - u*dt*(nextpoint - 2*point + point)/(dx**2)


import numpy as np
import matplotlib.pyplot as plt

xmin, xmax, xmeshnumber = 0, 1, 50
dx = (xmax - xmin) / xmeshnumber
tmin, tmax, tmeshnumber = 0, 8, 80
dt = (tmax - tmin) / tmeshnumber

u = 0.08

xmesh = np.arange(xmin, (xmax + dx), dx)
tmesh = np.arange(tmin, (tmax + dt), dt)
tempmesh = np.zeros(len(xmesh))

evaltimes = [0, 4, 8]

# Get the exact solution for all of the times that we're evaluating.
exacttempmesh = np.zeros(len(xmesh))
exacttempmesh = np.array(exacttempmesh)
exacteval = []
for time in evaltimes:
    for x in range(len(exacttempmesh)):
        exacttempmesh[x] = exactsolution(xmesh[x], time, u)
    exacteval.append(exacttempmesh)
    exacttempmesh = np.zeros(len(xmesh))
    print()

# Set up the initial mesh for the numerical solution
for x in range(len(tempmesh)):
    if xmesh[x] < 0.2:
        tempmesh[x] = 1 - (10*xmesh[x] - 1)**2

# Solve the equation using the first scheme
firstsolution = []
parameters = dx, dt, u
for t in range(len(tmesh)):
    oldmesh = tempmesh
    for x in range(1, len(xmesh) - 1):
        tempmesh[x] = firstscheme(tempmesh[x+1], tempmesh[x], tempmesh[x-1],
                                  parameters)
    firstsolution.append(tempmesh)

# Plot the exact solution
for i in range(len(evaltimes)):
    plt.title('Exact Soltuion', fontsize=16)
    plt.xlabel('x Value', fontsize=12)
    plt.ylabel('Temperature Value', fontsize=12)
    plt.plot(xmesh, exacteval[i], label='Time = {}'.format(evaltimes[i]))
    plt.legend(bbox_to_anchor=(1, 1), loc=2)
    plt.grid(b=True, which='both')
plt.show()