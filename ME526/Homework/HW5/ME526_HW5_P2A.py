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
    that the initial temperature is known throughout.  Function in 1D."""
    dx, dt, u = parameters
    return point - (0.5*u*dt/dx) * (nextpoint - prevpoint)


def secondscheme(nextpoint, prevpoint, parameters, oldpoint):
    """A function that finds the solution at the next time step for a given
    point using the leapfrog approximation for time advancement and second
    order central difference approximation for spacial derivative.  Assumes
    that the initial temperature is known throughout.  Function in 1D."""
    dx, dt, u = parameters
    return oldpoint - (u*dt/dx) * (nextpoint - prevpoint)


def laxwendroff(nextpoint, point, prevpoint, parameters):
    """A function given in the homework to solve the differential equation."""
    dx, dt, u = parameters
    gamma = u * dt / dx
    return point - gamma*0.5*(nextpoint - prevpoint) + 0.5*(gamma**2) *\
        (nextpoint - 2*point + prevpoint)


import numpy as np
import matplotlib.pyplot as plt

xmin, xmax, xmeshnumber = 0, 1, 50
dx = (xmax - xmin) / xmeshnumber
tmin, tmax, tmeshnumber = 0, 8, 1600
dt = (tmax - tmin) / tmeshnumber

u = 0.08

xmesh = np.arange(xmin, (xmax + dx), dx)
tmesh = np.arange(tmin, (tmax + dt), dt)
tempmesh = np.zeros(len(xmesh))

evaltimes = [0, 4, 8]

# Get the exact solution for all of the times that we're evaluating.
exacttempmesh = np.zeros(len(xmesh))
exacteval = []
for time in evaltimes:
    for x in range(len(exacttempmesh)):
        exacttempmesh[x] = exactsolution(xmesh[x], time, u)
    exacteval.append(exacttempmesh)
    exacttempmesh = np.zeros(len(xmesh))

# Set up the initial mesh for the numerical solution
for x in range(len(tempmesh)):
    if xmesh[x] < 0.2:
        tempmesh[x] = 1 - (10*xmesh[x] - 1)**2

# Solve the equation using the first scheme
firstsolution = []
parameters = dx, dt, u
firstmesh = np.array(tempmesh)
for t in range(len(tmesh)):
    oldmesh = np.array(firstmesh)
    for x in range(1, len(xmesh) - 1):
        firstmesh[x] = firstscheme(oldmesh[x+1], oldmesh[x], oldmesh[x-1],
                                  parameters)
    firstmesh = np.array(firstmesh)
    firstsolution.append(firstmesh)

# Solve the equation using the second scheme
secondsolution = []
secondmesh, oldermesh, halfmesh = np.array(tempmesh), np.array(tempmesh),\
    np.array(tempmesh)
# Do RK2 for one timestep
for x in range(1, len(tempmesh) - 1):
    halfmesh[x] = tempmesh[x] - (0.25*u*dt/dx)*(tempmesh[x+1] - tempmesh[x-1])
for x in range(1, len(halfmesh) - 1):
    oldermesh[x] = tempmesh[x] - (0.5*u*dt/dx)*(halfmesh[x+1] - halfmesh[x-1])
secondsolution.append(oldermesh)
# Use the function above to solve for the rest of the range
for t in range(1, len(tmesh)):
    oldmesh = np.array(secondmesh)
    for x in range(1, len(xmesh) - 1):
        secondmesh[x] = secondscheme(oldmesh[x+1], oldmesh[x-1], parameters,
                                     oldermesh[x])
    secondmesh = np.array(secondmesh)
    oldermesh = np.array(oldmesh)
    secondsolution.append(secondmesh)

# Solve the equation using the Lax-Wendroff scheme
lwsolution = []
lwmesh = np.array(tempmesh)
for t in range(len(tmesh)):
    oldmesh = np.array(lwmesh)
    for x in range(1, len(xmesh) - 1):
        lwmesh[x] = laxwendroff(oldmesh[x+1], oldmesh[x], oldmesh[x-1],
                                  parameters)
    lwmesh = np.array(lwmesh)
    lwsolution.append(lwmesh)

# Plot the exact solution
for i in range(len(evaltimes)):
    plt.title('Exact Soltuion (dt = {})'.format(dt), fontsize=16)
    plt.xlabel('x Value', fontsize=12)
    plt.ylabel('Temperature Value', fontsize=12)
    plt.plot(xmesh, exacteval[i], label='Time = {}'.format(evaltimes[i]))
    plt.legend(bbox_to_anchor=(1, 1), loc=2)
    plt.grid(b=True, which='both')
plt.show()

# Plot the solution obtained from the first scheme
for i in [0, int((tmeshnumber / 2)), tmeshnumber]:
    plt.title('First Scheme Soltuion (dt = {})'.format(dt), fontsize=16)
    plt.xlabel('x Value', fontsize=12)
    plt.ylabel('Temperature Value', fontsize=12)
    plt.plot(xmesh, firstsolution[i], label='Time = {}'.format(tmesh[i]))
    plt.legend(bbox_to_anchor=(1, 1), loc=2)
    plt.grid(b=True, which='both')
    plt.ylim(-1, 1.5)
plt.show()

# Plot the solution obtained from the second scheme
for i in [0, int((tmeshnumber / 2)), tmeshnumber]:
    plt.title('Second Scheme Soltuion (dt = {})'.format(dt), fontsize=16)
    plt.xlabel('x Value', fontsize=12)
    plt.ylabel('Temperature Value', fontsize=12)
    plt.plot(xmesh, secondsolution[i], label='Time = {}'.format(tmesh[i]))
    plt.legend(bbox_to_anchor=(1, 1), loc=2)
    plt.grid(b=True, which='both')
    plt.ylim(-1, 1.5)
plt.show()

# Plot the soltuion obtained from the Lax-Wendroff scheme
for i in [0, int((tmeshnumber / 2)), tmeshnumber]:
    plt.title('Lax-Wendroff Scheme Soltuion (dt = {})'.format(dt), fontsize=16)
    plt.xlabel('x Value', fontsize=12)
    plt.ylabel('Temperature Value', fontsize=12)
    plt.plot(xmesh, lwsolution[i], label='Time = {}'.format(tmesh[i]))
    plt.legend(bbox_to_anchor=(1, 1), loc=2)
    plt.grid(b=True, which='both')
plt.show()
