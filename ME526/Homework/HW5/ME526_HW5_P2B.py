#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 10:28:32 2016

@author: Andrew Alferman and Dan Magee
"""


def exactsolution(x, t, u):
    """A function that returns the exact solution for a given t and x."""
    if 0 <= (x - u*t) and (x - u*t) <= 0.2:
        temp = 1 - (10 * (x - u*t) -1)**2
    else:
        temp = 0
    return temp


def firstscheme(nextpoint, point, prevpoint, parameters, alpha):
    """A function that finds the solution at the next time step for a given
    point using the forward euler approximation for time advancement and second
    order central difference approximation for spacial derivative.  Assumes
    that the initial temperature is known throughout.  Function in 1D."""
    dx, dt, u = parameters
    return point - (0.5*u*dt/dx) * (nextpoint - prevpoint) +\
        (alpha*dt/(dx**2)) * (nextpoint - 2*point + prevpoint)


def secondscheme(nextpoint, point, prevpoint, parameters, oldpoint, alpha):
    """A function that finds the solution at the next time step for a given
    point using the leapfrog approximation for time advancement and second
    order central difference approximation for spacial derivative.  Assumes
    that the initial temperature is known throughout.  Function in 1D."""
    dx, dt, u = parameters
    return oldpoint - (u*dt/dx) * (nextpoint - prevpoint) +\
        (alpha*dt/(dx**2)) * (nextpoint - 2*point + prevpoint)


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
tmin, tmax, tmeshnumber = 0, 8, 32
dt = (tmax - tmin) / tmeshnumber

u = 0.08
alpha = 0.001

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
                                  parameters, alpha)
    firstmesh = np.array(firstmesh)
    firstsolution.append(firstmesh)

# Solve the equation using the second scheme
secondsolution = []
secondmesh, oldermesh, halfmesh = np.array(tempmesh), np.array(tempmesh),\
    np.array(tempmesh)
# Do RK2 for one timestep
for x in range(1, len(tempmesh) - 1):
    halfmesh[x] = tempmesh[x] + dt * 0.5 * (
        (alpha*dt/(dx**2))*(tempmesh[x+1] - 2*tempmesh[x] + tempmesh[x-1]) -\
        (0.5*u*dt/dx)*(tempmesh[x+1] - tempmesh[x-1]))
for x in range(1, len(halfmesh) - 1):
    oldermesh[x] = tempmesh[x] + dt * (
        (alpha*dt/(dx**2))*(halfmesh[x+1] - 2*halfmesh[x] + halfmesh[x-1]) -\
        (0.5*u*dt/dx)*(halfmesh[x+1] - halfmesh[x-1]))
secondsolution.append(oldermesh)
# Use the function above to solve for the rest of the range
for t in range(1, len(tmesh)):
    oldmesh = np.array(secondmesh)
    for x in range(1, len(xmesh) - 1):
        secondmesh[x] = secondscheme(oldmesh[x+1], oldmesh[x], oldmesh[x-1],
                                     parameters, oldermesh[x], alpha)
    secondmesh = np.array(secondmesh)
    oldermesh = np.array(oldmesh)
    secondsolution.append(secondmesh)

# Solve the equation using the Lax-Wendroff scheme
lwsolution = []
lwmesh = np.array(tempmesh)
errorlist = []
for t in range(len(tmesh)):
    oldmesh = np.array(lwmesh)
    error = [0]
    for x in range(1, len(xmesh) - 1):
        lwmesh[x] = laxwendroff(oldmesh[x+1], oldmesh[x], oldmesh[x-1],
                                  parameters)
        if tmesh[t] == 0. or tmesh[t] == 4. or tmesh[t] == 8.:
            appenderror = True
            error.append(np.abs(lwmesh[x] -
                                exactsolution(xmesh[x], tmesh[t], u)))
    if appenderror == True:
        appenderror = False
        error.append(0)
        errorlist.append(error)
    lwmesh = np.array(lwmesh)
    lwsolution.append(lwmesh)

# Plot the exact solution
for i in range(len(evaltimes)):
    plt.title('Exact Solution', fontsize=16)
    plt.xlabel('x Value', fontsize=12)
    plt.ylabel('Temperature Value', fontsize=12)
    plt.plot(xmesh, exacteval[i], label='Time = {}'.format(evaltimes[i]))
    plt.legend(bbox_to_anchor=(1, 1), loc=2)
    plt.grid(b=True, which='both')
plt.show()

# Plot the solution obtained from the first scheme
for i in [0, int((tmeshnumber / 2)), tmeshnumber]:
    plt.title('First Scheme Solution (dt = {})'.format(dt), fontsize=16)
    plt.xlabel('x Value', fontsize=12)
    plt.ylabel('Temperature Value', fontsize=12)
    plt.plot(xmesh, firstsolution[i], label='Time = {}'.format(tmesh[i]))
    plt.legend(bbox_to_anchor=(1, 1), loc=2)
    plt.grid(b=True, which='both')
plt.show()

# Plot the solution obtained from the second scheme
for i in [0, int((tmeshnumber / 2)), tmeshnumber]:
    plt.title('Second Scheme Solution (dt = {})'.format(dt), fontsize=16)
    plt.xlabel('x Value', fontsize=12)
    plt.ylabel('Temperature Value', fontsize=12)
    plt.plot(xmesh, secondsolution[i], label='Time = {}'.format(tmesh[i]))
    plt.legend(bbox_to_anchor=(1, 1), loc=2)
    plt.grid(b=True, which='both')
plt.show()

# Plot the soltuion obtained from the Lax-Wendroff scheme
for i in [0, int((tmeshnumber / 2)), tmeshnumber]:
    plt.title('Lax-Wendroff Scheme Solution (dt = {})'.format(dt), fontsize=16)
    plt.xlabel('x Value', fontsize=12)
    plt.ylabel('Temperature Value', fontsize=12)
    plt.plot(xmesh, lwsolution[i], label='Time = {}'.format(tmesh[i]))
    plt.legend(bbox_to_anchor=(1, 1), loc=2)
    plt.grid(b=True, which='both')
plt.show()

# Plot the error of the Lax-Wendroff scheme
for i in errorlist:
    plt.title('Error of Lax-Wendroff Scheme (dt = {:.3f})'.format(dt), fontsize=16)
    plt.xlabel('x Value', fontsize=12)
    plt.ylabel('Absolute Error Value', fontsize=12)
    plt.plot(xmesh, i)
    plt.legend(bbox_to_anchor=(1, 1), loc=2)
    plt.grid(b=True, which='both')
plt.show()