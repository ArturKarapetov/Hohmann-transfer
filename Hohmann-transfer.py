#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
# Computational Exercise 3 for PHYS*2310 (Mechanics) - Winter 2020
# orbits.py
# Artur Karapetov - March 14th 2020

"""

# -*- coding: utf-8 -*-
"""
Modified from Kelly Roos, Engineering Physics | Bradley University

This exercise looks at the orbits of the terrestrial planets (Mercury, Venus,
Earth and Mars) using a modified Euler method - the "Euler-Cromer" method.

The same method will be used to launch a satellite into an elliptical orbit
that sends it from Earth to Mars.  This is known as a Hohmann transfer.

"""


#### Section 1:  Start by importing relevant libraries
#---------------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
from  matplotlib import animation as animation
from matplotlib.patches import Circle
import time as time

###Section 2:  Define functions for the script
#---------------------------------------------------------------------------

def streak(t,x,y,tl):
    DT = max(0,t-tl)
    pathE.set_data(x[DT:t+1],y[DT:t+1])
    frontE.set_data(x[t],y[t])
    return pathE,frontE,


### Section 3a:  Main body - Define Constants
#---------------------------------------------------------------------------
#FOR EARTH                              #FOR MARS 
# perihelion distance x = 0.983 AU      # perihelion distance x = 1.381 AU
# perihelion speed 0.0175 AU/day        # perihelion speed 0.0153 AU/day
# aphelion distance xmax = 1.01671 AU   # aphelion distance xmax = 1.66634 AU
# period 365 days                       # period 686.97 days

#FOR VENUS                             #FOR MERCURY 
# perihelion distance x = 0.72 AU      # perihelion distance x = 0.3075 AU
# perihelion speed 0.0205 AU/day       # perihelion speed 0.04031 AU/day
# aphelion distance xmax = 0.73 AU     # aphelion distance xmax = 0.466584 AU
# period 224.65 days                       # period 87.97 days

# Define physical constants
    # work in distance units AU (1 AU = 1.49598e11 m)
    # work in time units of days (1 day = 86400 s)

M_Sun = 1.9891e30#kg
mE = 5.972e24#kg 
G = (6.67408e-11)*(24*3600)**2/((1.49598e11)**3)               #G in  AU^3/d^2 units (SI: 6.67408e-11)
GM = G*M_Sun       #G*M_Sun in AU^3/d^2 units


# trajectory setup (time in days)
dt = 0.1
T = np.arange(0.0, 1095., dt)

# MERCURY SETUP
# 2D arrays for position and velocity (x,y components)
rE = np.zeros((T.size,2)) #array is created like this because dealing with a 2D orbit so need x & y coord.
vE = np.zeros_like(rE) #shape of rE which is T.size by 2.
# define initial conditions: (sun is at zero)
rE[0,:] = [0.3075, 0.0] #Planet starts at its perihelion distance on the x-axis
vE[0,:] = [0.0,0.04031] #speed is only in the 'y' direction since starting on the x-axis - perihelion speed

# VENUS SETUP
# 2D arrays for position and velocity (x,y components)
rE = np.zeros((T.size,2)) #array is created like this because dealing with a 2D orbit so need x & y coord.
vE = np.zeros_like(rE) #shape of rE which is T.size by 2.
# define initial conditions: (sun is at zero)
rE[0,:] = [0.72, 0.0] #Planet starts at its perihelion distance on the x-axis
vE[0,:] = [0.0,0.0205] #speed is only in the 'y' direction since starting on the x-axis - perihelion speed

# EARTH SETUP
# 2D arrays for position and velocity (x,y components)
rE = np.zeros((T.size,2)) #array is created like this because dealing with a 2D orbit so need x & y coord.
vE = np.zeros_like(rE) #shape of rE which is T.size by 2.
# define initial conditions: (sun is at zero)
rE[0,:] = [0.983, 0.0] #Planet starts at its perihelion distance on the x-axis
vE[0,:] = [0.0,0.0175] #speed is only in the 'y' direction since starting on the x-axis - perihelion speed

# MARS SETUP
# 2D arrays for position and velocity (x,y components)
rE = np.zeros((T.size,2)) #array is created like this because dealing with a 2D orbit so need x & y coord.
vE = np.zeros_like(rE) #shape of rE which is T.size by 2.
# define initial conditions: (sun is at zero)
rE[0,:] = [1.381, 0.0] #Planet starts at its perihelion distance on the x-axis
vE[0,:] = [0.0,0.0153] #speed is only in the 'y' direction since starting on the x-axis - perihelion speed
 

### Section 3b:  Main body - Set up Euler's (Euler-Cromer) Method
#---------------------------------------------------------------------------
# the Euler Method updates rE and vE using information from the previous time
# with Euler-Cromer, now rE is updated using the NEW values for vE
clockStart = time.time()

for t in range(1,T.size):
    #calculate acceleration at current position
    R = np.sqrt(np.sum(rE[t-1,:]**2))
    aE = -GM/(R**3)*rE[t-1]
    
    #update velocity using Euler
    vE[t,:] = vE[t-1,:] + aE*dt
    #vEC[t,:] = vE[t-1,:] + aE*dt
    
    # update rE using Euler
    rE[t,:] = rE[t-1,:] + vE[t,:]*dt  #Euler has vE[t-1,:], Euler-Cromer has vE[t,:]
    #rEC[t,:] = rE[t-1,:] + vE[t,:]*dt

    # Instead, update rE using Euler-Cromer

    ####need perihelion position for Mars, convert to AU and days - initial boost speed & when launched to reach Mars - need calculation for how close you get to Mars and calculate how close you are
    #Can ask Dr. Massa for help with part C for fun after I try to do it - it's too cool to ignore this time (unless I reeally need time )
    
print("duration of for loop: %5f" % (time.time()-clockStart) )


### Section 3c:  Main body - calculate the energy of Earth
#---------------------------------------------------------------------------
# We need to calculate the kinetic and potential energy of Earth 
# at each position on its orbit


### Section 3d:  Main body - plot the trajectory IN UNIFORM B FIELD
#---------------------------------------------------------------------------

runAnim = True
if not(runAnim):
    # plot the kinetic and potential energy of Earth over a few orbits
    fig = plt.figure(figsize=(24,12))
    ax1 = fig.add_subplot(1, 2, 1)
    ax2 = fig.add_subplot(1, 2, 2)
    


if runAnim:
    #plot trajectory in the xy-plane
    fig = plt.figure(figsize=(12,12))
    ax1 = fig.add_subplot(1, 1, 1)
    ax1.set_aspect('equal')
    ax1.set_xlim(-2,2)
    ax1.set_ylim(-2,2)
    ax1.add_artist(Circle((0.,0.),0.1,color='#ffaa00'))
    ax1.add_artist(Circle((0.03,0.03),0.01,color='#000000'))
    ax1.add_artist(Circle((-0.03,0.03),0.01,color='#000000'))
    ax1.add_artist(Circle((0.,-0.03),0.015,color='#000000'))
    
    # Set up Earth's orbit animation
    pathE, = ax1.plot(rE[0,0],rE[0,1],'k')
    frontE, = ax1.plot(rE[0,0],rE[0,1],'bo')

    streakLength = 1000     #length of trail following planet (shows orbit)
    aniSpeed = 10           #run-speed of animation; in ms
    ani = animation.FuncAnimation(fig,streak,
                                  fargs=(rE[::10,0],rE[::10,1],streakLength),
                                  interval=aniSpeed,blit=True)

