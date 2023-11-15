#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 14:23:27 2023

@author: jason
"""

import numpy as np
from matplotlib import pyplot as plt

x = np.linspace(0,50e3,1001) # longitudinal coordinate [m]

aB = 3000 # height of gaussian peak 
bB = -20000 # location of peak of gaussian [m]
cB = 20000 # width of gaussian

zmin = -500 # max fjord depth

bedrock = aB*np.exp(-(x-bB)**2/(2*cB**2)) 

aS = 300 # height of sill
bS = 40000 # location of sill
cS = 500 # width of sill

sill = aS*np.exp(-(x-bS)**2/(2*cS**2))
 
zs = bedrock + sill + zmin

plt.plot(x,zs)