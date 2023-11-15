#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 11:22:31 2023

@author: jason
"""

import numpy as np
from matplotlib import pyplot as plt

#%%

dBdz = 0.01
ELA = 1000
Bmax = 5

k = 0.005 # smoothing factor


z = np.linspace(0,3000,1001)

zthreshold = Bmax/dBdz+ELA

heaviside = dBdz * (1- np.exp(2*k*(z-zthreshold))/(1+np.exp(2*k*(z-zthreshold))))

B = dBdz * (z - 1/(2*k)*np.log(1 + np.exp(2*k*(z-zthreshold)))) - 10

plt.plot(z, dBdz*z - 10 )
plt.plot(z, B, '.')
