# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 11:17:28 2022

@author: fkurt
"""


import importlib.machinery

pyMCDS = importlib.machinery.SourceFileLoader('pyMCDS','./analysis/pyMCDS.py').load_module()

import os.path
from os import path

from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.animation

import numpy as np
import pandas as pd
#from fury import window, actor, utils, primitive, io, ui
#from fury.data import read_viz_textures, fetch_viz_textures
import itertools
#import vtk
import glob
import time
import random
import scipy.io as sio
import xml.etree.ElementTree as ET
from mpl_toolkits.axes_grid1.inset_locator import inset_axes



saving_times = np.arange(0,2040,60)

main_path = Path(os.getcwd()).parent
out_path = os.path.join(main_path, "output")

os.chdir(out_path)


time_point = "output000000"
number_of_frames = len(saving_times)

biomass_fluxes = []
int_glc = []
int_gln = []
int_lac = []
time=[]

for i in range(number_of_frames):
        time_p = time_point + '%02d'%(i)
        if path.exists(time_p + "_cells.mat"):
            data = sio.loadmat(time_p + "_cells.mat")
            intracellular_glucose = data['cells'][107]
            intracellular_glutamine = data['cells'][108]
            intracellular_lactate = data['cells'][109]
            biomass_flux = data['cells'][98]
            biomass_fluxes.append(biomass_flux[0])
            int_glc.append(intracellular_glucose[0])
            int_gln.append(intracellular_glutamine[0])
            int_lac.append(intracellular_lactate[0])
            time.append(i)


plt.figure()
plt.plot(time[1:],biomass_fluxes[1:])
plt.title(' Biomass Fluxes over time')
plt.xlabel('time(hr)')
plt.ylabel('Biomass Flux (1/hr)')



plt.figure()
plt.plot(time,int_glc)
plt.title(' Intracellular Glucose Concentration')
plt.xlabel('time(hr)')
plt.ylabel('Intracellular Glucose Concentration (1/mM)')


plt.figure()
plt.plot(time,int_gln)
plt.title(' Intracellular Glutamine Concentration')
plt.xlabel('time(hr)')
plt.ylabel('Intracellular Glutamine Concentration (1/mM)')


plt.figure()
plt.plot(time,int_lac)
plt.title(' Intracellular Lactate Concentration')
plt.xlabel('time(hr)')
plt.ylabel('Intracellular Lactate Concentration (1/mM)')


