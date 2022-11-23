# -*- coding: utf-8 -*-
"""
Created on Fri Aug 27 14:41:55 2021

@author: Furkan and Kali
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

saving_times = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 20.0])


main_path = Path(os.getcwd()).parent
GT_out_path = os.path.join(main_path, "output_Ground_Truth")


os.chdir(GT_out_path)

time_point = "output000000"
number_of_frames = len(saving_times)

def get_subs_name ():
        tree = ET.parse("initial.xml")
        root = tree.getroot()
        
        subs_names = []
        for substrate in root.iter('variable'):
            subs_names.append(substrate.attrib['name'])
        
        return subs_names

subs_list = get_subs_name()


def GT_data_parser (time_point):
    fine_data = sio.loadmat(time_point + "_microenvironment0.mat")['multiscale_microenvironment']
    fine_x = np.unique(fine_data[0,:])
    fine_y = np.unique(fine_data[1,:])
    fine_X, fine_Y = np.meshgrid(fine_y, fine_x)
    fine_oxy = fine_data[4,np.where(fine_data[2,:] == 16)]
    fine_oxy = fine_oxy.reshape((len(fine_y),len(fine_x)))
    fine_oxy = np.transpose(fine_oxy)
    fine_glu = fine_data[5,np.where(fine_data[2,:] == 16)]
    fine_glu = fine_glu.reshape((len(fine_y),len(fine_x)))
    fine_glu = np.transpose(fine_glu)
    fine_chem = fine_data[6,np.where(fine_data[2,:] == 16)]
    fine_chem = fine_chem.reshape((len(fine_y),len(fine_x)))
    fine_chem = np.transpose(fine_chem)
    fine_oxy_tuple = (fine_X, fine_Y, fine_oxy)
    fine_glu_tuple = (fine_X, fine_Y, fine_glu)
    fine_chem_tuple = (fine_X, fine_Y, fine_chem)
    fine_tuple = (fine_oxy_tuple, fine_glu_tuple, fine_chem_tuple)
    return fine_tuple
    
fig, axs = plt.subplots()

# color bar
tp = "output00000020"
ft = GT_data_parser(tp)
fine_X, fine_Y, fine_oxy = ft[0]
w_X = fine_X
w_Y = fine_Y
w_O = fine_oxy
zmin = min([min(zl) for zl in w_O])
zmax = max([max(zl) for zl in w_O])
#levels = np.linspace(zmin, 0.28500001,41)
#kw = dict(levels=levels, vmin=zmin, vmax=0.28500001, origin='lower')
levels = np.linspace(0, 17.5,41)
kw = dict(levels=levels, vmin=0, vmax=17.5, origin='lower')
cp = axs.contourf(w_Y,w_X,w_O, **kw)
cbar = plt.colorbar(cp,format='%0.4f')
axs.clear()


def animate(i):
    time_p= time_point + '%02d'%(i)
    ft = GT_data_parser(time_p)
    fine_X, fine_Y, fine_oxy = ft[0]
    w_X = fine_X
    w_Y = fine_Y
    w_O = fine_oxy
    axs.clear()
    axs.contourf(w_X,w_Y,w_O, **kw)
    axs.set_title('Glucose, Z=16 um, time = ' +str(saving_times[i])+ ' minutes') 
    axs.invert_xaxis()
    axs.axis('scaled')
  

number_of_frames = len(saving_times)

ani = matplotlib.animation.FuncAnimation(fig,animate,blit=False, frames=number_of_frames,repeat=False)

plt.show()

ani.save('./glucose.gif', writer='imagemagick', fps=4)


# fig2, ax = plt.subplots()

# # color bar
# tp = "output00000020"
# ft = GT_data_parser(tp)
# fine_X, fine_Y, fine_glu = ft[1]
# w_X = fine_X
# w_Y = fine_Y
# w_O = fine_glu
# zmin2 = min([min(zl) for zl in w_O])
# zmax2 = max([max(zl) for zl in w_O])
# #levels = np.linspace(zmin, 0.28500001,41)
# #kw = dict(levels=levels, vmin=zmin, vmax=0.28500001, origin='lower')
# levels2 = np.linspace(0, 5.5,41)
# kw2 = dict(levels=levels2, vmin=0, vmax=5.5, origin='lower')
# cp2 = ax.contourf(w_Y,w_X,w_O, **kw2)
# cbar = plt.colorbar(cp2,format='%0.4f')
# ax.clear()


# def animate(i):
#     time_p= time_point + '%02d'%(i)
#     ft = GT_data_parser(time_p)
#     fine_X, fine_Y, fine_glu = ft[0]
#     w_X = fine_X
#     w_Y = fine_Y
#     w_G = fine_glu
#     ax.clear()
#     ax.contourf(w_X,w_Y,w_G, **kw2)
#     ax.set_title('Glutamine, Z=16 um, time = ' +str(saving_times[i])+ ' minutes') 
#     ax.invert_xaxis()
#     ax.axis('scaled')
  

# number_of_frames = len(saving_times)

# ani2 = matplotlib.animation.FuncAnimation(fig2,animate,blit=False, frames=number_of_frames,repeat=False)

# plt.show()

# ani2.save('./glutamine.gif', writer='imagemagick', fps=4)


