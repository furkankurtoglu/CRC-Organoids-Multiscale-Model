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
out_path = os.path.join(main_path, "output")

os.chdir(out_path)

time_point = "output000000"
number_of_frames = len(saving_times)

Temporospatial_Plotting = 'N'
Total_Amount_Analysis = 'Y'





if Temporospatial_Plotting == 'Y':
    oxy_glu_analysis = 'Y'
    chem_analysis = 'Y'
    def data_parser (time_point):
        # Fine MicroEnv Data Parsing
        fine_tuple = []
        coarse_tuple = []
        transfer_tuple = []
       
        if path.exists(time_point + "_microenvironment0.mat"):
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
            
        
        # Coarse MicroEnv Data Parsing
        if path.exists(time_point + "_microenvironment1.mat"):
            coarse_data = sio.loadmat(time_point + "_microenvironment1.mat")['multiscale_microenvironment']
            coarse_x = coarse_data[0,:]
            coarse_y = np.unique(fine_data[1,:])
            coarse_X, coarse_Y = np.meshgrid(coarse_y, coarse_x)
            coarse_oxy = coarse_data[4,:]
            coarse_oxy = np.transpose(np.tile(coarse_oxy,(90,1)))
            coarse_glu = coarse_data[5,:]
            coarse_glu = np.transpose(np.tile(coarse_glu,(90,1)))
            coarse_chem = coarse_data[6,:]
            coarse_chem = np.transpose(np.tile(coarse_chem,(90,1)))
            coarse_tuple = (coarse_X, coarse_Y, coarse_oxy, coarse_glu, coarse_chem)
            
            
        if path.exists(time_point + "_microenvironment2.mat"):
            transfer_region = sio.loadmat(time_point + "_microenvironment2.mat")['multiscale_microenvironment']
    
        return fine_tuple, coarse_tuple, transfer_tuple
        
    def get_subs_name ():
        tree = ET.parse("initial.xml")
        root = tree.getroot()
        
        subs_names = []
        for substrate in root.iter('variable'):
            subs_names.append(substrate.attrib['name'])
        
        return subs_names
            
    subs_list = get_subs_name()
    
    if oxy_glu_analysis == 'Y':
        fig, axs = plt.subplots()
        
        # color bar
        tp = "output00000020"
        ft, ct, tt = data_parser(tp)
        fine_X, fine_Y, fine_oxy = ft[0]
        cX, cY, cOxy, cGlu, cChem = ct
        w_X = np.concatenate((fine_X,cX),axis=0)
        w_Y = np.concatenate((fine_Y,cY),axis=0)
        w_O = np.concatenate((fine_oxy,cOxy),axis=0)
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
            ft, ct, tt = data_parser(time_p)
            fine_X, fine_Y, fine_oxy = ft[0]
            cX, cY, cOxy, cGlu, cChem = ct
            w_X = np.concatenate((fine_X,cX),axis=0)
            w_Y = np.concatenate((fine_Y,cY),axis=0)
            w_O = np.concatenate((fine_oxy,cOxy),axis=0)
            axs.clear()
            axs.contourf(w_X,w_Y,w_O, **kw)
            axs.set_title('Glucose, Z=16 um, time = ' +str(saving_times[i])+ ' minutes') 
            axs.invert_xaxis()
            axs.axis('scaled')
            
            
            
            
        number_of_frames = len(saving_times)
        
        ani = matplotlib.animation.FuncAnimation(fig,animate,blit=False, frames=number_of_frames,repeat=False)
        
        plt.show()
        
        ani.save('./glucose.gif', writer='imagemagick', fps=4)
    
    
        fig2, ax = plt.subplots()
        
        # color bar
        tp = "output00000020"
        ft, ct, tt = data_parser(tp)
        fine_X, fine_Y, fine_glu = ft[1]
        cX, cY, cOxy, cGlu, cChem = ct
        w_X = np.concatenate((fine_X,cX),axis=0)
        w_Y = np.concatenate((fine_Y,cY),axis=0)
        w_G = np.concatenate((fine_glu,cGlu),axis=0)
        zmin2 = min([min(zl) for zl in w_G])
        zmax2 = max([max(zl) for zl in w_G])
        #levels2 = np.linspace(zmin2, zmax2)
        #kw2 = dict(levels=levels2, vmin=zmin2, vmax=zmax2, origin='lower')
        levels2 = np.linspace(0, 5.5)
        kw2 = dict(levels=levels2, vmin=0, vmax=5.5, origin='lower')
        cp2 = ax.contourf(w_X,w_Y,w_G, **kw2)
        cbar2 = plt.colorbar(cp2,format='%0.2f')
        ax.clear()
        
        def animate2(i):
            time_p= time_point + '%02d'%(i)
            ft, ct, tt = data_parser(time_p)
            fine_X, fine_Y, fine_glu = ft[1]
            cX, cY, cOxy, cGlu, cChem = ct
            w_X = np.concatenate((fine_X,cX),axis=0)
            w_Y = np.concatenate((fine_Y,cY),axis=0)
            w_G = np.concatenate((fine_glu,cGlu),axis=0)
            ax.clear()
            ax.contourf(w_X,w_Y,w_G, **kw2)
            ax.set_title('Glutamine, Z=16 um, time = ' +str(saving_times[i])+ ' minutes') 
            ax.invert_xaxis()
            ax.axis('scaled')
            
        
        
        ani2 = matplotlib.animation.FuncAnimation(fig2,animate2,blit=False, frames=number_of_frames,repeat=False)
        
        plt.show()
        
        ani2.save('./glutamine.gif', writer='imagemagick', fps=4)
    
    
    
    
    
    if chem_analysis == 'Y':
        fig3, ax3 = plt.subplots()
        
        # color bar
        tp = "output00000020"
        ft, ct, tt = data_parser(tp)
        fine_X, fine_Y, fine_chem = ft[2]
        cX, cY, cOxy, cGlu, cChem = ct
        w_X = np.concatenate((fine_X,cX),axis=0)
        w_Y = np.concatenate((fine_Y,cY),axis=0)
        w_C = np.concatenate((fine_chem,cChem),axis=0)
        zmin3 = min([min(zl) for zl in w_C])
        zmax3 = max([max(zl) for zl in w_C])
        #levels3 = np.linspace(0, zmax3)
        #kw3 = dict(levels=levels3, vmin=zmin3, vmax=zmax3, origin='lower')
        levels3 = np.linspace(0, 16.0)
        kw3 = dict(levels=levels3, vmin=0, vmax=16, origin='lower')
        cp3 = ax3.contourf(w_X,w_Y,w_C, **kw3)
        cbar3 = plt.colorbar(cp3,format='%0.5f')
        ax3.clear()
        
        def animate3(i):
            time_p= time_point + '%02d'%(i)
            ft, ct, tt = data_parser(time_p)
            fine_X, fine_Y, fine_chem = ft[2]
            cX, cY, cOxy, cGlu, cChem = ct
            w_X = np.concatenate((fine_X,cX),axis=0)
            w_Y = np.concatenate((fine_Y,cY),axis=0)
            w_C = np.concatenate((fine_chem,cChem),axis=0)
            ax3.clear()
            ax3.contourf(w_X,w_Y,w_C, **kw3)
            ax3.set_title('lactate, Z=16 um, time = ' +str(saving_times[i])+ ' minutes') 
            ax3.invert_xaxis()
            ax3.axis('scaled')
            
        
        
        ani3 = matplotlib.animation.FuncAnimation(fig3,animate3,blit=False, frames=number_of_frames,repeat=False)
        
        plt.show()
        
        ani3.save('./lactate.gif', writer='imagemagick', fps=4)
    



#%%



if Total_Amount_Analysis == 'Y':
    o2_uptake_rate_per_cell = 0.005
    glu_uptake_rate_per_cell = 0.01
    chem_secretion_rate_per_cell_per_min = 0.01
    number_of_cells = 170278
    
    
    total_O2 = []
    total_glu = []
    total_chem = []
    
    initial_O2= 0;
    
    previous_data = np.array([0,0,0])
    previous_time = 0;
    
    for i in range(number_of_frames):
        time_p = time_point + '%02d'%(i)
        if path.exists(time_p + "_microenvironment0.mat"):
            fine_data = sio.loadmat(time_p + "_microenvironment0.mat")['multiscale_microenvironment']
            voxel_volume = 32768/(10**15) # liters
            cell_volume = 2494/(10**15) # liters
            fine_volume = voxel_volume * len(fine_data[0])
            coarse_data = sio.loadmat(time_p + "_microenvironment1.mat")['multiscale_microenvironment']
            c_voxel_volume = voxel_volume*8100
            coarse_volume = c_voxel_volume*len(coarse_data[0])
            micEnv_O2 =sum(np.multiply(np.asarray(fine_data[4,:]), voxel_volume)) # mmol
            micEnv_glu = sum(np.multiply(np.asarray(fine_data[5,:]), voxel_volume)) # mmol
            micEnv_chem = sum(np.multiply(np.asarray(fine_data[6,:]), voxel_volume)) # mmol
            cEnv_O2 =sum(np.multiply(np.asarray(coarse_data[4,:]), c_voxel_volume)) # mmol
            cEnv_glu = sum(np.multiply(np.asarray(coarse_data[5,:]), c_voxel_volume)) # mmol
            cEnv_chem = sum(np.multiply(np.asarray(coarse_data[6,:]), c_voxel_volume)) # mmol
            '''
            if i == 0:
                prev_oxy = micEnv_O2
                initial_O2 = micEnv_O2
                initial_glu = micEnv_glu
                initial_chem = micEnv_chem
                prev_oxy_density = (micEnv_O2.copy())/fine_volume # mmol/L ; mM
                O2_diff.append(micEnv_O2)
            if i > 0:
                #uptaken_glu = glu_uptake_rate_per_cell*dx*dx*dx * number_of_cells * saving_times[i]
                uptaken_O2 = cell_volume * number_of_cells * o2_uptake_rate_per_cell * prev_oxy_density * time_diff[i-1] # mmol
                #uptaken_O2 = o2_uptake_rate_per_cell* prev_oxy * number_of_cells * saving_times[i]
                #uptaken_O2 = o2_uptake_rate_per_cell * number_of_cells * saving_times[i]
                #uptaken_O2 = o2_uptake_rate_per_cell* prev_oxy * saving_times[i]
                O2_diff.append(prev_oxy_density*fine_volume - uptaken_O2)
                prev_oxy_density = (micEnv_O2.copy())/fine_volume
                
            '''    
    
            total_O2.append(micEnv_O2 + cEnv_O2)
            total_glu.append(micEnv_glu + cEnv_glu)
            total_chem.append(micEnv_chem + cEnv_chem)

            
    print("Total Oxygen Amount (mmol)")
    print(total_O2)        
    print("------------------")
    print("Total Glucose Amount (mmol)")
    print(total_glu)
    print("------------------")
    print("Total Chemokine Amount (mmol)")
    print(total_chem)
    plt.figure()
    plt.plot(saving_times, total_O2)
    plt.title('Total glucose is conserved during simulation (mmol)')
    plt.xlabel('time(min)')
    plt.ylabel('Amount (mmol)')
    plt.ylim((0., 0.0014))
    plt.show()
    plt.figure()
    plt.plot(saving_times, total_glu)
    plt.title('Total Glutamine (mmol)')
    plt.xlabel('time(min)')
    plt.ylabel('Amount (mmol)')
    plt.ylim((0, .0005))
    plt.show()
    plt.figure()
    plt.plot(saving_times, total_chem)
    plt.title('Total Lactate (mmol)')
    plt.xlabel('time(min)')
    plt.ylabel('Amount (mmol)')
    plt.ylim((0, .0001))
    plt.show() 
    
    
    
