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
import re #for regex

saving_times = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 20.0])
#number_of_frames = len(saving_times)

main_path = Path(os.getcwd()).parent

GT_list_of_files = glob.glob(str(main_path) + r"\output_Ground_Truth" + r"\*")
MD_lis_of_files = glob.glob(str(main_path) + r"\output_Multiple_Domains" + r"\*")


# def get_time_points(list_of_files):
#     pattern = re.compile('[0-9]{8}_microenvironment0.mat')
#     txt = [s for s in list_of_files if pattern.search(s)]
#     for i in txt:
        
#     return txt
#text = get_time_points(GT_list_of_files)    

    

    
def get_subs_name (): 
        tree = ET.parse(str(main_path) +"\output_Ground_Truth\initial.xml")
        root = tree.getroot()
        
        subs_names = []
        for substrate in root.iter('variable'):
            subs_names.append(substrate.attrib['name'])
        return subs_names


time_pointGT = str(main_path) +"\output_Ground_Truth\output000000"
time_pointMD = str(main_path) +"\output_Multiple_Domains\output000000"


oxy_glu_analysis = 'Y'
chem_analysis = 'Y'
def data_parserGT (time_pointGT):
    # Fine MicroEnv Data Parsing
    fine_tuple = []
    coarse_tuple = []
    transfer_tuple = []
   
    if path.exists(time_pointGT + "_microenvironment0.mat"):
        fine_data = sio.loadmat( time_pointGT + "_microenvironment0.mat")['multiscale_microenvironment']
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
    

def data_parserMD (time_point):
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
            
subs_list = get_subs_name()

if oxy_glu_analysis == 'Y':
    figGT, axsGT = plt.subplots()
    figMD, axsMD = plt.subplots()
    
    # color bar
    tpGT = str(main_path) + "\output_Ground_Truth\output00000020"
    ftGT= data_parserGT(tpGT)
    fine_X_GT, fine_Y_GT, fine_oxy_GT = ftGT[0]
    w_X_GT = fine_X_GT
    w_Y_GT = fine_Y_GT
    w_O_GT = fine_oxy_GT
    zmin_GT = min([min(zl) for zl in w_O_GT])
    zmax_GT = max([max(zl) for zl in w_O_GT])
    #levels = np.linspace(zmin, 0.28500001,41)
    #kw = dict(levels=levels, vmin=zmin, vmax=0.28500001, origin='lower')
    levels_GT = np.linspace(0, 17.5,41)
    kw_GT = dict(levels=levels_GT, vmin=0, vmax=17.5, origin='lower')
    cp_GT = axsGT.contourf(w_Y_GT,w_X_GT,w_O_GT, **kw_GT)
    cbar_GT = plt.colorbar(cp_GT,format='%0.4f')
    axsGT.clear()
    
    
        
    # color bar
    tpMD = str(main_path) + "\output_Multiple_Domains\output00000020"
    ft_MD, ct_MD, tt_MD = data_parserMD(tpMD)
    fine_X_MD, fine_Y_MD, fine_oxy_MD = ft_MD[0]
    cX_MD, cY_MD, cOxy_MD, cGlu_MD, cChem_MD = ct_MD
    w_X_MD = np.concatenate((fine_X_MD,cX_MD),axis=0)
    w_Y_MD = np.concatenate((fine_Y_MD,cY_MD),axis=0)
    w_O_MD = np.concatenate((fine_oxy_MD,cOxy_MD),axis=0)
    zmin_MD = min([min(zl) for zl in w_O_MD])
    zmax_MD = max([max(zl) for zl in w_O_MD])
    #levels = np.linspace(zmin, 0.28500001,41)
    #kw = dict(levels=levels, vmin=zmin, vmax=0.28500001, origin='lower')
    levels = np.linspace(0, 17.5,41)
    kw_MD = dict(levels=levels, vmin=0, vmax=17.5, origin='lower')
    cp_MD = axsMD.contourf(w_Y_MD,w_X_MD,w_O_MD, **kw_MD)
    cbar_MD = plt.colorbar(cp_MD,format='%0.4f')
    axsMD.clear()
        
    
    def animate(i):
        time_pGT= time_pointGT + '%02d'%(i)
        ftGT = data_parserGT(time_pGT)
        fine_XGT, fine_YGT, fine_oxyGT = ftGT[0]
        w_XGT = fine_XGT
        w_YGT = fine_YGT
        w_OGT = fine_oxyGT
        axsGT.clear()
        axsGT.contourf(w_XGT,w_YGT,w_OGT, **kw_GT)
        axsGT.set_title('Glucose Ground Truth, Z=16 um, time = ' +str(saving_times[i])+ ' minutes') 
        axsGT.invert_xaxis()
        axsGT.axis('scaled')
        time_pMD= time_pointMD + '%02d'%(i)
        ftMD, ctMD, ttMD = data_parserMD(time_pMD)
        fine_XMD, fine_YMD, fine_oxyMD = ftMD[0]
        cX_MD, cY_MD, cOxy_MD, cGlu_MD, cChem_MD = ctMD
        w_X_MD = np.concatenate((fine_XMD,cX_MD),axis=0)
        w_Y_MD = np.concatenate((fine_YMD,cY_MD),axis=0)
        w_O_MD = np.concatenate((fine_oxyMD,cOxy_MD),axis=0)
        axsMD.clear()
        axsMD.contourf(w_X_MD,w_Y_MD,w_O_MD, **kw_MD)
        axsMD.set_title('Glucose MD, Z=16 um, time = ' +str(saving_times[i])+ ' minutes') 
        axsMD.invert_xaxis()
        axsMD.axis('scaled')
        MD_tuple = [w_X_MD,w_Y_MD,w_O_MD]
        return MD_tuple
        
        
        
        
    number_of_frames = len(saving_times)
    
    ani = matplotlib.animation.FuncAnimation(figGT,animate,blit=False, frames=number_of_frames,repeat=False)
    
    plt.show()
    
    ani.save('./GT_glucose.gif', writer='imagemagick', fps=4)


#     fig2, ax = plt.subplots()
    
#     # color bar
#     tp = "\output_Ground_Truth\output00000020"
#     ft, ct, tt = data_parser(tp)
#     fine_X, fine_Y, fine_glu = ft[1]
#     cX, cY, cOxy, cGlu, cChem = ct
#     w_X = np.concatenate((fine_X,cX),axis=0)
#     w_Y = np.concatenate((fine_Y,cY),axis=0)
#     w_G = np.concatenate((fine_glu,cGlu),axis=0)
#     zmin2 = min([min(zl) for zl in w_G])
#     zmax2 = max([max(zl) for zl in w_G])
#     #levels2 = np.linspace(zmin2, zmax2)
#     #kw2 = dict(levels=levels2, vmin=zmin2, vmax=zmax2, origin='lower')
#     levels2 = np.linspace(0, 5.5)
#     kw2 = dict(levels=levels2, vmin=0, vmax=5.5, origin='lower')
#     cp2 = ax.contourf(w_X,w_Y,w_G, **kw2)
#     cbar2 = plt.colorbar(cp2,format='%0.2f')
#     ax.clear()
    
#     def animate2(i):
#         time_p= time_point + '%02d'%(i)
#         ft, ct, tt = data_parser(time_p)
#         fine_X, fine_Y, fine_glu = ft[1]
#         cX, cY, cOxy, cGlu, cChem = ct
#         w_X = np.concatenate((fine_X,cX),axis=0)
#         w_Y = np.concatenate((fine_Y,cY),axis=0)
#         w_G = np.concatenate((fine_glu,cGlu),axis=0)
#         ax.clear()
#         ax.contourf(w_X,w_Y,w_G, **kw2)
#         ax.set_title('Glutamine, Z=16 um, time = ' +str(saving_times[i])+ ' minutes') 
#         ax.invert_xaxis()
#         ax.axis('scaled')
        
    
    
#     ani2 = matplotlib.animation.FuncAnimation(fig2,animate2,blit=False, frames=number_of_frames,repeat=False)
    
#     plt.show()
    
#     ani2.save('./glutamine.gif', writer='imagemagick', fps=4)





# if chem_analysis == 'Y':
#     fig3, ax3 = plt.subplots()
    
#     # color bar
#     tp = "\output_Ground_Truth\output00000020"
#     ft, ct, tt = data_parser(tp)
#     fine_X, fine_Y, fine_chem = ft[2]
#     cX, cY, cOxy, cGlu, cChem = ct
#     w_X = np.concatenate((fine_X,cX),axis=0)
#     w_Y = np.concatenate((fine_Y,cY),axis=0)
#     w_C = np.concatenate((fine_chem,cChem),axis=0)
#     zmin3 = min([min(zl) for zl in w_C])
#     zmax3 = max([max(zl) for zl in w_C])
#     #levels3 = np.linspace(0, zmax3)
#     #kw3 = dict(levels=levels3, vmin=zmin3, vmax=zmax3, origin='lower')
#     levels3 = np.linspace(0, 16.0)
#     kw3 = dict(levels=levels3, vmin=0, vmax=16, origin='lower')
#     cp3 = ax3.contourf(w_X,w_Y,w_C, **kw3)
#     cbar3 = plt.colorbar(cp3,format='%0.5f')
#     ax3.clear()
    
#     def animate3(i):
#         time_p= time_point + '%02d'%(i)
#         ft, ct, tt = data_parser(time_p)
#         fine_X, fine_Y, fine_chem = ft[2]
#         cX, cY, cOxy, cGlu, cChem = ct
#         w_X = np.concatenate((fine_X,cX),axis=0)
#         w_Y = np.concatenate((fine_Y,cY),axis=0)
#         w_C = np.concatenate((fine_chem,cChem),axis=0)
#         ax3.clear()
#         ax3.contourf(w_X,w_Y,w_C, **kw3)
#         ax3.set_title('lactate, Z=16 um, time = ' +str(saving_times[i])+ ' minutes') 
#         ax3.invert_xaxis()
#         ax3.axis('scaled')
        
    
    
#     ani3 = matplotlib.animation.FuncAnimation(fig3,animate3,blit=False, frames=number_of_frames,repeat=False)
    
#     plt.show()
    
#     ani3.save('./lactate.gif', writer='imagemagick', fps=4)




#%%

