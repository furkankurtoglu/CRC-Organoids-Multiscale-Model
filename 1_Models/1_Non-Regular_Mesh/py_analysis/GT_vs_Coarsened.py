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

#saving_times = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 20.0])
#number_of_frames = len(saving_times)

saving_times = np.array([0.0, 0.5, 1.0, 2.0, 4.0, 8.0 , 16.0])

main_path = Path(os.getcwd()).parent

GT_list_of_files = glob.glob(str(main_path) + r"\output_Ground_Truth" + r"\*")
AM_lis_of_files = glob.glob(str(main_path) + r"\output_AM_dt_1_e2_dx_32" + r"\*")


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
time_pointAM = str(main_path) +"\output_AM_dt_1_e2_dx_32\output000000"

plot_2D = 'Y'
glc_gln_analysis = 'Y'
lac_analysis = 'N'
save_gif = 'Y'

GT_variance_analysis = 'Y'

GT_AM_convergence_test = 'N'


if plot_2D == 'Y':
    def data_parserGT (time_pointGT):
        # Fine MicroEnv Data Parsing
        fine_tuple = []
       
        if path.exists(time_pointGT + "_microenvironment0.mat"):
            fine_data = sio.loadmat( time_pointGT + "_microenvironment0.mat")['multiscale_microenvironment']
            fine_x = np.unique(fine_data[0,:])
            fine_y = np.unique(fine_data[1,:])
            fine_X, fine_Y = np.meshgrid(fine_y, fine_x)
            fine_glc = fine_data[4,np.where(fine_data[2,:] == 8)]
            fine_glc = fine_glc.reshape((len(fine_y),len(fine_x)))
            fine_glc = np.transpose(fine_glc)
            fine_gln = fine_data[5,np.where(fine_data[2,:] == 8)]
            fine_gln = fine_gln.reshape((len(fine_y),len(fine_x)))
            fine_gln = np.transpose(fine_gln)
            fine_lac = fine_data[6,np.where(fine_data[2,:] == 8)]
            fine_lac = fine_lac.reshape((len(fine_y),len(fine_x)))
            fine_lac = np.transpose(fine_lac)
            fine_glc_tuple = (fine_X, fine_Y, fine_glc)
            fine_gln_tuple = (fine_X, fine_Y, fine_gln)
            fine_lac_tuple = (fine_X, fine_Y, fine_lac)
            fine_tuple = (fine_glc_tuple, fine_gln_tuple, fine_lac_tuple)
        return fine_tuple
        
    
    def data_parserAM (time_point):
            # Fine MicroEnv Data Parsing
            fine_tuple = []
            coarse_tuple = []
            transfer_tuple = []
           
            if path.exists(time_point + "_microenvironment0.mat"):
                fine_data = sio.loadmat(time_point + "_microenvironment0.mat")['multiscale_microenvironment']
                fine_x = np.unique(fine_data[0,:])
                fine_y = np.unique(fine_data[1,:])
                fine_X, fine_Y = np.meshgrid(fine_y, fine_x)
                fine_glc = fine_data[4,np.where(fine_data[2,:] == 16)]
                fine_glc = fine_glc.reshape((len(fine_y),len(fine_x)))
                fine_glc = np.transpose(fine_glc)
                fine_gln = fine_data[5,np.where(fine_data[2,:] == 16)]
                fine_gln = fine_gln.reshape((len(fine_y),len(fine_x)))
                fine_gln = np.transpose(fine_gln)
                fine_lac = fine_data[6,np.where(fine_data[2,:] == 16)]
                fine_lac = fine_lac.reshape((len(fine_y),len(fine_x)))
                fine_lac = np.transpose(fine_lac)
                fine_glc_tuple = (fine_X, fine_Y, fine_glc)
                fine_gln_tuple = (fine_X, fine_Y, fine_gln)
                fine_lac_tuple = (fine_X, fine_Y, fine_lac)
                fine_tuple = (fine_glc_tuple, fine_gln_tuple, fine_lac_tuple)
                
            
            # Coarse MicroEnv Data Parsing
            if path.exists(time_point + "_microenvironment1.mat"):
                coarse_data = sio.loadmat(time_point + "_microenvironment1.mat")['multiscale_microenvironment']
                coarse_x = coarse_data[0,:]
                coarse_y = np.unique(fine_data[1,:])
                coarse_X, coarse_Y = np.meshgrid(coarse_y, coarse_x)
                coarse_glc = coarse_data[4,:]
                coarse_glc = np.transpose(np.tile(coarse_glc,(90,1)))
                coarse_gln = coarse_data[5,:]
                coarse_gln = np.transpose(np.tile(coarse_gln,(90,1)))
                coarse_lac = coarse_data[6,:]
                coarse_lac = np.transpose(np.tile(coarse_lac,(90,1)))
                coarse_tuple = (coarse_X, coarse_Y, coarse_glc, coarse_gln, coarse_lac)
                
                
            if path.exists(time_point + "_microenvironment2.mat"):
                transfer_region = sio.loadmat(time_point + "_microenvironment2.mat")['multiscale_microenvironment']
            return fine_tuple, coarse_tuple, transfer_tuple         
    subs_list = get_subs_name() 
    if glc_gln_analysis == 'Y':
        ## 
        figGTglc, axsGTglc = plt.subplots()
        figGTgln, axsGTgln = plt.subplots()
        figAMglc, axsAMglc = plt.subplots()
        figAMgln, axsAMgln = plt.subplots()
        
        # color bar
        tpGT = str(main_path) + "\output_Ground_Truth\output00000006"
        ftGT= data_parserGT(tpGT)
        fine_X_GT, fine_Y_GT, fine_glc_GT = ftGT[0]
        a,b,fine_gln_GT = ftGT[1]
        w_X_GT = fine_X_GT
        w_Y_GT = fine_Y_GT
        w_glc_GT = fine_glc_GT
        w_gln_GT = fine_gln_GT
        levels_GT = np.linspace(0, 17.501,41)
        levels_GT2 = np.linspace(0,5.501,41)
        kw_GT = dict(levels=levels_GT, vmin=0, vmax=17.501, origin='lower')
        kw_GT2 = dict(levels=levels_GT2, vmin=0, vmax=5.501, origin='lower')
        cp_GT = axsGTglc.contourf(w_Y_GT,w_X_GT,w_glc_GT, **kw_GT)
        cp_GT2 = axsGTgln.contourf(w_Y_GT,w_X_GT,w_gln_GT, **kw_GT2)
        cbar_GT = figGTglc.colorbar(cp_GT,format='%0.4f', ax=axsGTglc)
        cbar_GT2 = figGTgln.colorbar(cp_GT2,format='%0.4f', ax=axsGTgln)
        axsGTglc.clear()
        axsGTgln.clear()
            
        
        # color bar
        tpAM = str(main_path) + "\output_AM_dt_1_e2_dx_32\output00000006"
        ft_AM, ct_AM, tt_AM = data_parserAM(tpAM)
        fine_X_AM, fine_Y_AM, fine_glc_AM = ft_AM[0]
        fine_X_AM, fine_Y_AM, fine_gln_AM = ft_AM[1]
        cX_AM, cY_AM, c_glc_AM, c_gln_AM, c_lac_AM = ct_AM
        w_X_AM = np.concatenate((fine_X_AM,cX_AM),axis=0)
        w_Y_AM = np.concatenate((fine_Y_AM,cY_AM),axis=0)
        w_glc_AM = np.concatenate((fine_glc_AM,cY_AM),axis=0)
        w_gln_AM = np.concatenate((fine_gln_AM,cY_AM),axis=0)
        levels_AM = np.linspace(0, 17.5,41)
        levels_AM2 = np.linspace(0, 5.5,41)
        kw_AM = dict(levels=levels_AM, vmin=0, vmax=17.5, origin='lower')
        kw_AM2 = dict(levels=levels_AM2, vmin=0, vmax=5.5, origin='lower')
        cp_AM = axsAMglc.contourf(w_Y_AM,w_X_AM,w_glc_AM, **kw_AM)
        cp_AM2 = axsAMgln.contourf(w_Y_AM,w_X_AM,w_gln_AM, **kw_AM2)
        cbar_AM = figAMglc.colorbar(cp_AM,format='%0.4f',ax=axsAMglc)
        cbar_AM2 = figAMgln.colorbar(cp_AM2,format='%0.4f',ax=axsAMgln)
        axsAMglc.clear()
        axsAMgln.clear()
             
        #GT_glc = np.zeros()
        
        def animateGTglc(i):
            time_pGT= time_pointGT + '%02d'%(i)
            ftGT = data_parserGT(time_pGT)
            fine_XGT, fine_YGT, fine_glc_GT = ftGT[0]
            w_XGT = fine_XGT
            w_YGT = fine_YGT
            w_glc_GT = fine_glc_GT
            axsGTglc.clear()
            axsGTglc.contourf(w_XGT,w_YGT,w_glc_GT, **kw_GT)
            axsGTglc.set_title('Glucose Ground Truth, Z=8 um, time = ' +str(saving_times[i])+ ' minutes') 
            axsGTglc.invert_xaxis()
            axsGTglc.axis('scaled')
            return w_XGT,w_YGT,w_glc_GT
    
        def animateGTgln(i):
            time_pGT= time_pointGT + '%02d'%(i)
            ftGT = data_parserGT(time_pGT)
            fine_XGT, fine_YGT, fine_gln_GT = ftGT[1]
            w_XGT = fine_XGT
            w_YGT = fine_YGT
            w_gln_GT = fine_gln_GT
            axsGTgln.clear()
            axsGTgln.contourf(w_XGT,w_YGT,w_gln_GT, **kw_GT2)
            axsGTgln.set_title('Glutamine Ground Truth, Z=8 um, time = ' +str(saving_times[i])+ ' minutes') 
            axsGTgln.invert_xaxis()
            axsGTgln.axis('scaled')
    
        def animateAMglc(i):
            time_pAM= time_pointAM + '%02d'%(i)
            ftAM, ctAM, ttAM = data_parserAM(time_pAM)
            fine_XAM, fine_YAM, fine_glc_AM = ftAM[0]
            cX_AM, cY_AM, c_glc_AM, c_gln_AM, c_lac_AM = ctAM
            w_X_AM = np.concatenate((fine_XAM,cX_AM),axis=0)
            w_Y_AM = np.concatenate((fine_YAM,cY_AM),axis=0)
            w_glc_AM = np.concatenate((fine_glc_AM,c_glc_AM),axis=0)
            axsAMglc.clear()
            axsAMglc.contourf(w_X_AM,w_Y_AM,w_glc_AM, **kw_AM)
            axsAMglc.set_title('Glucose Adaptive Mesh, Z=16 um, time = ' +str(saving_times[i])+ ' minutes') 
            axsAMglc.invert_xaxis()
            axsAMglc.axis('scaled')
            
            
        def animateAMgln(i):
            time_pAM= time_pointAM + '%02d'%(i)
            ftAM, ctAM, ttAM = data_parserAM(time_pAM)
            fine_XAM, fine_YAM, fine_gln_AM = ftAM[1]
            cX_AM, cY_AM, c_glc_AM, c_gln_AM, c_lac_AM = ctAM
            w_X_AM = np.concatenate((fine_XAM,cX_AM),axis=0)
            w_Y_AM = np.concatenate((fine_YAM,cY_AM),axis=0)
            w_gln_AM = np.concatenate((fine_gln_AM,c_gln_AM),axis=0)
            axsAMgln.clear()
            axsAMgln.contourf(w_X_AM,w_Y_AM,w_gln_AM, **kw_AM2)
            axsAMgln.set_title('Glutamine Adaptive Mesh, Z=16 um, time = ' +str(saving_times[i])+ ' minutes') 
            axsAMgln.invert_xaxis()
            axsAMgln.axis('scaled')
           
        number_of_frames = len(saving_times)
        
        aniGTglc = matplotlib.animation.FuncAnimation(figGTglc,animateGTglc,blit=False, frames=number_of_frames,repeat=False)
        aniGTgln = matplotlib.animation.FuncAnimation(figGTgln,animateGTgln,blit=False, frames=number_of_frames,repeat=False)
    
        aniAMglc = matplotlib.animation.FuncAnimation(figAMglc,animateAMglc,blit=False, frames=number_of_frames,repeat=False)
        aniAMgln = matplotlib.animation.FuncAnimation(figAMgln,animateAMgln,blit=False, frames=number_of_frames,repeat=False)
        plt.show()
    
        if save_gif == 'Y':
            aniGTglc.save('./GT_glucose_1e4_16.gif', writer='imagemagick', fps=4)
            aniGTgln.save('./GT_glutamine_1e4_16.gif', writer='imagemagick', fps=4)
            aniAMglc.save('./AM_glucose_1e2_32.gif', writer='imagemagick', fps=4)
            aniAMgln.save('./AM_glutamine_1e2_32.gif', writer='imagemagick', fps=4)
            


if GT_variance_analysis == 'Y':
    TS_GT_data = np.zeros((7,1360800,len(saving_times)))
    for t in range(0,len(saving_times)):
        time_pGT= time_pointGT + '%02d'%(t)
        time_pAM = time_pointAM + '%02d'%(t)
        GTdata = sio.loadmat(time_pGT + "_microenvironment0.mat")['multiscale_microenvironment']
        TS_GT_data[:,:,t] = GTdata
        
    x_centers = np.unique(TS_GT_data[0,:,0])
    glc_variances = np.zeros((len(saving_times),len(np.unique(GTdata[0,:]))))
    gln_variances = np.zeros((len(saving_times),len(np.unique(GTdata[0,:]))))
    lac_variances = np.zeros((len(saving_times),len(np.unique(GTdata[0,:]))))
    for x in x_centers:
        for t in range(0,len(saving_times)):
            layer_specific_data_indices = np.where(TS_GT_data[0,:,t] == x)
            glucose_data = TS_GT_data[4,layer_specific_data_indices,t]
            glutamine_data = TS_GT_data[5,layer_specific_data_indices,t]
            lac_data = TS_GT_data[6,layer_specific_data_indices,t]
            glc_var = np.var(glucose_data)
            glc_variances[t,:] = glc_var
            gln_var = np.var(glutamine_data)
            gln_variances[t,:] = gln_var
            lac_var = np.var(lac_data)
            lac_variances[t,:] = np.var(lac_var)
        
    plt.plot(np.sum(glc_variances,axis=1))
    plt.title('Ground Truth Variances at X-dimension')
    plt.xlabel('time (min)')
    plt.ylabel('Total Variances in domain')
            
if GT_AM_convergence_test == 'Y':
    TS_GT_data = np.zeros((7,1360800,len(saving_times)))
    TS_f_AM_data = np.zeros((7,129600,len(saving_times)))
    TS_c_AM_data = np.zeros((7,152,len(saving_times)))
    for t in range(0,len(saving_times)):
        time_pGT= time_pointGT + '%02d'%(t)
        time_pAM = time_pointAM + '%02d'%(t)
        GTdata = sio.loadmat(time_pGT + "_microenvironment0.mat")['multiscale_microenvironment']
        AM_fine_data = sio.loadmat(time_pAM + "_microenvironment0.mat")['multiscale_microenvironment']
        AM_coarse_data = sio.loadmat(time_pAM + "_microenvironment1.mat")['multiscale_microenvironment']
        TS_GT_data[:,:,t] = GTdata
        TS_f_AM_data[:,:,t] = AM_fine_data
        TS_c_AM_data[:,:,t] = AM_coarse_data
    
    x_centers = np.unique(TS_GT_data[0,:,0])
    GT_glucose = np.zeros((len(saving_times),len(np.unique(GTdata[0,:]))))
    AM_glucose = np.zeros((len(saving_times),len(np.unique(GTdata[0,:]))))
    euclidean_dists = np.zeros((len(saving_times),len(np.unique(GTdata[0,:]))))
    
    for t in range(0,len(saving_times)):
        for x in range(0,len(x_centers)):
            GT_layer_specific_data_indices = np.where(TS_GT_data[0,:,t] == x_centers[x])
            AM_f_layer_specific_data_indices = np.where(TS_f_AM_data[0,:,t] == x_centers[x])
            AM_c_layer_specific_data_indices = np.where(TS_c_AM_data[0,:,t] == x_centers[x])
            GT_glucose_data = TS_GT_data[4,GT_layer_specific_data_indices,t]
            GT_glc_mean = np.mean(GT_glucose_data)
            GT_glucose[t,x] = GT_glc_mean
            
            # glutamine_data = TS_GT_data[5,GT_layer_specific_data_indices,t]
            # lac_data = TS_GT_data[6,GT_layer_specific_data_indices,t]
            # GT_glc_mean = np.mean(GT_glucose_data)
            if len(AM_f_layer_specific_data_indices[0]) != 0:
                AM_glucose_data = TS_f_AM_data[4,AM_f_layer_specific_data_indices,t]
                AM_glc_mean = np.mean(AM_glucose_data)
                AM_glucose[t,x] = AM_glc_mean
            
            if len(AM_c_layer_specific_data_indices[0]) != 0:
                AM_glc_mean = TS_c_AM_data[4,AM_c_layer_specific_data_indices,t]
                AM_glucose[t,x] = AM_glc_mean
            
            euclidean_dist = (GT_glc_mean - AM_glc_mean)**2
            euclidean_dists[t,x] = euclidean_dist

    a = euclidean_dists[5,:]
    b = euclidean_dists[10:,:]
    c = a.reshape((-1,1))
    trimmed_euc_dist = np.concatenate((c.T,b),axis=0)
    log_trim_euc = np.log(trimmed_euc_dist)
    log_trim_euc[np.where(log_trim_euc < -80)] = -80
    
    plt.figure()
    for i in range(len(log_trim_euc)):
        plt.plot(log_trim_euc[i,:])
    plt.legend(['0.5-min','1-min','2-min','3-min','4-min','5-min','6-min','7-min','8-min','9-min','10-min','20-min'],loc='lower left')
    plt.title('L2Norm Between Ground Truth and Adaptive Mesh Approach')
    plt.xlabel('Domain location in X-dimension (micrometer)')
    plt.ylabel('log-Euclidean Distance')
    plt.show()