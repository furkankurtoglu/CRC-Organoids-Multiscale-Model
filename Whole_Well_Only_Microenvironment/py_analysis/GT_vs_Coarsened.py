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
MD_out_path = os.path.join(main_path, "output_Multiple_Domains")

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
    
ft = GT_data_parser()