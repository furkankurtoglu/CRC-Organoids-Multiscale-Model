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

main_path = Path(os.getcwd()).parent

GT_list_of_files = glob.glob(str(main_path) + r"\output_Ground_Truth" + r"\*")
MD_lis_of_files = glob.glob(str(main_path) + r"\output_Multiple_Domains" + r"\*")


def get_time_points(list_of_files):
    pattern = re.compile('[0-9]{8}_microenvironment0.mat')
    txt = [s for s in list_of_files if pattern.search(s)]
    print(txt.group())
    return txt
    
text = get_time_points(GT_list_of_files)    
    

    
# def get_subs_name (): 
#         tree = ET.parse("initial.xml")
#         root = tree.getroot()
        
#         subs_names = []
#         for substrate in root.iter('variable'):
#             subs_names.append(substrate.attrib['name'])
        
#         return subs_names

# subs_list = get_subs_name()
