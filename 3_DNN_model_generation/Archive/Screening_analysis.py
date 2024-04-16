# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 14:47:59 2023

@author: Furkan
"""
from numpy import loadtxt
from keras.models import Sequential
from keras.layers import Dense
import numpy as np
import matplotlib.pyplot as plt
from keras2cpp import export_model
from sklearn.metrics import r2_score


dataname = 'Tracked_Intracellular_Metabolites_Complete_Screening.csv'

dataset = loadtxt(dataname, delimiter=',')
#%%
pep = dataset[:,13]
pep = np.array(pep)

pep_check = np.where(pep < 0)



data_no_pep1 = dataset[:,10:13]
data_no_pep1 = np.array(data_no_pep1)
no_pep_check1 = np.where(data_no_pep1 > 0)

data_no_pep2 = dataset[:,14:]
data_no_pep2 = np.array(data_no_pep2)
no_pep_check2 = np.where(data_no_pep2 > 0)