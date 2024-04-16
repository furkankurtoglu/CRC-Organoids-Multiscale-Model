# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 09:41:37 2024

@author: Furkan
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import csv


cell_type = 'WT'
train_model = 'Y'
save_model = 'Y'
draw_convergence = 'N'
Regression_Testing = 'N'
save_test_CSV = 'N'






# load the dataset
dataname = cell_type + '_in_silico_data_GLC_GLN_and_Seven_Metabolites_5_even.csv'

dataset = []
with open(dataname,encoding="utf-8-sig") as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    for row in csv_reader:
        # print(row[:6])
        if (row[0] != 'Glucose_exchange_input'):
            dataset.append(row[:10])
        
np_dataset = np.asarray(dataset,dtype=np.float32)
#%%

biomass = np_dataset[:,-1]

plt.figure()
plt.hist(biomass, bins=30, color='skyblue', edgecolor='black')





column_values = ['Glucose_exchange', 'Glutamine_exchange', 'G6P_intra', 'FBP_intra','G3P_intra','PEP_intra','LAC_intra','GLN_intra','GLU_intra','Biomass'] 

df_dataset = pd.DataFrame(data = np_dataset, columns = column_values)


biomass_between_0_031_and_0_037 = df_dataset[ (df_dataset['Biomass'] < 0.038) & (df_dataset['Glucose_exchange'] > 0.110) & (df_dataset['Glutamine_exchange'] > 0.0013)]
#biomass_between_0_031_and_0_037 = df_dataset[(df_dataset['Biomass'] > 0.031) & (df_dataset['Biomass'] < 0.038)]
plt.figure()
plt.hist(biomass_between_0_031_and_0_037['Biomass'], bins=30, color='skyblue', edgecolor='black')



plt.figure()
G6P_n, G6P_bins, G6P_patches = plt.hist(biomass_between_0_031_and_0_037['G6P_intra'], bins=3, color='skyblue', edgecolor='black')
plt.title('G6P_intra')
G6P_n = G6P_n/ 33787

plt.figure()
FBP_n, FBP_bins, FBP_patches = plt.hist(biomass_between_0_031_and_0_037['FBP_intra'], bins=5, color='skyblue', edgecolor='black')


plt.title('FBP_intra')


plt.figure()
G3P_n, G3P_bins, G3P_patches = plt.hist(biomass_between_0_031_and_0_037['G3P_intra'], bins=5, color='skyblue', edgecolor='black')
plt.title('G3P_intra')
G3P_n = G3P_n/ 33787


plt.figure()
PEP_n, PEP_bins, PEP_patches = plt.hist(biomass_between_0_031_and_0_037['PEP_intra'], bins=5, color='skyblue', edgecolor='black')
plt.title('PEP_intra')
PEP_n = PEP_n / 33787

plt.figure()
LAC_n, LAC_bins, LAC_patches = plt.hist(biomass_between_0_031_and_0_037['LAC_intra'], bins=3, color='skyblue', edgecolor='black')
plt.title('LAC_intra')

plt.figure()
GLN_n, GLN_bins, GLN_patches = plt.hist(biomass_between_0_031_and_0_037['GLN_intra'], bins=5, color='skyblue', edgecolor='black')
plt.title('GLN_intra')
GLN_n = GLN_n / 33787


plt.figure()
GLU_n, GLU_bins, GLU_patches =plt.hist(biomass_between_0_031_and_0_037['GLU_intra'], bins=5, color='skyblue', edgecolor='black')
plt.title('GLU_intra')
GLU_n = GLU_n / 33787
