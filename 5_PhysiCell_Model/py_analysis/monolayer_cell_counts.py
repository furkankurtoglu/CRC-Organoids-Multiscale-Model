# -*- coding: utf-8 -*-
"""
Created on Mon May 22 15:19:28 2023

@author: Furkan
"""


import numpy as np
import matplotlib.pyplot as plt



M4_time_points = np.array([0, 12,24, 36, 48, 60, 72])
M4_reps = np.array([[1200,1215,2400,4612,5143,9586,16510],[1200,1215,2400,4544,4987,9538,14488],[1200,1215,2400,4542,4996,9555,14395]])

M4_means = np.array([np.mean(M4_reps[:,0]),np.mean(M4_reps[:,1]),np.mean(M4_reps[:,2]),np.mean(M4_reps[:,3]),np.mean(M4_reps[:,4]),np.mean(M4_reps[:,5]),np.mean(M4_reps[:,6])])
M4_std = np.array([np.std(M4_reps[:,0]),np.std(M4_reps[:,1]),np.std(M4_reps[:,2]),np.std(M4_reps[:,3]),np.std(M4_reps[:,4]),np.std(M4_reps[:,5]),np.std(M4_reps[:,6])])


Exp_data = np.array([[1205,1225,1281,1083,1091],[2699,2729,2927,3027,3061],[12408,16888,16835,12730,18499]])
Exp_time_points = np.array([[0,0,0,0,0],[24,24,24,24,24],[72,72,72,72,72]])







plt.figure()
plt.plot(Exp_time_points[0,:],Exp_data[0,:],'b^')
plt.plot(Exp_time_points[1,:],Exp_data[1,:],'b^')
plt.plot(Exp_time_points[2,:],Exp_data[2,:],'b^', label = 'Experimental Cell Counts')


plt.plot(M4_time_points,M4_means,'r',label = 'Computational Simulation Cell Counts')
plt.fill_between(M4_time_points, (M4_means - M4_std),(M4_means + M4_std), color = 'r', alpha=.2)
plt.ylabel('Number of Cells')
plt.xlabel('Time (hour)')
plt.legend(loc = 'upper left')