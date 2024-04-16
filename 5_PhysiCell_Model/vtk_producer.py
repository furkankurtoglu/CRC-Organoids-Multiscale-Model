# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 06:57:28 2024

@author: Furkan
"""
import pathlib
import pyMCDS
import matplotlib.pyplot as plt
import numpy as np


time_stamp = 'output00000001'
s_path_3d = str(pathlib.Path().parent.resolve()/'output')
s_file_3d = time_stamp+'.xml'
s_pathfile_3d = f'{s_path_3d}/{s_file_3d}'

mcds = pyMCDS.pyMCDS(xmlfile=s_file_3d, output_path=s_path_3d, microenv=True, graph=False, settingxml='PhysiCell_Settings.xml')
a= mcds.get_conc_df()

#%%
a = a[a["mesh_center_m"] == 496]
y = a["mesh_center_n"].unique()
z = a["mesh_center_p"].unique()

YY,ZZ = np.meshgrid(y,z)

lac = a["lactate"]
lac = lac.values.reshape([180,180])

fig1, ax2 = plt.subplots()
CS = ax2.contourf(YY, ZZ, lac)

cbar = fig1.colorbar(CS)


print(mcds.get_concentration_at(x=-496,y=0,z=-610)[2])

#%%
#mcds.make_conc_vtk