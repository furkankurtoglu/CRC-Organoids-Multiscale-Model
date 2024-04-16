import numpy as np

import pathlib
import pcdl
import pandas as pd


output_list = range(10,32)



time_stamp = 'output00000001'
s_path_3d = str(pathlib.Path().parent.resolve()/'output')
s_file_3d = time_stamp+'.xml'
s_pathfile_3d = f'{s_path_3d}/{s_file_3d}'

# Reading pcdl
mcds = pcdl.pyMCDS(xmlfile=s_file_3d, output_path=s_path_3d, microenv=True, graph=False, settingxml='PhysiCell_Settings.xml')
mcds.make_conc_vtk()