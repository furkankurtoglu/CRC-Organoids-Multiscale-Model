import sys
import os
import glob
import numpy as np
#from pyMCDS_cells import pyMCDS_cells
import scipy.io
import matplotlib.pyplot as plt

argc = len(sys.argv)-1
print("# args=",argc)

if (argc < 1):
#  data_dir = int(sys.argv[kdx])
  print("Usage: provide output subdir")
  sys.exit(-1)

#data_dir = 'output'
kdx = 1
data_dir = sys.argv[kdx]
print('data_dir = ',data_dir)
os.chdir(data_dir)
#xml_files = glob.glob('output/output*.xml')
#xml_files = glob.glob('output*.xml')
mat_files = glob.glob('output*ment0.mat')                                                                             
os.chdir('..')
mat_files.sort()
#print('xml_files = ',xml_files)

file_count = len(mat_files)
print("----- file_count = ",file_count)
#mcds = [pyMCDS_cells(xml_files[i], data_dir) for i in range(file_count)]

#tval = np.linspace(0, mcds[-1].get_time(), file_count)
#print('tval= ',tval)

info_dict = {}
    scipy.io.loadmat(fullname, info_dict)
    M = info_dict['multiscale_microenvironment']
    # type(M)
    f=M[8,:]

collagen_val = np.array( [np.floor(mcds[idx].data['discrete_cells']['assembled_virion']).sum()  for idx in range(ds_count)] ).astype(int)
print(y_load)

plt.plot(tval,collagen_val)
plt.title(data_dir + ": collagen")
#plt.savefig(data_dir + '.png')
plt.show()
