import sys
import os
import glob
import numpy as np
from pyMCDS_cells import pyMCDS_cells
import matplotlib.pyplot as plt

argc = len(sys.argv)-1
print("# args=",argc)

#data_dir = 'output'
if (argc < 1):
#  data_dir = int(sys.argv[kdx])
  print("Usage: provide output subdir")
  sys.exit(-1)

kdx = 1
data_dir = sys.argv[kdx]
print('data_dir = ',data_dir)
os.chdir(data_dir)
xml_files = glob.glob('output*.xml')
os.chdir('..')
xml_files.sort()
#print('xml_files = ',xml_files)

ds_count = len(xml_files)
print("----- ds_count = ",ds_count)
mcds = [pyMCDS_cells(xml_files[i], data_dir) for i in range(ds_count)]

tval = np.linspace(0, mcds[-1].get_time(), ds_count)
print('tval= ',tval)

# # mac,neut,cd8,DC,cd4,Fib
# Macs
#yval4 = np.array( [(np.count_nonzero((mcds[idx].data['discrete_cells']['cell_type'] == 4) & (mcds[idx].data['discrete_cells']['cycle_model'] < 100.) == True)) for idx in range(ds_count)] )

# count Neuts
#yval5 = np.array( [(np.count_nonzero((mcds[idx].data['discrete_cells']['cell_type'] == 5) & (mcds[idx].data['discrete_cells']['cycle_model'] < 100.) == True)) for idx in range(ds_count)] )

# count CD8
yval6 = np.array( [(np.count_nonzero((mcds[idx].data['discrete_cells']['cell_type'] == 3) & (mcds[idx].data['discrete_cells']['cycle_model'] < 100.) == True)) for idx in range(ds_count)] )

# count DC
#yval7 = np.array( [(np.count_nonzero((mcds[idx].data['discrete_cells']['cell_type'] == 6) & (mcds[idx].data['discrete_cells']['cycle_model'] < 100.) == True)) for idx in range(ds_count)] )

# count CD4
#yval8 = np.array( [(np.count_nonzero((mcds[idx].data['discrete_cells']['cell_type'] == 7) & (mcds[idx].data['discrete_cells']['cycle_model'] < 100.) == True)) for idx in range(ds_count)] )

# count Fibroblasts
yval9 = np.array( [(np.count_nonzero((mcds[idx].data['discrete_cells']['cell_type'] == 8) & (mcds[idx].data['discrete_cells']['cycle_model'] < 100.) == True)) for idx in range(ds_count)] )

#plt.plot(tval, yval4, label='Mac', linewidth=1, color='lime')
#plt.plot(tval, yval5, linestyle='dashed', label='Neut', linewidth=1, color='cyan')
plt.plot(tval, yval6, label='CD8', linewidth=1, color='red')
#plt.plot(tval, yval7, linestyle='dashed', label='DC', linewidth=1, color='fuchsia')
#plt.plot(tval, yval8, label='CD4', linewidth=1, color='orange')
plt.plot(tval, yval9, linestyle='dashed',  label='Fib', linewidth=1, color='orange')

plt.xticks([1440,2*1440, 3*1440, 4*1440, 5*1440, 6*1440, 7*1440, 8*1440, 9*1440], ('1', '2','3','4','5','6','7','8','9'))

plt.legend(loc='center left', prop={'size': 10})

#plt.xlabel('Time (mins)')
plt.xlabel('Time (days)')

plt.ylabel('Number of cells')

#plt.xlim(0,11520)  # 8 days
plt.xlim(0,12960)  # 9 days

plt.title(data_dir)
plt.savefig(data_dir + '.png')
plt.show()
