# plot_cell_data.py:  plot a cell's custom data value over time
#
# Usage:
#  python plot_cell_data.py <...>
#    i.e., the arguments <...> are optional and have defaults.
# 
# Dependencies include matplotlib and numpy. We recommend installing the Anaconda Python3 distribution.
#
# Examples (run from directory containing the .xml and .mat files):
#  python 
#
# Author: Randy Heiland 
#
__author__ = "Randy Heiland"

import sys
import glob
import os
import xml.etree.ElementTree as ET
import math
from pyMCDS import pyMCDS
#join_our_list = "(Join/ask questions at https://groups.google.com/forum/#!forum/physicell-users)\n"
join_our_list = "( Submit a ticket at https://sourceforge.net/p/physicell/tickets/)\n"
try:
  import matplotlib
  import matplotlib.colors as mplc
  from matplotlib.patches import Circle, Ellipse, Rectangle
  from matplotlib.collections import PatchCollection
except:
  print("\n---Error: cannot import matplotlib")
  print("---Try: python -m pip install matplotlib")
  print(join_our_list)
#  print("---Consider installing Anaconda's Python 3 distribution.\n")
  raise
try:
  import numpy as np  # if mpl was installed, numpy should have been too.
except:
  print("\n---Error: cannot import numpy")
  print("---Try: python -m pip install numpy\n")
  print(join_our_list)
  raise
from collections import deque
try:
  # apparently we need mpl's Qt backend to do keypresses 
#  matplotlib.use("Qt5Agg")
  matplotlib.use("TkAgg")
  import matplotlib.pyplot as plt
except:
  print("\n---Error: cannot use matplotlib's TkAgg backend")
  print(join_our_list)
#  print("Consider installing Anaconda's Python 3 distribution.")
  raise


print("# args=",len(sys.argv)-1)

xml_files = glob.glob('output*.xml')
xml_files.sort()

ds_count = len(xml_files)
mcds = [pyMCDS(xml_files[i], '.') for i in range(ds_count)]

def cell_data_plot(xname, yname_list, t):
  tname = "time"
  discrete_cells_names = ['virion', 'assembled_virion']
  tval = np.linspace(0, mcds[-1].get_time(), len(xml_files))

  if xname == tname:
     xval = tval
  elif xname in discrete_cells_names:
     xval = np.array([mcds[i].data['discrete_cells'][xname].sum() for i in range(ds_count)])
  else:
    if xname == 'susceptible_cells':
      xval = np.array([(mcds[i].data['discrete_cells']['assembled_virion'] <= 1).sum() for i in range(ds_count)])
      + np.array([(mcds[i].data['discrete_cells']['cycle_model'] < 6).sum() for i in range(ds_count)])
    elif xname == 'infected_cells':
      xval = np.array([(mcds[i].data['discrete_cells']['assembled_virion'] > 1).sum() for i in range(ds_count)]) \
      + np.array([(mcds[i].data['discrete_cells']['cycle_model'] < 6).sum() for i in range(ds_count)])
    elif xname == 'dead_cells':
      xval = np.array([len(mcds[0].data['discrete_cells']['ID']) - len(mcds[i].data['discrete_cells']['ID']) for i in range(ds_count)]) \
      + np.array([(mcds[i].data['discrete_cells']['cycle_model'] >= 6).sum() for i in range(ds_count)])

  for yname in yname_list:
    if yname in discrete_cells_names:
      yval = np.array([mcds[i].data['discrete_cells'][yname].sum() for i in range(ds_count)])
    else:
      if yname == 'susceptible_cells':
        yval = np.array([(mcds[i].data['discrete_cells']['assembled_virion'] <= 1).sum() for i in range(ds_count)])
        + np.array([(mcds[i].data['discrete_cells']['cycle_model'] < 6).sum() for i in range(ds_count)])
      elif yname == 'infected_cells':
        yval = np.array([(mcds[i].data['discrete_cells']['assembled_virion'] > 1).sum() for i in range(ds_count)])
        + np.array([(mcds[i].data['discrete_cells']['cycle_model'] < 6).sum() for i in range(ds_count)])
      elif yname == 'dead_cells':
        yval = np.array([len(mcds[0].data['discrete_cells']['ID']) - len(mcds[i].data['discrete_cells']['ID']) for i in range(ds_count)]) \
        + np.array([(mcds[i].data['discrete_cells']['cycle_model'] >= 6).sum() for i in range(ds_count)])
    p = plt.plot(xval, yval, label=yname)
    if (t >= 0):
      plt.plot(xval[t], yval[t], p[-1].get_color(), marker='o')
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.margins(0)

  plt.xlabel('total ' * (xname != tname) + xname)
  plt.ylabel('total ' + (yname_list[0] if len(yname_list) == 1 else ', '.join(yname_list)))
  plt.legend()
  plt.tight_layout()
  plt.show()

# cell_data_plot('susceptible_cells', ['infected_cells'], 20)
#cell_data_plot('susceptible_cells', ['infected_cells'], -1)

# time series
cell_data_plot('time', ['assembled_virion'], -1)
# multiple time series
cell_data_plot('time', ['susceptible_cells', 'infected_cells', 'dead_cells'], -1)
# multiple time series with current time
cell_data_plot('time', ['susceptible_cells', 'infected_cells', 'dead_cells'], 20)
# phase diagram
cell_data_plot('susceptible_cells', ['infected_cells'], -1)
