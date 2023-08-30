#
# anim_substrate_cmap.py:  render/animate PhysiCell substrate (output%08d_microenvironment0.mat files), using left/right arrows on keyboard
#
# Usage:
#  python anim_substrate_cmap.py <show_nucleus start_index axes_min axes_max>
#    i.e., the arguments <...> are optional and have defaults.
# 
# Keyboard arrows: right/left arrows will single step forward/backward; up/down will increment/decrement step size
#
# Dependencies include matplotlib and numpy. We recommend installing the Anaconda Python3 distribution.
#
# Author: Randy Heiland (except for the circles() function)
#
#
__author__ = "Randy Heiland"

import sys,pathlib
import xml.etree.ElementTree as ET
import math
import scipy.io
import matplotlib
#import matplotlib.pyplot as plt  # NB! do this AFTER the TkAgg line below!
#import matplotlib.colors as mplc
import numpy as np
from numpy.random import randn


try:
  # apparently we need mpl's Qt backend to do keypresses 
#  matplotlib.use("Qt5Agg")
  matplotlib.use("TkAgg")
  import matplotlib.pyplot as plt
except:
  print("\n---Error: cannot use matplotlib's TkAgg backend")
#  print("Consider installing Anaconda's Python 3 distribution.")
  raise


current_idx = 0

if (len(sys.argv) < 3):
  frame_idx = 0
  field_index = 4
else:
  kdx = 1
  frame_idx = int(sys.argv[kdx])
  kdx += 1
  field_index = int(sys.argv[kdx])

print('frame, field = ',frame_idx, field_index)

fig = plt.figure(figsize=(7,5.8))
#ax = fig.gca()

time_delay = 0.1
count = -1

cbar = None

#def plot_substrate(FileId):
def plot_substrate():
    global current_idx, axes_max, cbar

    # select whichever substrate index you want, e.g., for one model:
    # 4=tumor cells field, 5=blood vessel density, 6=growth substrate
#    field_index = 4  
#    field_index = 5  

    xml_file = "output%08d.xml" % current_idx
    tree = ET.parse(xml_file)
    root = tree.getroot()
#    print('time=' + root.find(".//current_time").text)
    mins = float(root.find(".//current_time").text)
    hrs = mins/60.
    days = hrs/24.
    title_str = '%d days, %d hrs, %d mins' % (int(days),(hrs%24), mins - (hrs*60))

    xml_file = "output%08d.xml" % current_idx
    tree = ET.parse(xml_file)
    root = tree.getroot()
    print('time=' + root.find(".//current_time").text)

#     print('debug> plot_substrate: idx=',FileId)
#    fname = "output%08d_microenvironment0.mat" % FileId
    fname = "output%08d_microenvironment0.mat" % current_idx
    output_dir_str = '.'
    fullname = output_dir_str + "/" + fname
    if not pathlib.Path(fullname).is_file():
        print("file not found",fullname)
        return

    info_dict = {}
    scipy.io.loadmat(fullname, info_dict)
    M = info_dict['multiscale_microenvironment']
#     global_field_index = int(mcds_field.value)
    print('plot_substrate: field_index=',field_index)
    f = M[field_index,:]   # 
    #plt.clf()
    #my_plot = plt.imshow(f.reshape(400,400), cmap='jet', extent=[0,20, 0,20])
    
#    fig = plt.figure(figsize=(7,7))
#    fig = plt.figure(figsize=(7,5.8))

    N = int(math.sqrt(len(M[0,:])))
    grid2D = M[0,:].reshape(N,N)
    xvec = grid2D[0,:]
    #xvec.size
    #xvec.shape
    num_contours = 30
    num_contours = 10
#    my_plot = plt.contourf(xvec,xvec,M[field_index,:].reshape(100,100), num_contours, cmap='viridis') #'viridis'
    my_plot = plt.contourf(xvec,xvec,M[field_index,:].reshape(N,N), num_contours, cmap='viridis') #'viridis'
#    cbar.remove()
    if cbar == None:
      cbar = plt.colorbar(my_plot)
    else:
      cbar = plt.colorbar(my_plot, cax=cbar.ax)
    axes_min = 0
    axes_min = -2000
    axes_max = 2000
#    plt.xlim(axes_min,axes_max)
#    plt.ylim(axes_min,axes_max)
    # plt.title(fname)
    plt.title(title_str)
#     ax.set_title(fname)
#    plt.axis('equal')

#    plt.show()
    plt.pause(time_delay)



step_value = 1
def press(event):
  global current_idx, step_value
#    print('press', event.key)
  sys.stdout.flush()
  if event.key == 'escape':
    sys.exit(1)
  elif event.key == 'h':  # help
    print('esc: quit')
    print('right arrow: increment by step_value')
    print('left arrow:  decrement by step_value')
    print('up arrow:   increment step_value by 1')
    print('down arrow: decrement step_value by 1')
    print('0: reset to 0th frame')
    print('h: help')
  elif event.key == 'left':  # left arrow key
#    print('go backwards')
#    fig.canvas.draw()
    current_idx -= step_value
    if (current_idx < 0):
      current_idx = 0
    plot_substrate()
  elif event.key == 'right':  # right arrow key
#        print('go forwards')
#        fig.canvas.draw()
    current_idx += step_value
    plot_substrate()
  elif event.key == 'up':  # up arrow key
    step_value += 1
    print('step_value=',step_value)
  elif event.key == 'down':  # down arrow key
    step_value -= 1
    if (step_value <= 0):
      step_value = 1
    print('step_value=',step_value)
  elif event.key == '0':  # reset to 0th frame/file
    current_idx = 0
    plot_substrate()
  else:
    print('press', event.key)

plot_substrate()
print("\nNOTE: click in plot window to give it focus before using keys.")

fig.canvas.mpl_connect('key_press_event', press)
#plot_substrate(frame_idx)
plt.show()

