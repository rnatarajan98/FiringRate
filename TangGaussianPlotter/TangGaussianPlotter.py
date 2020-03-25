# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from scipy.io import loadmat
from cell import cell


x = loadmat('C:\\Uni\\TangGaussianPlotter\\cells.mat')

cells = []

i=0

#1. Iterate through the list of cells in x
for tempCellMat in range(x["cells"].size):
    x_tolist = x["cells"][0,tempCellMat].tolist()
    x_tolist[0]

    #2. Form a temporary cell from the fields
    #TODO: simplify this?
    tempCellPy = cell()
    tempCellPy.name = x_tolist[0]
    tempCellPy.Spont_mean = x_tolist[1]
    tempCellPy.Spont_SEM = x_tolist[2]
    tempCellPy.rf_sta_lat = x_tolist[3]
    tempCellPy.rf_xdeg = x_tolist[4]
    tempCellPy.rf_ydeg = x_tolist[5]
    tempCellPy.rf_center = x_tolist[6]
    tempCellPy.rf_sigma = x_tolist[7]
    tempCellPy.sf_tuning = x_tolist[8]
    tempCellPy.sf_preferred = x_tolist[9]
    tempCellPy.sf_cutoff =  x_tolist[10]
    tempCellPy.sf_latency = x_tolist[11]
    tempCellPy.sf_spikes = x_tolist[12]
    tempCellPy.tf_tuning = x_tolist[13]
    tempCellPy.tf_spikes = x_tolist[14]
    tempCellPy.c_tuning = x_tolist[15]
    tempCellPy.c_gain = x_tolist[16]
    tempCellPy.c_c50 = x_tolist[17]
    tempCellPy.c_fit = x_tolist[18]
    tempCellPy.c_spikes = x_tolist[19]
    tempCellPy.d_tuning = x_tolist[20]
    tempCellPy.d_preferred = x_tolist[21]
    tempCellPy.d_DSI = x_tolist[22]
    tempCellPy.d_OSI = x_tolist[23]
    tempCellPy.d_spikes = x_tolist[24]
    tempCellPy.f_transient_index = x_tolist[25]
    tempCellPy.f_spikes = x_tolist[26]
    tempCellPy.sf_curvefit = x_tolist[27]
    tempCellPy.c_curvefit = x_tolist[28]

    #Add the temporary cell to the full list of cells
    cells.append(tempCellPy)

rf_sigma_temp = cells[0].rf_sigma[0]
rf_sigma_avg = (rf_sigma_temp[0]+rf_sigma_temp[1])/2
print(rf_sigma_avg)