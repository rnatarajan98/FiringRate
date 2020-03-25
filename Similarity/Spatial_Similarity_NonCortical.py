# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 12:36:27 2020
@author: Richard
"""
import csv
import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
import pylgn
import quantities as pq
import numpy as np
import pylgn.kernels as kernel
import my_kernels as mk
from scipy.io import loadmat
from scipy.stats import pearsonr
class Reading:
    def __init__(self, mouse, id, response):
        self.mouse = mouse
        self.id = id
        self.response = response
        self.valid = 0
        self.exp_response_X = []
        self.exp_response_Y = []
@dataclass
class mu_Reading:
    mouse: str = None
    id: str = None
    dsi: bool = False
    osi: float = None
    dir_fr: float = None
    sf_X: float = None
    sf_Y: float = None
    sf: float = None
    sf_cutoff: float = None
    tf_X: float = None
    tf_Y: float = None
    tf: float = None
    tf_bandwidth: float = None
def similarity(cell, plot):
    array1 = []
    array2 = []
    stretch = 18
    for i in range(50):
        array1.append(np.squeeze(cell.exp_response_X)[i*stretch])
        array2.append(np.squeeze(cell.exp_response_Y)[i*stretch])
    corr, _ = pearsonr(cell.response, array2)
    if plot == 1:
        plt.figure(figsize=(9, 4), tight_layout=True)
        plt.suptitle(readings[x].id +'           '+ readings[x].mouse + '\nPearsons correlation: %.3f \n' % corr, fontsize=15)
        ax = plt.subplot(121)
        plt.title('\n\n\npyLGN Model Spatial Tuning Curve')
        #plt.plot(spatial_angular_freqs, readings[x].response)
        plt.semilogx(spatial_angular_freqs, readings[x].response, '-o', marker = "None")
        plt.xticks([2**1*0.01, 2**2*0.01, 2**3*0.01,2**4*0.01, 2**5*0.01,2**6*0.01, 1.28, 2.56], ['0.02', '0.04', '0.08', '0.16', '0.32', '0.64', '1.28', '2.56'], rotation=20)
        plt.xlabel('Wavenumber\n(1/deg) ')
        plt.ylabel('Relay Response', multialignment='center')
        plt.subplot(122)
        plt.title('\n\n\nExperimental Spatial Tuning Curve')
        #plt.plot(np.squeeze(readings[x].exp_response_X), np.squeeze(readings[x].exp_response_Y))
        #plt.plot(array1, array2)
        plt.semilogx(array1, array2, '-o', marker = "None")
        plt.xticks([2**1*0.01, 2**2*0.01, 2**3*0.01,2**4*0.01, 2**5*0.01,2**6*0.01, 1.28], ['0.02', '0.04', '0.08', '0.16', '0.32', '0.64', '1.28'], rotation=20)
        plt.xlabel('Wavenumber\n(1/deg) ')
        plt.show
    return corr
def hist(readings):
    similarities = []
    for reading in readings:
        if len(reading.exp_response_X) != 0:
            similarities.append(similarity(reading, 0))
    fig, ax1 = plt.subplots(figsize=(14,5))
    ax1.set_ylabel('Quantity\n', fontsize=20)
    ax1.grid(True)
    ax1.hist(similarities, 80, facecolor='g', alpha=1, fc='k', ec='k')
    ax1.set_title('Similarity distribution\nSpatial, Non-Cortical', fontsize=30)
    ax1.set_xlabel('Pearson Correlation', fontsize=20, color='k')

    ax2 = ax1.twinx()
    plt.hist(similarities, 4, fc='skyblue', alpha=0.3)
    ax2.tick_params(axis='y', labelcolor='tab:blue')


    
#__________________________________________
#MU EXPERIMENTAL
#read lines into list
lines = []
with open('cell_responses_cut_new.csv', newline='') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count=0
    for row in csv_reader:
        if len(row) != 0:
            lines.append(row)
        line_count += 1
    print("Processed", line_count, "lines")
#initialise list for readings
readings = []
#add readings to list
line = 0
for i in range(0, 417):
    readings.append( Reading(lines[line], lines[line+1], lines[line+2]))
    line+=3
#format data
responsedata_temp = []
responsedata_string_temp = ""
for reading in readings:
    #Format Mouse
    reading.mouse = reading.mouse[1].replace('[','').replace(']','').replace('\'', '')
    #Format ID
    reading.id = reading.id[1].replace('[','').replace(']','').replace('\'', '')
    #Format Readings into float ndarray
    responsedata_string_temp = ("".join(reading.response[1])).replace('\r', '').replace('[', '').replace(']', '').replace(' 1/s', '').replace(' ', '')
    reading.response = responsedata_string_temp.split('\n')
    reading.response = np.array(reading.response, dtype=float)
#get x
network = pylgn.Network()

# create integrator
integrator = network.create_integrator(nt=5, nr=7, dt=1*pq.ms, dr=1*pq.deg)
spatial_angular_freqs = integrator.spatial_angular_freqs[:len(readings[0].response)]

x = 0
b = 0
'''
for i in range(3,8):
    #BODGED IN, used the x axis from pyLGN, put in proper x axis
    a[b].plot(spatial_angular_freqs, readings[i].response)
    b+=1
    print('correlation between 0 and ', i, "=", np.corrcoef(readings[0].response, readings[i].response))
plt.show
'''
#_________________________________
#MU
mu_x = loadmat('N:\\Data\\cells_char.mat')
mu_readings = []
for mu_cell_no in range(1190):
    mu_x_tolist = mu_x["cells_char"][mu_cell_no].tolist()
    #2. Form a temporary cell from the fields
    #TODO: simplify this?
    #if(mu_x_tolist[2]==1): #checks whether the cell was fitted with dog
    tempMuCellPy = mu_Reading()
    tempMuCellPy.mouse = mu_x_tolist[0][0]
    tempMuCellPy.id = mu_x_tolist[0][1]
    tempMuCellPy.sf_X = mu_x_tolist[0][5] #x1
    tempMuCellPy.sf_Y = mu_x_tolist[0][6] #x2
    #Add the temporary cell to the full list of cells
    mu_readings.append(tempMuCellPy)
#plt.plot(np.squeeze(mu_readings[0].sf_X), np.squeeze(mu_readings[0].sf_Y))
#plt.show
for reading in readings:
    line = 0
    found = 0
    while found == 0 and line < len(mu_readings):
        if(reading.id == mu_readings[line].id and reading.mouse == mu_readings[line].mouse):
            reading.exp_response_X = mu_readings[line].sf_X
            reading.exp_response_Y = mu_readings[line].sf_Y
            reading.valid = 1
        line += 1
        
count = 0
for reading in readings:
    #if similarity(readings[x], 0) < 1 and similarity(readings[x], 0) < -0.9:
    if len(reading.exp_response_X) > 0:
        if similarity(reading,0) > 0.5:
            count+=1
print(count)
#returns the correlation value
#pass in 1 to plot, 0 to not
similarity(readings[x],1)
hist(readings)
#Plot similarity histogram
'''
similarities = []
for reading in readings:
    if len(reading.exp_response_X) != 0:
        similarities.append(similarity(reading, 0))
fig, ax = plt.subplots(figsize=(14,5))
n, bins, patches = plt.hist(similarities, 50, density=True, facecolor='g', alpha=0.75)
'''