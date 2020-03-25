"""
Created on Fri Mar 13 12:36:27 2020
@author: Richard
"""
import csv
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
    stretch = 14
    for i in range(66):
        array1.append(np.squeeze(cell.exp_response_X)[(i+1)*stretch])
        array2.append(np.squeeze(cell.exp_response_Y)[(i+1)*stretch])
    corr, _ = pearsonr(cell.response, array2)
    
    if plot == 1:
        plt.figure(figsize=(9, 4), tight_layout=True)
        plt.suptitle(readings[x].id +'           '+ readings[x].mouse + '\nPearsons correlation: %.3f \n' % corr, fontsize=15)
        ax = plt.subplot(121)
        plt.title('\n\n\npyLGN Model Temporal Tuning Curve')
        plt.plot(temporal_freqs, readings[x].response)
        plt.xlabel('Frequency\n(Hz) ')
        plt.ylabel('Relay Response', multialignment='center')
        plt.subplot(122)
        plt.title('\n\n\nExperimental Temporal Tuning Curve')
        #plt.plot(np.squeeze(readings[x].exp_response_X), np.squeeze(readings[x].exp_response_Y))
        plt.plot(array1, array2)
        plt.xlabel('Frequency\n(Hz) ')
        plt.show
    return corr
def hist(readings):
    similarities = []
    for reading in readings:
        if len(np.squeeze(reading.exp_response_X)) != 0:
            similarities.append(similarity(reading, 0))
    fig, ax1 = plt.subplots(figsize=(14,5))
    ax1.hist(similarities, 50, facecolor='g', alpha=1, fc='k', ec='k')
    ax1.set_title('Similarity distribution:\nTemporal Tuning Curves with Cortical Feedback', fontsize=30)
    ax1.set_xlabel('Pearson Correlation', fontsize=20)
    ax1.set_ylabel('Quantity\n', fontsize=20)
    ax1.grid(True)

    ax2 = ax1.twinx()
    plt.hist(similarities, 4, fc='skyblue', alpha=0.3)
    ax2.tick_params(axis='y', labelcolor='tab:blue')
    plt.tight_layout()
    plt.savefig('Temporal_tc_w_Cfeedback.jpg')
    sum_sim = 0
    count = 0
    for i in range(len(similarities)):
        if similarities[i]>=0.5:
            sum_sim = sum_sim + similarities[i]
            count = count+1
    print('sum of similarities = ', sum_sim, 'summed over ', count)
    
#__________________________________________
#MU EXPERIMENTAL
#read lines into list
lines = []
with open("cell_responses_temporal.csv", newline='') as csv_file:
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
for i in range(0, 415):
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
integrator = network.create_integrator(nt=8, nr=5, dt=0.1*pq.ms, dr=1*pq.deg)
temporal_freqs = integrator.temporal_angular_freqs[0:66]
x = 0
b = 0
'''
for i in range(3,8):
    #BODGED IN, used the x axis from pyLGN, put in proper x axis
    a[b].plot(temporal_freqs, readings[i].response)
    b+=1
    print('correlation between 0 and ', i, "=", np.corrcoef(readings[0].response, readings[i].response))
plt.show
'''
#_________________________________
#MU
mu_x = loadmat('C:\\Users\\John\Downloads\\Paul_Mu_Rar_extract\\mres_neurotech_2019_zmu-004\\mres_neurotech_2019\\data\\processed_data\\cells_char.mat')
mu_readings = []
for mu_cell_no in range(1190):
    mu_x_tolist = mu_x["cells_char"][mu_cell_no].tolist()
    #2. Form a temporary cell from the fields
    #TODO: simplify this?
    #if(mu_x_tolist[2]==1): #checks whether the cell was fitted with dog
    tempMuCellPy = mu_Reading()
    tempMuCellPy.mouse = mu_x_tolist[0][0]
    tempMuCellPy.id = mu_x_tolist[0][1]
    tempMuCellPy.tf_X = mu_x_tolist[0][9]
    tempMuCellPy.tf_Y = mu_x_tolist[0][10]
    #Add the temporary cell to the full list of cells
    mu_readings.append(tempMuCellPy)
#plt.plot(np.squeeze(mu_readings[0].sf_X), np.squeeze(mu_readings[0].sf_Y))
#plt.show
for reading in readings:
    line = 0
    found = 0
    while found == 0 and line < len(mu_readings):
        if(reading.id == mu_readings[line].id and reading.mouse == mu_readings[line].mouse):
            reading.exp_response_X = mu_readings[line].tf_X
            reading.exp_response_Y = mu_readings[line].tf_Y
            reading.valid = 1
        line += 1
x=1
#returns the correlation value
#pass in 1 to plot, 0 to not
similarity(readings[0],1)
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