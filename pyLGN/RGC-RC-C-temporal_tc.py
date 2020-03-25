# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 20:04:36 2020

@author: John
"""

import quantities as pq
import numpy as np
import matplotlib.pyplot as plt
import pylgn
import pylgn.kernels as kernel
import my_kernels as mk
from scipy.io import loadmat
from dataclasses import dataclass
import csv

@dataclass
class cell:
    name: str = None
    Spont_mean: float = None
    Spont_SEM: float = None
    rf_sta_lat: float = None
    rf_xdeg: float = None
    rf_ydeg: float = None
    rf_center: float = None
    rf_sigma: float = None
    sf_tuning: float = None
    sf_preferred: float = None
    sf_cutoff: float = None
    sf_latency: float = None
    sf_spikes: float = None
    tf_tuning: float = None
    tf_spikes: float = None
    c_tuning: float = None
    c_gain: float = None
    c_c50: float = None
    c_fit: float = None
    c_spikes: float = None
    d_tuning: float = None
    d_preferred: float = None
    d_DSI: float = None
    d_OSI: float = None
    d_spikes: float = None
    f_transient_infex: float = None
    f_spikes: float = None
    sf_curvefit: float = None
    c_curvefit: float = None

@dataclass
class mu_cell:
    mouse: str = None
    ID: str = None
    include: bool = False
    g1_xloc: float = None #x1
    g1_yloc: float = None #x2
    g2_xloc: float = None #x3
    g2_yloc: float = None #x4
    g1_xvar: float = None #x5
    g1_rot: float = None #x6 (clockwise)
    g1_yvar: float = None #x7
    g2_xvar: float = None #x8
    g2_rot: float = None #x9 (clockwise)
    g2_yvar: float = None #x10
    g1_tcen: float = None #x11 ("ie centre at x11 ms before")
    g2_tcen: float = None #x12 ("ie centre at x11 ms before")
    g1_tvar: float = None #x13
    g2_tvar: float = None #x14
    g1_tamp: float = None #x15
    g2_tamp: float = None #x16
    predicted_tfmax: float = None #save pyLGN predicted peak response temporal frequency

class Reading:
    def __init__(self, mouse, id, response):
        self.mouse = mouse
        self.id = id
        self.response = response
        self.valid = 0
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
    tf: float = 0
    tf_bandwidth: float = None

#MU pyLGN
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
    tempMuCellPy.sf_X = mu_x_tolist[0][5] #x1
    tempMuCellPy.sf_Y = mu_x_tolist[0][6] #x2
    if mu_x_tolist[0][11]>=0:
        tempMuCellPy.tf = mu_x_tolist[0][11] #peak temporal freq. response
    #Add the temporary cell to the full list of cells
    mu_readings.append(tempMuCellPy)
    
#Load the matlab file into x, type dict
mu_x = loadmat('C:\\Users\\John\Downloads\\Paul_Mu_Rar_extract\\mres_neurotech_2019_zmu-004\\mres_neurotech_2019\\data\\processed_data\\qmi_rfvs_gauss_fit-extracted.mat')

#Initialise the empty list to load cell data into
mu_cells = []

for mu_cell_no in range(1297):
    mu_x_tolist = mu_x["rfv_d2g"][mu_cell_no].tolist()
    
    #2. Form a temporary cell from the fields
    #TODO: simplify this?
    if(mu_x_tolist[2]==1): #checks whether the cell was fitted with dog
        tempMuCellPy = mu_cell()
        tempMuCellPy.mouse = mu_x_tolist[0]
        tempMuCellPy.ID = mu_x_tolist[1]
        tempMuCellPy.g1_xloc = mu_x_tolist[4] #x1
        tempMuCellPy.g1_yloc = mu_x_tolist[5] #x2
        tempMuCellPy.g2_xloc = mu_x_tolist[6] #x3
        tempMuCellPy.g2_yloc = mu_x_tolist[7] #x4
        tempMuCellPy.g1_xvar = mu_x_tolist[8] #x5
        tempMuCellPy.g1_rot = mu_x_tolist[9] #x6 (clockwise)
        tempMuCellPy.g1_yvar = mu_x_tolist[10] #x7
        tempMuCellPy.g2_xvar = mu_x_tolist[11] #x8
        tempMuCellPy.g2_rot = mu_x_tolist[12] #x9 (clockwise)
        tempMuCellPy.g2_yvar = mu_x_tolist[13] #x10
        tempMuCellPy.g1_tcen = mu_x_tolist[14] #x11 ("ie centre at x11 ms before")
        tempMuCellPy.g2_tcen = mu_x_tolist[15] #x12 ("ie centre at x11 ms before")
        tempMuCellPy.g1_tvar = mu_x_tolist[16] #x13
        tempMuCellPy.g2_tvar = mu_x_tolist[17] #x14
        tempMuCellPy.g1_tamp = mu_x_tolist[18] #x15
        tempMuCellPy.g2_tamp = mu_x_tolist[19] #x16
        
        #Add the temporary cell to the full list of cells
        mu_cells.append(tempMuCellPy)
    
#PyLGN MODEL

#cell_responses = []
#w = csv.writer(open("cell_responses_temporal.csv", "w")) #file location to save responses to

for cell_no in range((len(mu_cells))):#100):#
    #select a cell
    selected_cell_no = cell_no
    #if mu_cells[selected_cell_no].ID == ['ID:2. ch1-3']:
    X = mu_cells[selected_cell_no]

    k_max_id = 2
    feedback_strength = np.array([0.0])#, 0.5, 0.8
    
    for j, d in enumerate(feedback_strength):
        # create network
        network = pylgn.Network()
        
        # create integrator
        integrator = network.create_integrator(nt=8, nr=5, dt=0.1*pq.ms, dr=1*pq.deg)#nt=8, nr=9, dt=1*pq.ms, dr=0.1*pq.deg
        temporal_freqs = np.absolute(integrator.temporal_angular_freqs[0:66])/1000
        temporal_freqs.units = 'Hz'
        #print(temporal_freqs)
        spatial_angular_freqs = integrator.spatial_angular_freqs[:k_max_id]
        angular_freq = temporal_freqs * 2 * np.pi
        
        # create kernels
        G = mk.create_gauss_ft(A=1, a=1*pq.deg,dx=0*pq.deg, dy=0*pq.deg,xscale=X.g1_xvar,yscale=X.g1_yvar, theta=X.g1_rot)#A=1, a=1*pq.deg,dx=X.g1_xloc*pq.deg, dy=X.g1_yloc*pq.deg,xscale=X.g1_xvar,yscale=X.g1_yvar, theta=X.g1_rot
        G1 = mk.create_gauss_ft(A=1, a=1*pq.deg,dx=(X.g2_xloc-X.g1_xloc)*pq.deg, dy=(X.g2_yloc-X.g1_yloc)*pq.deg,xscale=X.g2_xvar,yscale=X.g2_yvar, theta=X.g2_rot)#A=1, a=1*pq.deg,dx=X.g2_xloc*pq.deg, dy=X.g2_yloc*pq.deg,xscale=X.g2_xvar,yscale=X.g2_yvar, theta=X.g2_rot
        T = mk.create_temp_gauss_ft(A=X.g1_tamp, a=X.g1_tvar*pq.ms,delay=X.g1_tcen*pq.ms)
        T1 = mk.create_temp_gauss_ft(A=X.g2_tamp, a=X.g2_tvar*pq.ms,delay=X.g2_tcen*pq.ms)
        Kd_t = kernel.temporal.create_delta_ft()
        Kd_s = kernel.spatial.create_delta_ft()

        # create neuron
        ganglion = network.create_ganglion_cell(kernel=(Kd_s, Kd_t))
        relay = network.create_relay_cell()
        cortical = network.create_cortical_cell()
        
        #connect neurons
        network.connect(ganglion, relay, (G, T),weight=1) #KRG
        network.connect(ganglion, relay, (G1, T1),weight=-1) #KRG
        
        wavenumber = network.integrator.spatial_angular_freqs[int(k_max_id)]
        tuning = relay.evaluate_irf_ft(angular_freq.rescale(1/pq.ms), wavenumber, 0/pq.deg)
        tuning = np.absolute(tuning)
    
    #plotting response
    print(X.mouse)
    print(X.ID)
    #plt.plot(temporal_freqs,tuning.T)
    plt.semilogx(temporal_freqs,tuning.T)
    plt.xlabel("Frequency")
    plt.ylabel("Relay Response")
    plt.show()
    
    
    #saving response to csv
    #tuning = np.array(tuning, dtype=None)
    tuning_split = np.zeros([len(temporal_freqs),1])
    for i, j in enumerate(tuning[0,:]):
        tuning_split[i,0] = j
        #print('i=',i,'j=',j)
    cell_dict = {'Mouse' : X.mouse, 'ID' : X.ID, 'response' : tuning_split}
    #for key, val in cell_dict.items():
        #w.writerow([key, val])
        
    #max response
    max_r = max(tuning[0,:]);
    for c in range(len(tuning[0,:])):
        if tuning[0,c]==max_r:
            freq_index = c
            break
    max_r_f = temporal_freqs[freq_index]
    mu_cells[selected_cell_no].predicted_tfmax = max_r_f
    #print('maximum response = ', max_r, 'at ',X.predicted_tfmax)

x = np.array([0.]*408)
y = np.array([0.]*408)
count = 0
for cell_no in range(len(mu_cells)): #count shows 407 cells show up in both mu cells and cells_char
    for char_no in range(len(mu_readings)):
        if((mu_cells[cell_no].mouse == mu_readings[char_no].mouse) and (mu_cells[cell_no].ID == mu_readings[char_no].id)):
            x[count] = mu_cells[cell_no].predicted_tfmax
            y[count] = mu_readings[char_no].tf
            count = count+1
plt.scatter(x,y)
print(np.corrcoef(x,y))