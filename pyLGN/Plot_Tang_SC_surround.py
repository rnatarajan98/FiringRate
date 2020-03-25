#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 09:55:45 2020

@author: shreyachawla
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 13:52:44 2020

@author: John
"""

"""
Created on Wed Jan 29 11:54:25 2020
@author: John
"""


import quantities as pq
import numpy as np
import matplotlib.pyplot as plt
import pylgn
import pylgn.kernels as kernel



k_max_id = 10
surround_strength = np.array([0, -1, -2, -3, -4, -5])
surround_width = np.array([7, 10,20,25])

C_feedback_width = np.array([5, 10, 20, 30, 40, 50, 60])
C_feedback_strength = np.array([0, 0.3, 0.4,0.5])


response = np.zeros([k_max_id, len(surround_strength)]) /pq.s

for c in C_feedback_width:#range(1):#C_feedback_width:
    for j, d in enumerate(C_feedback_strength):
        # create network
        network = pylgn.Network()
    
        # create integrator --> computes fast fourier transfotm
        integrator = network.create_integrator(nt=1, nr=10, dt=1*pq.ms, dr=0.1*pq.deg)
        
        #half t --> 
        spatial_angular_freqs = integrator.spatial_angular_freqs[:k_max_id]
    
        # create temporal kernels
        Kd_t = kernel.temporal.create_delta_ft()
        
        
        #create spatial kernels
        Kd_s = kernel.spatial.create_delta_ft()        
    
        
        Kg_s = kernel.spatial.create_gauss_ft(A=1, a= python_rfsigav[8]*pq.deg)#
        Kg_s1 = kernel.spatial.create_gauss_ft(A=1, a= 20*pq.deg)
        Kg_s2 = kernel.spatial.create_gauss_ft(A=1, a= c*pq.deg)
        
    
    
        # create neuron
        ganglion = network.create_ganglion_cell(kernel=(Kd_s, Kd_t))
        relay = network.create_relay_cell(background_response =0)
        cortical = network.create_cortical_cell()
        
        #connect neurons
        network.connect(ganglion, relay, (Kg_s, Kd_t),weight= 7) #KRG
        network.connect(ganglion, relay, (Kg_s1, Kd_t),weight= -5) #KRIG #addng weight turns this into difference of gaussian 
        network.connect(relay, cortical, (Kd_s, Kd_t), weight = 1) #
        network.connect(cortical, relay, (Kg_s2, Kd_t), weight=d) #KexRCR
        
     
        
        
    
        #Changed to full fielld grating and wave_number = k_d because it runs through 
        #different spatial frequencies
        
        for i, k_d in enumerate(spatial_angular_freqs):
            # create stimulus
            stimulus = pylgn.stimulus.create_fullfield_grating_ft(wavenumber= k_d,                                 
                                                               contrast=0.98,
                                                               orient = 0 *pq.deg)
            
            network.set_stimulus(stimulus)
            
            
            # compute
            network.compute_response(relay, recompute_ft=True)
            response[i, j] = relay.center_response[0]
        
     
    
            
    #visualize Relay cell response
    for d, R in zip(C_feedback_strength, response.T):
        print('C feedback strength=',d,'C feedback width =',c)
        plt.semilogx(spatial_angular_freqs, R, '-o', label="C_feedback_strength={}".format(d), marker = "None")
        plt.xlabel("Wavenumber (1/deg)")
        plt.ylabel("Relay Response")
        plt.legend(loc="upper right")
        #plt.axis([0,1, -1,15])
        plt.xticks([2**1*0.01, 2**2*0.01, 2**3*0.01,2**4*0.01, 2**5*0.01,2**6*0.01, 0.96], ['0.02', '0.04', '0.08', '0.16', '0.32', '0.64', '0.96'], rotation=20)
        plt.title("Cortical Feedback Width={}".format(c) + ", Surround Width = 20")
        plt.savefig("Cortical Feedback Width={}".format(c) + "_Surround Width = 20.png")
    
    plt.show()
    
    
   