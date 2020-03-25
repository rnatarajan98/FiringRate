#import matlab
#import matlab.engine
from scipy.io import loadmat
from cell import cell
import numpy as np
import matplotlib.pyplot as plt


#Load the matlab file into x, type dict
x = loadmat('C:\\Uni\\yr3_neuro\\Receptive_Field_Plotter\\cells.mat')

#Initialise the empty list to load cell data into
cells = []

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


#Choose which cell to plot
print("Which cell would you like to plot the receptive field of?")
print("There are currently", len(cells), "number of cells loaded")
plottedcell = int(input())

#Plot the receptive field using matplotlib
rfield = cells[plottedcell].rf_sta_lat.tolist()

#normalise the matrix
plt.subplot2grid((1,3), (0,0), rowspan=1, colspan=1)
imgdata = (rfield - np.min(rfield))/np.ptp(rfield)
plt.axis('off')
imgplot = plt.imshow(imgdata)


plt.subplot2grid((1,3), (0, 1), rowspan=1, colspan=2)
parameter = cells[plottedcell].sf_tuning[:,0]
tuning = cells[plottedcell].sf_tuning[:,1]
tuningSEM = cells[plottedcell].sf_tuning[:,2]
parameterfit = cells[plottedcell].sf_curvefit[:,0]
tuningfit = cells[plottedcell].sf_curvefit[:,1]
spont = [cells[plottedcell].Spont_mean, cells[plottedcell].Spont_mean]
spontSEM = [cells[plottedcell].Spont_SEM, cells[plottedcell].Spont_SEM]
plt.semilogx(parameterfit, tuningfit, 'xkcd:black', np.exp(-0.02 / 5.0))
plt.xlim(parameter[1], parameter[-1])
plt.xticks(parameter, parameter)
plt.grid()
(_, caps, _) = plt.errorbar(parameter.astype(float),tuning.astype(float),yerr=tuningSEM.astype(float), fmt='k.',  mfc='none', ms=5, capsize=2, elinewidth=1.5, markeredgewidth=2)
for cap in caps:
    cap.set_markeredgewidth(1)



plt.tight_layout()
plt.show()

#Plot the graph on the matlab engine
#Deprecated in favour of matplotlib

    #Start the matlab engine
    #eng = matlab.engine.start_matlab()

    #eng.workspace['fig1_data'] = matlab.double(cells[plottedcell].rf_sta_lat.tolist())
    #eng.figure(nargout=0)
    #eng.eval("imshow(rescale(fig1_data))", nargout=0)
    #eng.hold("on",nargout=0)
    #eng.box("on",nargout=0)
    #eng.hold("off",nargout=0)

    #eng.workspace['parameter'] = matlab.double(cells[plottedcell].sf_tuning[:,0].tolist())
    #eng.workspace['tuning'] = matlab.double(cells[plottedcell].sf_tuning[:,1].tolist())
    #eng.workspace['tuningSEM'] = matlab.double(cells[plottedcell].sf_tuning[:,2].tolist())
    #eng.workspace['parameterfit'] = matlab.double(cells[plottedcell].sf_curvefit[:,0].tolist())
    #eng.workspace['tuningfit'] = matlab.double(cells[plottedcell].sf_curvefit[:,1].tolist())
    #eng.workspace['spont'] = [matlab.double(cells[plottedcell].Spont_mean.tolist()), matlab.double(cells[plottedcell].Spont_mean.tolist())] 
    #eng.workspace['spontSEM'] = [matlab.double(cells[plottedcell].Spont_SEM.tolist()), matlab.double(cells[plottedcell].Spont_SEM.tolist())]
    #eng.eval("hf = figure;", nargout=0)
    #eng.eval("heb = errorbar(parameter,tuning,tuningSEM,'ok');", nargout=0)
    #eng.hold("on", nargout=0)
    #eng.eval("hp = plot(parameterfit,tuningfit,'k');", nargout=0)
    #eng.eval("h = gca;",nargout=0)
    #eng.eval("h.XScale = 'log';",nargout=0)
    #eng.eval("h.XLim = parameter([1 end]);",nargout=0)
    #eng.eval("h.XTick = parameter;",nargout=0)
    #eng.eval("h.XLabel.String = 'Spatial Frequency (cy/deg)';",nargout=0)
    #eng.eval("h.YLabel.String = 'Spikes/s';",nargout=0)
    #eng.eval("h.Layer = 'top';", nargout=0)

















plt.show()
