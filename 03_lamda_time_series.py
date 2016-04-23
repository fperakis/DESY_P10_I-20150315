'''
Code for opening and ploting the lamda detector data .nxs
Uses: nxs_handling.py
'''
# -- import packages
import numpy as np
#from nexusformat.nexus import *
from matplotlib import pyplot as plt
from nxs_handling import *

# -- folder and filenames
#base_dir = '/Volumes/Seagate Backup Plus Drive/DESY_july2015/'
#data_dir ='P10_XPCS_amorphous_ice/current/raw/'
#dataname= 'depKW_05'
base_dir = '/Volumes/Seagate Backup Plus Drive/DESY_april2016/'
data_dir ='data/raw/'
dataname= 'align'
detector='l02' # l01 fornt, l02 back
start_number,end_number = 17,616 #batch number
#it = 1

mean_photon_density=np.zeros(end_number-start_number+1)

# -- load file
file_path =format_filepath(base_dir,data_dir,dataname,detector,start_number,end_number)
data = np.array(load_nxs_data_series(file_path))
avg_data=np.zeros((data[0].shape))

for it in range(end_number-start_number+1):
    print it
    # -- threshold mask
    threshold =10 # some threshold to filter hot pixels
    data[data>threshold]=0

    # -- mean photon count
    num_pixels = len(data[it].flatten())
    num_photons = np.sum(data[it].flatten())
    mean_photon_density[it] = float(num_photons)/float(num_pixels)
    avg_data+=data[it]

avg_data/=float(end_number-start_number)

# -- photon histogram
bi,bf,db = 1,100,1
hy,hx = np.histogram(avg_data,bins = np.arange(bi,bf,db)) 

# -- Plots
plt.figure(figsize=(9,7))

# 2D plot
plt.subplot(2,1,1)
plt.imshow(avg_data,vmax= 2,interpolation='nearest')
plt.title('%s_%d'%(dataname,start_number))

# vertical average
plt.subplot(2,2,3)
plt.plot(np.mean(avg_data,axis=0))
plt.title('vertical average')
plt.ylim([0,0.1])

# histogram
plt.subplot(2,2,4)
plt.bar(hx[:-1]-db/2.,hy,width = db)
plt.title('mean intensity: %.3f photons/pixel'%(np.average(mean_photon_density)))
plt.yscale('log',nonposy='clip')
plt.xlabel('photons per pixel')
plt.ylabel('number of pixels')
plt.xlim([db/2.,5])

plt.tight_layout()

# mean photon number over 
plt.figure()
plt.plot(mean_photon_density,'go')
plt.xlabel('frame number')
plt.ylabel('mean photon density (photons/pixel')
plt.ylim([0.019,0.023])
plt.axhline(y=np.average(mean_photon_density),c='red',ls='--')
plt.title('mean intensity: %.3f photons/pixel'%(np.average(mean_photon_density)))

plt.show()


