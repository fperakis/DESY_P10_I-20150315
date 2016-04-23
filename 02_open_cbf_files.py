'''
Code for opening and ploting the lamda detector data .nxs
Uses: nxs_handling.py
'''
# -- import packages
import numpy as np
#from nexusformat.nexus import *
from matplotlib import pyplot as plt
#from nxs_handling import *
from cbf_handling import readcbf



# -- folder and filenames
base_dir = '/Volumes/Seagate Backup Plus Drive/DESY_july2015/'
data_dir ='P10_XPCS_amorphous_ice/current/raw/'
dataname= 'aerogel1'
#base_dir = '/Volumes/Seagate Backup Plus Drive/DESY_april2016/'
#data_dir ='data/raw/'
#dataname= 'test_00041'
detector='p300' # p300 front, l01 back
start_number,end_number = 1,20 #batch number
it = 0

filesep = '/'
#filenum = '%.5d-%.5d'%(int(start_number),int(end_number))
filenum = '%.5d'%(int(start_number))
dataext = '.cbf'
filename = dataname+'_'+filenum+dataext
filepath = base_dir+data_dir+dataname+filesep+ detector+filesep+filename

'''
# -- load file
file_path =format_filepath(base_dir,data_dir,dataname,detector,start_number,end_number)
data = np.array(load_nxs_data_series(file_path))[it]
#count_time,shutter_time = load_nxs_metadata(file_path)
#print count_time, shutter_time
'''
pilatus_shape = [619,487]
content = readcbf(filepath,pilatus_shape[0],pilatus_shape[1])
data = np.reshape(content,(pilatus_shape[0],pilatus_shape[1]))

#cbf.read(filepath)
#numpy_array_with_data = content.data
#header_metadata = content.metadata

print data.shape#numpy_array_with_data

# -- threshold mask
threshold =1e5 # some threshold to filter hot pixels
data[data>threshold]=0

# -- mean photon count
num_pixels = len(data.flatten())
num_photons = np.sum(data.flatten())
mean_photon_density = float(num_photons)/float(num_pixels)

# -- photon histogram
bi,bf,db = 1,100,1
hy,hx = np.histogram(data,bins = np.arange(bi,bf,db)) 

# -- Plots
plt.figure(figsize=(12,7))

# 2D plot
plt.subplot(1,2,1)
plt.imshow(np.log10(data),vmax= 4,vmin=2.5,interpolation='nearest')
plt.title('%s_%d'%(dataname,start_number))

# vertical average
plt.subplot(2,2,2)
plt.plot(np.mean(data,axis=0))
plt.title('vertical average')
plt.ylim([0,25])

# histogram
plt.subplot(2,2,4)
plt.bar(hx[:-1]-db/2.,hy,width = db)
plt.title('mean intensity: %.3f photons/pixel'%(mean_photon_density))
plt.yscale('log',nonposy='clip')
plt.xlabel('photons per pixel')
plt.ylabel('number of pixels')
plt.xlim([db/2.,100])

plt.tight_layout()

plt.show()


