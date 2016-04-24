'''
Code for opening and ploting the lamda detector data .nxs
Uses: nxs_handling.py
'''
# -- import packages
import numpy as np
from matplotlib import pyplot as plt
from nxs_handling import *
from masking_lambda import *
from negative_binomial import *

# -- folder and filenames
base_dir = '/Volumes/Seagate Backup Plus Drive/DESY_april2016/'
data_dir ='data/raw/'
mask_dir = '/Users/fivos/Documents/Projects/03 XPCS of Amorhous ices/02_second_experiment/05_data_analysis/01_python/'
dataname= 'ehda_dp32'
detector='l02'

# -- parameters
n_frames = 1000#100
start_number,end_number = 8603,12202 
nx,ny = 516,1556 

# -- definition
mean_intensity=np.zeros(n_frames)

# -- mask
mask=np.array(np.load(mask_dir+'lambda_mask.npy'),dtype=bool)

# -- load file
file_path =format_filepath(base_dir,data_dir,dataname,detector,start_number,end_number)
data = load_nxs_big_data_series(file_path,np.array([n_frames,nx,ny]))

# -- mean intensity
for i in range(n_frames):
    idata=np.array(data[i,:,:])
    mean_intensity[i]=np.average(idata[np.array(mask)])#data[i,:,:])

# -- plots
plt.figure()
plt.plot(mean_intensity,'ro',alpha=0.3)
plt.xlabel('Frames')
plt.ylabel('Mean intensity (photons/pixel)')
plt.title('%s_%d_%d'%(dataname,start_number,end_number))
plt.axhline(y=np.average(mean_intensity),ls='--',c='red',label='%.4f'%np.average(mean_intensity))
plt.legend()

plt.show()

