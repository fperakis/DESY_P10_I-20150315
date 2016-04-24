'''
Code for opening and ploting the lamda detector data .nxs
Uses: nxs_handling.py
'''
# -- import packages
import numpy as np
#from nexusformat.nexus import *
from matplotlib import pyplot as plt
from nxs_handling import *
from masking_lambda import *
from negative_binomial import *

# -- folder and filenames
base_dir = '/Volumes/Seagate Backup Plus Drive/DESY_april2016/'
data_dir ='data/raw/'
mask_dir = '/Users/fivos/Documents/Projects/03 XPCS of Amorhous ices/02_second_experiment/05_data_analysis/01_python/'
dataname= 'ehda_dp035'
detector='l02'

# -- parameters
n_frames = 1000#100
start_number,end_number = 7220,8219#1001,4600#7601,8600#4601,5600#12203,15802#3620,7219 
nx,ny = 516,1556 
max_photon_count=6

# -- definition
mean_intensity=np.zeros(n_frames)
hist = np.zeros((n_frames,max_photon_count-1))
contrast = np.zeros(n_frames)
# -- mask
mask = np.array(np.load(mask_dir+'lambda_mask.npy'),dtype=bool)

# -- load file
file_path =format_filepath(base_dir,data_dir,dataname,detector,start_number,end_number)
data = load_nxs_big_data_series(file_path,np.array([n_frames,nx,ny]))

# -- mean intensity
for i in range(n_frames):
    print i
    idata = data[i,:,:]
    mean_intensity[i] = np.average(idata[mask])
    # photon histogram
    bi,bf,db = 0,max_photon_count,1
    hy,hx = np.histogram(idata[mask],bins = np.arange(bi,bf,db))
    contrast[i] = 1./fit_negative_binomial_scipy(hx[:-1],hy,np.average(mean_intensity),p0=10)
    #contrast[i]=1./nb_par

# -- plots
plt.figure()
plt.subplot(1,2,1)
plt.plot(mean_intensity,'ro',alpha=0.3)
plt.xlabel('Frames')
plt.ylabel('Mean intensity (photons/pixel)')
plt.title('%s_%d_%d'%(dataname,start_number,end_number))
plt.axhline(y=np.average(mean_intensity),ls='--',c='red',label='%.2f'%np.average(mean_intensity))
plt.legend()

plt.subplot(1,2,2)
plt.plot(mean_intensity,contrast,'bo',alpha=0.5,label='%.2f'%np.average(contrast[contrast>0]))
plt.legend()
plt.ylim([0,1])
plt.show()

'''
#sum_data=np.sum(data,axis=0)

# -- photon histogram
bi,bf,db = 1,1e3,1#max(sum_data.flatten()),1
hy,hx = np.histogram(sum_data,bins = np.arange(bi,bf,db)) 

# -- apply mask
mask[sum_data>max_thr]=0
mask[sum_data<min_thr]=0
sum_data*=mask

# -- 2d plots
plt.figure(figsize=(12,5))
plt.subplot(2,2,1)
plt.title('masked:%.2f percent'% ((1-float(len(mask[mask==1]))/float(len(mask.flatten())))*100))
plt.imshow(sum_data,interpolation='nearest',vmax=100,vmin=50)

plt.subplot(2,2,3)
plt.imshow(mask,interpolation='nearest')

# -- photon histograms
plt.subplot(1,2,2)
plt.title('min=%.2f,max=%.2f'%(min_thr,max_thr))
plt.plot(hx[:-1],hy,'r-')
plt.yscale('log',nonposy='clip')
plt.xscale('log',nonposy='clip')
plt.xlabel('Intensity (photons)')
plt.ylabel('Number of pixels')
plt.xlim([1e0,1e3])
plt.axvline(x=max_thr,ls='--',c='gray')
plt.axvline(x=min_thr,ls='--',c='gray')

plt.show()
np.save(save_dir+'lambda_mask.npy',mask)
'''

