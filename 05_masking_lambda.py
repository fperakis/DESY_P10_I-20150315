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

# -- folder and filenames
base_dir = '/Volumes/Seagate Backup Plus Drive/DESY_april2016/'
data_dir ='data/raw/'
mask_dir ='mask/'
save_dir = '/Users/fivos/Documents/Projects/03 XPCS of Amorhous ices/02_second_experiment/05_data_analysis/01_python/'
dataname= 'ehda_dp035'
detector='l02'

# -- parameters
n_frames = 3600#100
start_number,end_number = 3620,7219 
#n_2theta = 1#20
nx,ny = 516,1556 
min_thr,max_thr = 40,110

# -- mask
mask=np.ones((nx,ny))
mask = np.rot90(np.loadtxt(base_dir+mask_dir+'mask.txt'))
mask *= static_mask()
mask *=circular_mask(center_x=230,center_y=790,radius=630)

#all_data=np.zeros((n_2theta*n_frames,nx,ny))
sum_data = np.zeros((nx,ny))

# -- load file
file_path =format_filepath(base_dir,data_dir,dataname,detector,start_number,end_number)
data = load_nxs_big_data_series(file_path,np.array([n_frames,nx,ny]))
sum_data=np.sum(data,axis=0)
'''
for i_theta in range(n_2theta):
    print i_theta
    start_number,end_number = 1+i_theta*n_frames,(i_theta+1)*n_frames   
    file_path =format_filepath(base_dir,data_dir,dataname,detector,start_number,end_number)
    data = np.array(load_nxs_data_series(file_path))
    sum_data+=np.sum(data,axis=0)
'''
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
#plt.figure()
#plt.bar(hx[:-1]-db/2.,hy,width = db)
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


