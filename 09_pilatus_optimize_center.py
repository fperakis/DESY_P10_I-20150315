#!/usr/bin/env python

'''
Code for opening and plotting the pilatus 300k detector data .cbf
Uses: cbf_handling.py
'''
# -- import packages
import numpy as np
import h5py as hp
import glob as g
import sys, os, re, time
from matplotlib import pyplot as plt
from scipy import optimize
from cbf_handling import readcbf
import hough_transform
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-f", "--file", action="store", type="string", dest="file",
                  help="File you wish to view", metavar="FILE.npy", default="")
parser.add_option("-m", "--mask", action="store", type="string", dest="mask",
                  help="Mask to apply to the loaded file", metavar="MASK.npy", default="")
parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
                  help="Print additional stuff", default=False)

(options, args) = parser.parse_args()

# -- folder and filenames
base_dir = '/Users/sellberg/kth/experiments/PETRA-III/P10-2016_Amorphous_ice_speckles/'
data_dir ='data/raw/'
dataname= 'water_0001'
#base_dir = '/Volumes/Seagate Backup Plus Drive/DESY_april2016/'
#data_dir ='data/raw/'
#dataname= 'test_00041'
dataext = '.cbf'
detector='p300' # p300 front, l01 back
#start_number,end_number = 1,20 #batch number
#it = 0

filesep = '/'
#filenum = '%.5d-%.5d'%(int(start_number),int(end_number))
#filenum = '%.5d'%(int(start_number))
#filename = dataname+'_'+filenum+dataext
#filepath = base_dir+data_dir+dataname+filesep+ detector+filesep+filename
filepath = base_dir+data_dir+dataname+filesep+detector+filesep

pilatus_shape = [619, 487] # (y, x)

if (options.file == ''):
    # -- sum over all files in directory
    sum = np.zeros(pilatus_shape)
    sum_squared = np.zeros(pilatus_shape)
    files = g.glob(filepath + dataname + '*' + dataext)
    print "found %d data files in %s" % (len(files), filepath)
    print "summing images..."
    i = 0
    for f in files:
        content = readcbf(f,pilatus_shape[0],pilatus_shape[1])
        data = np.reshape(content,(pilatus_shape[0],pilatus_shape[1]))
        sum += data
        sum_squared += data*data
        i += 1
        if (options.verbose and i % 10 == 0):
            print "\tsummed %d images" % i
else:
    sum = np.load(options.file)

# -- mask
if (options.mask == ''):
    mask = np.ones(pilatus_shape).astype(np.bool)
    # -- threshold mask
    thresholds = [1450, 5800] # thresholds to filter shadowed/hot pixels
    mask = (mask & (sum >= thresholds[0]) & (sum <= thresholds[1]))
else:
    mask = np.load(options.mask)

inverted_mask = mask.astype(np.int)*-1 + 1
masked_sum = np.ma.masked_array(sum, mask=inverted_mask)

# -- center fit optimizing ring intensity and width
#c_estimate = (-800, 310)
c_estimate = (-700, 310)
#c_estimate = (0, 310)
print "original center:", c_estimate
x, y = hough_transform.calc_XY(*c_estimate)
hx, hy = hough_transform.circular_hough(masked_sum, x, y, r_min=800, r_max=1100, c_min=-100, c_max=100, c_delta=10)
#c_estimate = (c_estimate[0] - hx, c_estimate[1] - hy) # this appears to be wrong sign, since calc_XY and calc_R treats shifts the same way
c_estimate = (c_estimate[0] + hx, c_estimate[1] + hy)
print "updated center:", c_estimate
x, y = hough_transform.calc_XY(*c_estimate)
hx, hy = hough_transform.circular_hough(masked_sum, x, y, r_min=800, r_max=1100)
#c_estimate = (c_estimate[0] - hx, c_estimate[1] - hy)
c_estimate = (c_estimate[0] + hx, c_estimate[1] + hy)
print "updated center:", c_estimate
# run minimization -- downhill simplex
opt_params = optimize.fmin_powell(hough_transform.objective_func, c_estimate, args=(masked_sum,), xtol=1e-5, ftol=1e-5, disp=0)
#[c, success] = optimize.leastsq(hough_transform.objective_func, c_estimate, args=masked_sum)
#if (success):
#    print "The optimized center is: ", c
if (opt_params.any()):
    print "The optimized parameters are: ", opt_params
else:
    "The fit failed, aborting..."
    sys.exit(1)

# plot stuff
c_opt = (opt_params[0], opt_params[1])
x, y = hough_transform.calc_XY(*c_opt)
r = hough_transform.calc_R(x, y)
I, Ir = hough_transform.calc_I(r, masked_sum)

fig = plt.figure()

canvas = fig.add_subplot(122)
canvas.set_title('radial average - center (%.0f, %.0f)' % c_opt)
plt.plot(Ir, I)
plt.plot(Ir[I.argmax()], I.max(), 'ko')
plt.plot(Ir[I.argmax() - 100], I[I.argmax() - 100], 'ko')
plt.plot(Ir[I.argmax() + 100], I[I.argmax() + 100], 'ko')
plt.xlabel('radius (pixels)')
plt.ylabel('intensity (photons/pixel)')

canvas = fig.add_subplot(121)
canvas.set_title(options.file)
axes = plt.imshow(masked_sum, interpolation='nearest')
colbar = plt.colorbar(axes, pad=0.01)
# approximate, center
circ = plt.Circle(c_opt, radius=Ir[I.argmax()])
circ.set_fill(False)
circ.set_edgecolor('k')
canvas.add_patch(circ)
circ = plt.Circle(c_opt, radius=Ir[I.argmax() - 100])
circ.set_fill(False)
circ.set_edgecolor('k')
canvas.add_patch(circ)
circ = plt.Circle(c_opt, radius=Ir[I.argmax() + 100])
circ.set_fill(False)
circ.set_edgecolor('k')
canvas.add_patch(circ)

plt.show()
