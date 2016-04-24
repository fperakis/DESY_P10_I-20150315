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
from plot_functions import img_class
import hough_transform
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-f", "--file", action="store", type="string", dest="file",
                  help="File you wish to view", metavar="FILE.npy", default="")
parser.add_option("-d", "--data", action="store", type="string", dest="data",
                  help="Data you wish to view", metavar="water_xxxx", default="")
parser.add_option("-m", "--mask", action="store", type="string", dest="mask",
                  help="Mask to apply to the loaded file", metavar="MASK.npy", default="")
parser.add_option("-s", "--save", action="store", type="string", dest="save",
                  help="Save file to filename (will add ending)", metavar="FILE", default="")
parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
                  help="Print additional stuff", default=False)

(options, args) = parser.parse_args()

# -- folder and filenames
base_dir = '/Users/sellberg/kth/experiments/PETRA-III/P10-2016_Amorphous_ice_speckles/'
data_dir ='data/raw/'
dataname = 'water_0001'
if (options.data != ''):
    dataname = options.data
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
            
        # -- normalize
        sum_norm = sum/i
        sum_squared_norm = sum_squared/i

        # -- standard deviation
        stdev = np.sqrt(sum_squared_norm - sum_norm*sum_norm)
else:
    sum_norm = np.load(options.file)


# -- save sum
if (options.save != ''):
    np.save(options.save + '_sum.npy', sum)
    print "saved file: %s_sum.npy" % options.save
    np.save(options.save + '_sum_squared.npy', sum_squared)
    print "saved file: %s_sum_squared.npy" % options.save
    np.save(options.save + '_mean.npy', sum_norm)
    print "saved file: %s_mean.npy" % options.save

# -- mask
if (options.mask == ''):
    mask = np.ones(pilatus_shape).astype(np.bool)
    # -- threshold mask
    thresholds = [1450/604.0, 5800/604.0] # thresholds to filter shadowed/hot pixels
    mask = (mask & (sum_norm >= thresholds[0]) & (sum_norm <= thresholds[1]))
    if (options.save != ''):
        np.save(options.save + '_mask.npy', mask)
        print "saved file: %s_mask.npy" % options.save

else:
    mask = np.load(options.mask)

inverted_mask = mask.astype(np.int)*-1 + 1
if (options.file == ''):
    masked_sum = np.ma.masked_array(sum, mask=inverted_mask)
    masked_stdev = np.ma.masked_array(stdev, mask=inverted_mask)
masked_sum_norm = np.ma.masked_array(sum_norm, mask=inverted_mask)

# -- center from water fit
c_opt = (-621, 342) # (x, y) center in pixels from start of array
x, y = hough_transform.calc_XY(*c_opt)
r = hough_transform.calc_R(x, y)
I, Ir = hough_transform.calc_I(r, masked_sum_norm)

# plot stuff
if (options.file != ''):
    curr_img = img_class(options.file, masked_sum_norm, (Ir, I), c_opt)
else:
    curr_img = img_class(dataname, masked_sum_norm, (Ir, I), c_opt)
    curr_std = img_class("standard deviation", masked_stdev, (Ir, I), c_opt)
    curr_std.draw_img()
curr_img.plot_img_and_radavg()

"""
fig = plt.figure(num=None, figsize=(11.5, 5), dpi=100, facecolor='w', edgecolor='k')

canvas = fig.add_subplot(122)
canvas.set_title('radial average - (%.0f, %.0f)' % c_opt)
plt.plot(Ir, I)
#plt.plot(Ir[I.argmax()], I.max(), 'ko')
#plt.plot(Ir[I.argmax() - 100], I[I.argmax() - 100], 'ko')
#plt.plot(Ir[I.argmax() + 100], I[I.argmax() + 100], 'ko')
plt.xlabel('radius (pixels)')
plt.ylabel('mean intensity (photons/pixel/shot)')

canvas = fig.add_subplot(121)
if (options.file != ''):
    canvas.set_title(options.file + ' - mean intensity (photons/pixel/shot)')
else:
    canvas.set_title(dataname + ' - mean intensity (photons/pixel/shot)')
axes = plt.imshow(masked_sum_norm, interpolation='nearest')
colbar = plt.colorbar(axes, pad=0.01)
# approximate, center
circ = plt.Circle(c_opt, radius=Ir[I.argmax()])
circ.set_fill(False)
circ.set_edgecolor('k')
canvas.add_patch(circ)
#circ = plt.Circle(c_opt, radius=Ir[I.argmax() - 100])
#circ.set_fill(False)
#circ.set_edgecolor('k')
#canvas.add_patch(circ)
#circ = plt.Circle(c_opt, radius=Ir[I.argmax() + 100])
#circ.set_fill(False)
#circ.set_edgecolor('k')
#canvas.add_patch(circ)

if (options.file != ''):
    fig = plt.figure(num=None, figsize=(6.5, 5), dpi=100, facecolor='w', edgecolor='k')
    canvas = fig.add_subplot(111)
    canvas.set_title(dataname + ' - standard deviation (photons/pixel/shot)')
    axes = plt.imshow(masked_stdev, interpolation='nearest')
    colbar = plt.colorbar(axes, pad=0.01)
    # approximate, center
    circ = plt.Circle(c_opt, radius=Ir[I.argmax()])
    circ.set_fill(False)
    circ.set_edgecolor('k')
    canvas.add_patch(circ)

plt.show()
"""
