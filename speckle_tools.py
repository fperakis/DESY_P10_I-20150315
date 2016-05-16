#!/usr/bin/python
# -*- coding: utf-8 -*-
import numpy as np
import h5py as h
from matplotlib import pyplot as plt
import os, sys, re, time
import skbeam.core.roi as roi
import skbeam.core.correlation as corr
import skbeam.core.utils as utils
# --------------------------------------------
# speckle_tools.py
# --------------------------------------------
# Created by Jonas A. Sellberg on 2016-04-26.
# --------------------------------------------
# Collects functions for speckle analysis by
# Fivos Perakis.
# --------------------------------------------
# Last modified by JAS on 2016-04-26.
# --------------------------------------------

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-f", "--file", action="store", type="string", dest="filename", help="input file containing data you wish to analyze", metavar="FILENAME", default=None)
parser.add_option("-p", "--playlist", action="store", type="string", dest="playlist", help="name of playlist that should be exported to Spotify", metavar="PLAYLIST", default=None)
parser.add_option("-u", "--user", action="store", type="string", dest="username", help="name of Spotify user to whom the playlist should be exported", metavar="USERNAME", default=None)
parser.add_option("-v", "--verbose", action="store_true", dest="verbose", help="prints out additional information", default=False)

(options, args) = parser.parse_args()

class SpecklePattern(object):
	"""
	--------------------------------------------
	Class: SpecklePattern
	--------------------------------------------
	Makes a playlist of song titles out of a
	message read from text file or terminal.
	The playlist can be exported to Spotify.
	--------------------------------------------
	Initialization: mySpeckles = PlaylistMessage(filename, verbosity)
	--------------------------------------------
	Arguments: (optional)
	filename  - string that contains relative or absolute
                    path to text file with message
	verbosity - if true enables additional output for object
	--------------------------------------------
	Methods:   parseInput
	           makePlaylist
		   findTrack
		   findTracksToUse
		   optimizePlaylist
		   minimizeUnmatchedWords
		   exportPlaylist
	--------------------------------------------
	Instance variables:
	message     - input message stored as a list of strings
	triedTitles - set that includes all tried titles
                      that have been queried at Spotify
	foundTracks - dict of tracks that were successfully
	              queried at Spotify. The key is the track 
		      name and the value contains a nested
		      dict with all of the track information
	tracksToUse - list of strings that contains the tracks
                      to use in the playlist message in
		      sequential order
	spotify     - the spotipy.Spotify() object that
                      communicates through the metadata API
	searchLimit - maximum number of tracks from search
	              results at Spotify through the metadata API
	--------------------------------------------
	"""
	
	
	def __init__(self, data=None, mask=None, filename=None, verbose=False):
		"""
		--------------------------------------------
		Method: SpecklePattern.__init__
		--------------------------------------------
		Initializes a SpecklePattern object
		--------------------------------------------
		Arguments: (optional)
		filename  - string that contains relative or
		            absolute path to text file with
			    message
		verbosity - if true enables additional output
		            for object
		--------------------------------------------

		"""
		
		# toggles verbosiy of object
		self.verbose = verbose
		# input data from file if given
		self.data = self.parseData(filename)
		if data:
			self.data = data
		if mask:
			self.mask = mask
		else:
			# - create unitary mask (in case it's none)
			self.mask = np.ones(self.data.shape)

	
	def parseData(self, filename=None):
		"""
		--------------------------------------------
		Method: SpecklePattern.parseData
		--------------------------------------------
		Parses input for SpecklePattern
		either from text file or terminal input.
		Input is converted to UTF-8.
		--------------------------------------------
		Usage: message = mySpeckles.parseInput(filename)
		--------------------------------------------
		Arguments: (optional)
		filename - string that contains relative or
		           absolute path to text file with
			   message. If no filename is passed
			   the message will be generated from
			   the terminal prompt
		--------------------------------------------
		Returns:
		message - list of strings that contains
		          the input message
		--------------------------------------------
		"""
		
		data = []
		if filename:
			self.filename = filename
			if self.verbose:
				print "Reading data from '%s'..." % (self.filename)
                return data
	
	
	def calculateG2(self, g2q=180, center=[265,655]):
		"""
		--------------------------------------------
        	Method: SpecklePattern.calculateG2
        	--------------------------------------------
		Calculates the intensity-intensity temporal autocorrelation using
		scikit-beam packages on the back detector
		for a q-range (g2q) and width (width, num_rings,spacing)
        	--------------------------------------------
        	Usage: mySpeckles.calculateG2()
        	--------------------------------------------
                Prerequisites:
		data - data stored as numpy 3D array
                --------------------------------------------
		Arguments: (optional)
		g2q - at which q (pixels) you wanna calculate g2
		center - beam position in the data array (y, x)
		--------------------------------------------
		Returns tuple with:
		qt - q (in Ang-1) at which you calculated g2
		lag_steps[1:] - time array
		g2_avg[1:] - g2 average array
		g2_err[1:] - g2 uncertainty array
		--------------------------------------------
		"""
    		# -- parameters
    		inner_radius = g2q#180 # radius of the first ring
    		width = 1
    		spacing = 0
    		num_rings = 10
    		#center = (273, 723) # center of the speckle pattern
    		dpix = 0.055 # The physical size of the pixels
    		energy = 8.4 #keV
    		h = 4.135667516*1e-18#kev*sec
    		c = 3*1e8 #m/s
    		lambda_ = h*c/energy*1e10 # wavelength of the X-rays
    		Ldet = 5080. # # detector to sample distance
    		
    		# -- average and mask data
    		avg_data = np.average(self.data,axis=0)#*mask
    		
    		# -- ring array
    		edges = roi.ring_edges(inner_radius, width, spacing, num_rings)
    		
    		# -- convert ring to q
    		two_theta = utils.radius_to_twotheta(Ldet, edges*dpix)
    		q_val = utils.twotheta_to_q(two_theta, lambda_)
    		q_ring = np.mean(q_val, axis=1)
    		
    		# -- ring mask
    		rings = roi.rings(edges, center, self.data[0].shape)
    		ring_mask = rings*mask
    		
    		# -- calulate g2
    		num_levels = 1 #7
    		num_bufs = 100 #2
    		g2, lag_steps = corr.multi_tau_auto_corr(num_levels,num_bufs,ring_mask,self.mask*self.data)
    		
    		# -- average
    		g2_avg = np.average(g2,axis=1)
    		qt = np.average(q_ring)
    		
    		# -- standard error
    		g2_err = np.std(g2,axis=1)/np.sqrt(len(g2[0,:]))
    		
    		return qt, lag_steps[1:], g2_avg[1:], g2_err[1:]

	
	
	def radialAverage(calibrated_center, threshold=0, nx=100, pixel_size=(1, 1),  min_x=None, max_x=None):
		"""
		--------------------------------------------
        	Method: SpecklePattern.radialAverage
        	--------------------------------------------
		Radial average of the the image data
		The radial average is also known as the azimuthal integration
		(adapted from scikit-beam.roi)
        	--------------------------------------------
        	Usage: tracksToUse = mySpeckles.findTracksToUse(message)
        	--------------------------------------------
		Arguments:
    		image : array
    		    Image to compute the average as a function of radius
    		calibrated_center : tuple
    		    The center of the image in pixel units
    		    argument order should be (row, col)
    		mask  : arrayint, optional
    		    Boolean array with 1s (and 0s) to be used (or not) in the average
    		threshold : int, optional
    		    Ignore counts above `threshold`
    		    default is zero
    		nx : int, optional
    		    number of bins in x
    		    defaults is 100 bins
    		pixel_size : tuple, optional
    		    The size of a pixel (in a real unit, like mm).
    		    argument order should be (pixel_height, pixel_width)
    		    default is (1, 1)
    		min_x : float, optional number of pixels
    		    Left edge of first bin defaults to minimum value of x
    		max_x : float, optional number of pixels
    		    Right edge of last bin defaults to maximum value of x
		--------------------------------------------
		Returns:
    		bin_centers : array
    		    The center of each bin in R. shape is (nx, )
    		phi_averages : array
    		    Radial average of the image. shape is (nx, ).

		--------------------------------------------
		"""
		
    		# - create angular grid
    		phi_val = utils.angle_grid(calibrated_center, self.data.shape,pixel_size)#*180./np.pi+180.
    		
    		# - bin values of self.data based on the angular coordinates
    		bin_edges, sums, counts = utils.bin_1D(np.ravel(phi_val*self.mask),
    		                                       np.ravel(self.data*self.mask), nx,
    		                                       min_x=min_x,
    		                                       max_x=max_x)
    		th_mask = counts > threshold
    		phi_averages = sums[th_mask]/ counts[th_mask]
    		bin_centers = utils.bin_edges_to_centers(bin_edges)[th_mask]
    		
    		return bin_centers, phi_averages

	
	
        def radialIntegrationBackDetector(self, n_bins=300, center=(265,655), min_bin=-np.pi, max_bin=np.pi,threshold=(0.001,10000)):
                """
                --------------------------------------------
                Method: SpecklePattern.radialIntegrationBackDetector
                --------------------------------------------
                Wraps radialAverage to perform circular integration and q-calibration on the back
		detector, with an additional threshold masking with parameters
		for Lambda l02 detector at P10, PETRA III
                --------------------------------------------
                Usage: if mySpeckles.radialIntegrationBackDetector()
		--------------------------------------------
		Arguments (optional):
		title - string that contains track title
		        to search for, preferentially unicode
			object with UTF-8 encoding
                --------------------------------------------
                Returns: boolean value if the query was
		         successful (True) or not (False)
                --------------------------------------------
                """
		
    		# -- parameters
    		Ldet = 5080. # detector to sample distance
    		dpix = 0.055 # um (pixel size)
    		energy = 8.4 #keV
    		inner_radius, width, spacing, num_rings = 180, 120, 0, 1
    		
    		# -- constants
    		h = 4.135667516*1e-18 #kev*sec
    		c = 3*1e8 #m/s
    		lambda_ = h*c/energy*1e10 # wavelength of the X-rays
    		
    		# -- ring mask
    		edges = roi.ring_edges(inner_radius, width, spacing, num_rings)
    		ring_mask = roi.rings(edges, center, data.shape)
    		
    		# -- apply threshold mask data
    		data[data<threshold[0]]=0
    		data[data>threshold[1]]=0
    		
    		# -- radial average
    		Iphi= self.radialAverage(self.data, calibrated_center=center, mask=self.mask*ring_mask, nx=n_bins, min_x=min_bin, max_x=max_bin)
    		
    		return Iphi[0]*180./np.pi+180.,Iphi[1]

	def circular_average(calibrated_center, threshold=0, nx=100,
	                     pixel_size=(1, 1),  min_x=None, max_x=None):
	    """Circular average of the the image data
	    The circular average is also known as the radial integration
	    (adapted from scikit-beam.roi)
	    Parameters
	    ----------
	    image : array
	        Image to compute the average as a function of radius
	    calibrated_center : tuple
	        The center of the image in pixel units
	        argument order should be (row, col)
	    mask  : arrayint, optional
	        Boolean array with 1s (and 0s) to be used (or not) in the average
	    threshold : int, optional
	        Ignore counts above `threshold`
	        default is zero
	    nx : int, optional
	        number of bins in x
	        defaults is 100 bins
	    pixel_size : tuple, optional
	        The size of a pixel (in a real unit, like mm).
	        argument order should be (pixel_height, pixel_width)
	        default is (1, 1)
	    min_x : float, optional number of pixels
	        Left edge of first bin defaults to minimum value of x
	    max_x : float, optional number of pixels
	        Right edge of last bin defaults to maximum value of x
	    Returns
	    -------
	    bin_centers : array
	        The center of each bin in R. shape is (nx, )
	    ring_averages : array
	        Radial average of the image. shape is (nx, ).
	    """
	    # - create radial grid
	    radial_val = utils.radial_grid(calibrated_center, self.data.shape, pixel_size)
	
	    # - bin values of image based on the radial coordinates
	    bin_edges, sums, counts = utils.bin_1D(np.ravel(radial_val*self.mask),
	                                           np.ravel(self.data*self.mask), nx,
	                                           min_x=min_x,
	                                           max_x=max_x)
	    th_mask = counts > threshold
	    ring_averages = sums[th_mask] / counts[th_mask]
	    ring_averages = sums[th_mask] / counts[th_mask]
	
	    bin_centers = utils.bin_edges_to_centers(bin_edges)[th_mask]
	
	    return bin_centers, ring_averages
	
	
	def angularIntegrationBackDetector(n_bins=100,center=(265,655),min_bin=40,max_bin=800,threshold=(0.001,10000)):
	    ''' 
	    Does the circular integration and q-calibration on the back
	    detector, with an additional threshold masking
	    '''
	
	    # -- parameters
	    Ldet = 5080. # detector to sample distance 
	    dpix = 0.055 # um (pixel size)
	    energy = 8.4 #keV
	
	    # -- constants
	    h = 4.135667516*1e-18 #kev*sec
	    c = 3*1e8 #m/s
	    lambda_ = h*c/energy*1e10 # wavelength of the X-rays
	
	    # -- calulate q-range
	    width,spacing = float(max_bin-min_bin)/float(n_bins),0
	    edges = roi.ring_edges(min_bin, width, spacing, n_bins)
	    two_theta = utils.radius_to_twotheta(Ldet, edges*dpix)
	    q_val = utils.twotheta_to_q(two_theta, lambda_)
	    q = np.mean(q_val, axis=1)
	
	    # -- apply threshold mask data
	    self.data[data<threshold[0]]=0
	    self.data[data>threshold[1]]=0
	
	    # -- angular average 
	    Iq= circularAverage(self.data,calibrated_center=center,nx=n_bins,min_x=min_bin,max_x=max_bin)
	    print q.shape,Iq[1].shape
	    return Iq[0],Iq[1]
	
	
	
	def angularIntegrationFrontDetector(n_bins=150, center=(180,2095), min_bin=600, max_bin=2000, threshold=(0.3,3)):
	    ''' 
	    Does the circular integration and q-calibration on the front
	    detector, with an additional threshold masking
	    '''
	
	
	    # -- constants
	    q_scale=0.97 #dirty fix - why do we need that? MUST FIX
	
	    # ice Ih peaks position in pixel units (from fa027_05)  
	    max_bins = np.array([927,990,1058,1374,1621,1762,1870,1897,1931])
	
	    # ice Ih peaks in q (from literature)
	    ice_q = self.iceIhPeaks()
	
	    # -- apply threshold mask data
	    self.data[self.data<threshold[0]] = 0
	    self.data[self.data>threshold[1]] = 0
	
	    # -- angular average 
	    Iq = self.circularAverage(calibrated_center=center, nx=n_bins, min_x=min_bin, max_x=max_bin)
	
	    # -- calculate q-range
	    dq_i=np.zeros(len(max_bins)-1)
	    for i in range(len(max_bins)):
	        for j in range(i+1,len(max_bins)):
	             dq_i[i] = (ice_q[j]-ice_q[i])/(max_bins[j]-max_bins[i])
	    q =(Iq[0]-max_bins[0])*np.mean(dq_i)*q_scale+ice_q[0]
	
	    return q,Iq[1]
	
	
	def iceIhPeaks():
	    '''
	    returns the ice Ih peak (hexagonal) positions T= 98 K (in Angst-1) 
	    taken from Nature 1960, 188, 1144)
	    http://www.nature.com/nature/journal/v188/n4757/abs/1881144a0.html
	    '''
	     # ice Ih peaks positions in degrees 
	    ice_2theta =np.array([22.82,24.26,25.89,33.55,40.09,43.70,46.62,47.41,48.34,53.24])
	    xrd_energy = 8.05 #keV using a copper anode
	    h = 4.135667516*1e-18 #kev*sec
	    c = 3*1e8 #m/s
	    lambda_ = h*c/xrd_energy*1e10 #1.5498 Angst # wavelength of the X-rays
	
	    
	    # -- convert degrees to theta
	    ice_q = 4.*np.pi*np.sin(ice_2theta/2.*np.pi/180.)/lambda_
	
	    return ice_q
	
	def iceIcPeaks():
	    '''
	    returns the ice Ic (cubic) peak positions at T = 88 K (in Angst-1) 
	    taken from Nature 1960, 188, 1144)
	    http://www.nature.com/nature/journal/v188/n4757/abs/1881144a0.html
	    '''
	     # ice Ih peaks positions in degrees 
	    ice_2theta =np.array([24.26,40.11,47.43])
	    xrd_energy = 8.05 #keV at Cu K-alpha (using a copper anode)
	    h = 4.135667516*1e-18 #kev*sec
	    c = 3*1e8 #m/s
	    lambda_ = h*c/xrd_energy*1e10 #1.5498 Angst # wavelength of the X-rays
	
	    
	    # -- convert degrees to theta
	    ice_q = 4.*np.pi*np.sin(ice_2theta/2.*np.pi/180.)/lambda_
	
	    return ice_q
	
	def diamondPeaks():
	
	    '''
	    Returns the peaks of the diamond (cvd) sustrate (in Angst-1)
	    '''
	    cvd_q =np.array([3.081]) 
	    return cvd_q
	
	



if __name__ == '__main__':
	mySpeckles = SpecklePattern(filename=options.filename, verbose=options.verbose)
	#mySpeckles.calculateContrast()
	mySpeckles.calculateG2()

