#!/usr/bin/env python

import os, sys
import numpy as N
import h5py as H
from  matplotlib import pyplot as P

########################################################
# Edit this variable accordingly
# Files are read for source_dir and
# written to write_dir.
# Be careful of the trailing "/";
# ensure you have the necessary read/write permissions.
########################################################
source_dir = "/Users/sellberg/kth/experiments/PETRA-III/P10-2016_Amorphous_ice_speckles/data/"
write_dir = "/Users/sellberg/kth/experiments/PETRA-III/P10-2016_Amorphous_ice_speckles/analysis/figures/"

########################################################
# Imaging class copied from Ingrid Ofte's pyana_misc code
########################################################
class img_class (object):
	"""
	Imaging class to plot 2D data interactively.
	
	Input:
	- filename
	- input array
	- input radial average (optional)
	- input center (optional)
	
	Main methods:
	- draw_img() draws 2D image without explicitly plotting it with plt.show()
	- plot_img() draws and plots 2D image
	- draw_img_and_radavg() draws 2D image and radial average without explicitly plotting it with plt.show()
	- plot_img_and_radavg() draws and plots 2D image and radial average
	- plot_instructions() plots instructions for interactive commands
	"""
	def __init__(self, filename, inarr, inrad=None, center=(0, 0)):
		self.inarr = inarr*(inarr>0)
		self.inrad = inrad
		self.center = center
		#for i in range(len(inarr)):
		#	self.inarr[i] = self.inarr[i][::-1]
		self.filename = filename
		self.cmax = self.inarr.max()
		self.cmin = self.inarr.min()
	
	def on_keypress(self,event):
		if event.key == 'p':
			if not os.path.exists(write_dir):
				os.mkdir(write_dir)
			pngtag = write_dir + "/%s.png" % (self.filename)	
			print "saving image as " + pngtag 
			P.savefig(pngtag)
		if event.key == 'r':
			colmin, colmax = self.orglims
			P.clim(colmin, colmax)
			P.draw()

	def on_click(self, event):
		if event.inaxes:
			lims = self.axes.get_clim()
			colmin = lims[0]
			colmax = lims[1]
			range = colmax - colmin
			value = colmin + event.ydata * range
			if event.button is 1 :
				if value > colmin and value < colmax :
					colmin = value
			elif event.button is 2 :
				colmin, colmax = self.orglims
			elif event.button is 3 :
				if value > colmin and value < colmax:
					colmax = value
			P.clim(colmin, colmax)
			P.draw()
				
	def draw_img(self):
		fig = P.figure(num=None, figsize=(6.5, 5), dpi=100, facecolor='w', edgecolor='k')
		cid1 = fig.canvas.mpl_connect('key_press_event', self.on_keypress)
		cid2 = fig.canvas.mpl_connect('button_press_event', self.on_click)
		canvas = fig.add_subplot(111)
		canvas.set_title(self.filename)
		self.axes = P.imshow(self.inarr, vmax = self.cmax, interpolation='nearest', origin='lower')
		self.colbar = P.colorbar(self.axes, pad=0.01)
		self.orglims = self.axes.get_clim()
		if (self.inrad):
			Ir = N.array(self.inrad[0])
			I = N.array(self.inrad[1])
			# approximate, center
			circ = P.Circle(self.center, radius=Ir[I.argmax()])
			circ.set_fill(False)
			circ.set_edgecolor('k')
			canvas.add_patch(circ)
	
	def plot_img(self):
		self.draw_img()
		self.print_instructions()
		P.show()

	def draw_img_and_radavg(self):
		fig = P.figure(num=None, figsize=(11.5, 5), dpi=100, facecolor='w', edgecolor='k')
		cid1 = fig.canvas.mpl_connect('key_press_event', self.on_keypress)
		cid2 = fig.canvas.mpl_connect('button_press_event', self.on_click)
		canvas = fig.add_subplot(121)
		canvas.set_title(self.filename)
		self.axes = P.imshow(self.inarr, vmax = self.cmax, interpolation='nearest', origin='lower')
		self.colbar = P.colorbar(self.axes, pad=0.01)
		self.orglims = self.axes.get_clim()
		if (self.inrad):
			Ir = N.array(self.inrad[0])
			I = N.array(self.inrad[1])
			# approximate, center
			circ = P.Circle(self.center, radius=Ir[I.argmax()])
			circ.set_fill(False)
			circ.set_edgecolor('k')
			canvas.add_patch(circ)
			canvas = fig.add_subplot(122)
			canvas.set_title('radial average - (%.0f, %.0f)' % self.center)
			P.plot(Ir, I)
			P.xlabel('radius (pixels)')
			P.ylabel('intensity (photons/pixel/frame)')
	
	def plot_img_and_radavg(self):
		self.draw_img_and_radavg()
		self.print_instructions()
		P.show()
	
	def print_instructions(self):
		print "Right-click on colorbar to set maximum scale."
		print "Left-click on colorbar to set minimum scale."
		print "Center-click on colorbar (or press 'r') to reset color scale."
		print "Interactive controls for zooming at the bottom of figure screen (zooming..etc)."
		print "Press 'p' to save PNG of image (with the current colorscales) in the appropriately named folder."
		print "Hit Ctl-\ or close all windows (Alt-F4) to terminate viewing program."


