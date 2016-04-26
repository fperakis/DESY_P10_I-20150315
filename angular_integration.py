import numpy as np
from matplotlib import pyplot as plt

# analysis tools from scikit-beam
# (https://github.com/scikit-beam/scikit-beam/tree/master/skbeam/core)
import skbeam.core.roi as roi
#import skbeam.core.correlation as corr
import skbeam.core.utils as utils
#from RingData import radial_profile

def circular_average(image, calibrated_center,mask,threshold=0, nx=100,
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
    radial_val = utils.radial_grid(calibrated_center, image.shape, pixel_size)

    # - create unitary mask (in case it's none)
    if mask is None:
        mask=np.ones(image.shape) 
  
    # - bin values of image based on the radial coordinates
    bin_edges, sums, counts = utils.bin_1D(np.ravel(radial_val*mask),
                                           np.ravel(image*mask), nx,
                                           min_x=min_x,
                                           max_x=max_x)
    th_mask = counts > threshold
    ring_averages = sums[th_mask] / counts[th_mask]
    ring_averages = sums[th_mask] / counts[th_mask]

    bin_centers = utils.bin_edges_to_centers(bin_edges)[th_mask]

    return bin_centers, ring_averages


def angular_integration_back_detector(data,n_bins=100,mask=None,center=(265,655),min_bin=40,max_bin=800,threshold=(0.001,10000)):
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
    data[data<threshold[0]]=0
    data[data>threshold[1]]=0

    # -- angular average 
    Iq= circular_average(data,calibrated_center=center,mask=mask,nx=n_bins,min_x=min_bin,max_x=max_bin)
    print q.shape,Iq[1].shape
    return Iq[0],Iq[1]



def angular_integration_front_detector(data,n_bins=150,mask=None,center=(180,2095),min_bin=600,max_bin=2000,threshold=(0.3,3)):
    ''' 
    Does the circular integration and q-calibration on the front
    detector, with an additional threshold masking
    '''


    # -- constants
    q_scale=0.97 #dirty fix - why do we need that? MUST FIX

    # ice Ih peaks position in pixel units (from fa027_05)  
    max_bins = np.array([927,990,1058,1374,1621,1762,1870,1897,1931])

    # ice Ih peaks in q (from literature)
    ice_q = iceIh_peaks()

    # -- apply threshold mask data
    data[data<threshold[0]]=0
    data[data>threshold[1]]=0

    # -- angular average 
    Iq= circular_average(data,calibrated_center=center,mask=mask,nx=n_bins,min_x=min_bin,max_x=max_bin)

    # -- calculate q-range
    dq_i=np.zeros(len(max_bins)-1)
    for i in range(len(max_bins)):
        for j in range(i+1,len(max_bins)):
             dq_i[i] = (ice_q[j]-ice_q[i])/(max_bins[j]-max_bins[i])
    q =(Iq[0]-max_bins[0])*np.mean(dq_i)*q_scale+ice_q[0]

    return q,Iq[1]


def iceIh_peaks():
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

def iceIc_peaks():
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

def diamond_peaks():

    '''
    Returns the peaks of the diamond (cvd) sustrate (in Angst-1)
    '''
    cvd_q =np.array([3.081]) 
    return cvd_q


