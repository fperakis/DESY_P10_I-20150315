import numpy as np
from matplotlib import pyplot as plt

# analysis tools from scikit-beam
# (https://github.com/scikit-beam/scikit-beam/tree/master/skbeam/core)
import skbeam.core.roi as roi
#import skbeam.core.correlation as corr
import skbeam.core.utils as utils
#from RingData import radial_profile

def radial_average(image, calibrated_center,mask,threshold=0, nx=100,
                     pixel_size=(1, 1),  min_x=None, max_x=None):
    """Radial average of the the image data
    The radial average is also known as the azimuthal integration
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
    phi_averages : array
        Radial average of the image. shape is (nx, ).
    """
    # - create angular grid
    phi_val = utils.angle_grid(calibrated_center, image.shape,pixel_size)#*180./np.pi+180.
    # - create unitary mask (in case it's none)
    if mask is None:
        mask=np.ones(image.shape)

    # - bin values of image based on the angular coordinates
    bin_edges, sums, counts = utils.bin_1D(np.ravel(phi_val*mask),
                                           np.ravel(image*mask), nx,
                                           min_x=min_x,
                                           max_x=max_x)
    th_mask = counts > threshold
    phi_averages = sums[th_mask]/ counts[th_mask]
    bin_centers = utils.bin_edges_to_centers(bin_edges)[th_mask]

    return bin_centers, phi_averages



def radial_integration_back_detector(data,n_bins=300,mask=None,center=(265,655),min_bin=-np.pi,max_bin=np.pi,threshold=(0.001,10000)):
    ''' 
    Does the circular integration and q-calibration on the back
    detector, with an additional threshold masking
    '''

    # -- parameters
    Ldet = 5080. # detector to sample distance 
    dpix = 0.055 # um (pixel size)
    energy = 8.4 #keV
    inner_radius,width,spacing,num_rings=180,120,0,1

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
    Iphi= radial_average(data,calibrated_center=center,mask=mask*ring_mask,nx=n_bins,min_x=min_bin,max_x=max_bin)

    return Iphi[0]*180./np.pi+180.,Iphi[1]


