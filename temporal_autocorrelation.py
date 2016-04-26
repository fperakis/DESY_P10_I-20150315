import numpy as np
from nexusformat.nexus import *
from matplotlib import pyplot as plt
from angular_integration import *

# analysis tools from scikit-beam
# (https://github.com/scikit-beam/scikit-beam/tree/master/skbeam/core)
import skbeam.core.roi as roi
import skbeam.core.correlation as corr
import skbeam.core.utils as utils

def calculate_g2(data,mask=None,g2q=180,center = [265,655]):
    '''
    Calculates the intensity-intensity temporal autocorrelation using
    scikit-beam packages on the back detector
    for a q-range (g2q) and width (width, num_rings,spacing) 
    '''

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
    avg_data = np.average(data,axis=0)#*mask

    # -- ring array
    edges = roi.ring_edges(inner_radius, width, spacing, num_rings)

    # -- convert ring to q
    two_theta = utils.radius_to_twotheta(Ldet, edges*dpix)
    q_val = utils.twotheta_to_q(two_theta, lambda_)
    q_ring = np.mean(q_val, axis=1)
 
    # -- ring mask
    rings = roi.rings(edges, center, data[0].shape)
    ring_mask = rings*mask

    # -- calulate g2
    num_levels = 1#7
    num_bufs = 100#2
    g2, lag_steps = corr.multi_tau_auto_corr(num_levels,num_bufs,ring_mask,mask*data)
    
    # -- average
    g2_avg = np.average(g2,axis=1)
    qt = np.average(q_ring)

    # -- standard error
    g2_err = np.std(g2,axis=1)/np.sqrt(len(g2[0,:]))

    return qt,lag_steps[1:],g2_avg[1:],g2_err[1:]



"""
def g2_calculate(data,mask,gq2,snapshots_per_batch,center_x,center_y,nx=516,ny=1556,dq=10,dphi=10):
    '''
    Calculates the intensity-intensity temporal autocorrelation
    function: g2(q,t) = <I(q,t)*I(q,0)>/<I(q,0)*I(q,0)>
    '''

    rho = radial_map(center_x,center_y)
    g2_radius_min,g2_radius_max = gq2,gq2+dq
    g2_ROI = take_ring(nx,ny,rho,g2_radius_min,g2_radius_max)*mask

    g2_0_data = data[0,:,:]*g2_ROI
    g2_norm = g2_0_data*g2_0_data
    g2 = np.zeros(snapshots_per_batch)
    dg2 = np.zeros(snapshots_per_batch)

    phi = angular_map(center_x,center_y,nx,ny)
    phimax,phimin = max(phi.flatten()),min(phi.flatten())
    phi_values = np.arange(phimin,phimax,dphi)
    g2_phi = np.zeros((snapshots_per_batch,len(phi_values)))

    for it in range(snapshots_per_batch):
        g2_data = data[it,:,:]*g2_ROI*g2_0_data
        for iphi in range(len(phi_values)):
            g2_phi_ROI = take_pie_slice(nx,ny,phi,phi_values[iphi],phi_values[iphi]+dphi)*g2_ROI
            g2_phi[it,iphi] = np.average(g2_data[g2_phi_ROI>0])/np.average(g2_norm[g2_phi_ROI>0])
        g2[it],dg2[it] = np.average(g2_phi[it,:]),np.std(g2_phi[it,:])
    return g2,dg2
"""
