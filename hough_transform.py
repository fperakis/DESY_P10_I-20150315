import numpy as np
import scipy
from scipy import optimize

#Gaussian peak fitting functions
Gauss_fitfunc = lambda p, x: p[0]*np.exp(-(x-p[1])**2/(2*p[2]**2))
Gauss_errfunc = lambda p, x, y: Gauss_fitfunc(p, x) - y

def calc_R(x, y, xc=0, yc=0):
    """ calculate the distance of XY maps from the center (xc, yc) """
    return np.sqrt((x-xc)**2 + (y-yc)**2)

def calc_I(r, arr, r_delta=1, r_bins=[]):
    """ calculate the radial average from distance map and intensity array """
    if len(r_bins) == 0:
        r_bins = np.arange(r.min(), r.max()+2*r_delta, r_delta) - r_delta/2
    hist_norm, r_bins = np.histogram(r, bins=r_bins)
    hist, r_bins = np.histogram(r, weights=arr, bins=r_bins)
    r_bins_center = np.array([(r_bins[i] + r_bins[i+1])/2 for i in np.arange(len(r_bins) - 1)])
    hist[np.where(hist_norm > 0)] /= hist_norm[np.where(hist_norm > 0)]
    assert hist.shape == r_bins_center.shape
    return (hist, r_bins_center)

def calc_XY(xc, yc):
    """ calculate the new XY map of each 2D points from the center (xc, yc) for the Pilatus 300k detector """
    vx = np.arange(0, 487) - xc
    vy = np.arange(0, 619) - yc
    x, y = np.meshgrid(vx, vy)
    return (x, y)

def objective_func(c, arr):
    """ optimize the peak intensity and peak width """
    x, y = calc_XY(*c)
    r = calc_R(x, y)
    I, Ir = calc_I(r, arr)
    p0 = [3.0E6, 1000, 100]
    [p1, success] = optimize.leastsq(Gauss_errfunc, p0[:], args=(Ir, I))
    if (success):
        #return p1[2]/p1[0]
        print "(%.0f, %.0f) I: %.2e, FWHM: %.2e => %.2e" % (c[0], c[1], I.max(), p1[2], p1[2]/I.max())
        return np.abs(p1[2])/I.max()
    else:
        return 1.0e300

def circular_hough(arr, x, y, r_min=0, r_max=786, r_delta=1, c_min=-10, c_max=10, c_delta=1):
    """ calculate the circular hough transform """
    print "looking for center between %.0f < c < %.0f for %.0f < r < %.0f" % (c_min, c_max, r_min, r_max)
    r_range = np.arange(r_min, r_max+r_delta, r_delta)
    #r_n = np.ceil((r_max-r_min)/r_delta) + 1
    r_n = len(r_range)
    c_range = np.arange(c_min, c_max + c_delta, c_delta)
    #c_n = np.ceil((c_max-c_min)/c_delta) + 1
    c_n = len(c_range)
    hough = np.zeros((c_n, c_n, r_n))
    print "allocating 3D hough array: %d x %d x %d = %d elements" % (c_n, c_n, r_n, hough.size)
    a_i = 0
    for a in np.arange(c_min, c_max + c_delta, c_delta):
        b_i = 0
        for b in np.arange(c_min, c_max + c_delta, c_delta):
            r = calc_R(x, y, a, b)
            r_bins = np.arange(r_min, r_max+2*r_delta, r_delta) - r_delta/2
            #hist, r_bins = np.histogram(r, weights=arr, bins=r_bins)
            hist, r_bins = calc_I(r, arr, r_bins=r_bins)
            hough[b_i, a_i, :] = hist
            b_i += 1
        a_i += 1
        print "%d/%d (%.0f%%)" % (a_i, c_n, a_i*100.0/c_n)
    hmax = hough.argmax()
    hr = hmax % 301
    hx = (hmax / 301) % 21
    hy = (hmax / 301) / 21
    print "found maximum at r = %.0f for (%.0f, %.0f)" % (r_min + hr*r_delta, c_min + hx*c_delta, c_min + hy*c_delta)
    return (c_min + hx*c_delta, c_min + hy*c_delta)

