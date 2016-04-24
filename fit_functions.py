import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

def fit_single_exponential(x,y,sigma=None,p0=None):
    '''
    fits a single exponential and return the fit resulting parameters and curve
    '''
    def function(x,a,b,c):
        return a*np.exp(-x*2*b)+c

    popt,pcov = curve_fit(function,x,y,sigma=sigma,p0=p0)
    curve = function(x,*popt)
    return popt,curve

def fit_double_exponential(x,y):
    '''
    fits a double exponential and return the fit resulting parameters and curve
    '''
    def function(x,a,b,c,d,e):
        return a*np.exp(-2*x*b)+c*np.exp(-2*x*d)+e

    popt,pcov = curve_fit(function,x,y)
    curve = function(x,*popt)
    return popt,curve

def fit_stretched_exponential(x,y,sigma=None,offset=1.0,p0=None):
    '''
    fits a single exponential and return the fit resulting parameters and curve
    '''
    def function(x,a,b,c,d):
        return a*np.exp(-(2*x*b)**c)+offset#d

    popt,pcov = curve_fit(function,x,y,sigma=sigma,p0=p0)
    curve = function(x,*popt)
    perr = np.sqrt(np.diag(pcov))
    return popt,curve

def fit_lorentzian(x,y):
    '''
    fits a lorentzian and return the fit resulting parameters and curve
    '''
    def function(x,a,b):
        return a/(b**(-2.)+(x-0.001)**2.)#a*np.exp(-(x*b)**c)+1#+d

    popt,pcov = curve_fit(function,x,y)
    curve = function(x,*popt)
    return popt,curve

def fit_gaussian(x,y):
    '''
    fits a lorentzian and return the fit resulting parameters and curve
    '''
    def function(x,a,b):
        return a*np.exp(-b**2*x**2.)#a*np.exp(-(x*b)**c)+1#+d

    popt,pcov = curve_fit(function,x,y)
    curve = function(x,*popt)
    return popt,curve

def fit_q2(x,y):

    '''
    fits a q^2 and return the fit resulting parameters and curve
    '''
    def function(x,a):
        return a*x**(-2.)#+b#np.exp(-b**2*x**2.)#a*np.exp(-(x*b)**c)+1#+d

    popt,pcov = curve_fit(function,x,y)
    curve = function(x,*popt)
    return popt,curve

def fit_q4(x,y):

    '''
    fits a lorentzian and return the fit resulting parameters and curve
    '''
    def function(x,a,b):
        return a*x**(-b)#+b#np.exp(-b**2*x**2.)#a*np.exp(-(x*b)**c)+1#+d

    popt,pcov = curve_fit(function,x,y)
    curve = function(x,*popt)
    return popt,curve


'''
def str_exp_fit(x,a,b,c,d):
        return c*np.exp(-a*x**b)+d#*np.exp(-d*x)+1#+c*x+d
def double_exp_fit(x,a,b,c,d):
        return c*np.exp(-x/abs(a))+d*np.exp(-x/abs(b))+offset
def exp_fit(x,a,b,c):
def linear_fit(x,a,b):
'''
