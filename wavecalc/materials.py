# -*- coding: utf-8 -*-

import wavecalc
import numpy
from wavecalc.classes import medium


'''

A library of common materials

Dielectric properties are given for 1064 nm light

'''

def sellmeier1(lamb,A1,A2,A3,B2,B3):
    return A1+(A2/(lamb**2 - B2))+(A3/(lamb**2 - B3))




def ktp(lamb = None):
    if lamb is None:
        lamb = 1.064
    elif isinstance(lamb,(int,float)):
        lamb = lamb/1000.   
    else:
        raise Exception("Wavelength must be specified as an int or a float")
    
    '''
    Source for Sellmeier coefficients: K. Kato and E. Takaoka, 'Sellmeier and thermo-optic dispersion formulas for KTP,' Applied Optics, Vol. 41, No. 24, 20 Aug 2002
    '''
    epx = sellmeier1(lamb,3.291,0.0414,9.35522,0.03978,31.45571)
    epy = sellmeier1(lamb,3.45018,0.04341,16.98825,0.04597,39.43799)
    epz = sellmeier1(lamb,4.59423,0.06206,110.80672,0.04763,86.12171)
    ep = numpy.array([[epx,0,0],[0,epy,0],[0,0,epz]])
    return medium(epsilon=ep)
