# -*- coding: utf-8 -*-

import wavecalc
import numpy
from wavecalc.functions import aux_rotvec as rotvec
#from wavecalc.classes import wave, surface, medium, bundle, chain




class prism:
    
    def __init__(self,medium = None, wedge = None, axis = None, everything = None):
        
        # material attribute
        #-------------------------------------------------------------------------------------
        if medium is False or everything is False:
            self.medium = None
        else:
            if isinstance(medium,wavecalc.classes.medium):
                self.medium = medium.epsilon
            elif isinstance(medium,numpy.ndarray) and numpy.shape(medium)==(3,3):
                self.medium = medium                
            elif medium is None:
                self.medium = None
            else:
                raise Exception("Material must be a wavecalc medium or (3,3) numpy ndarray")
                
            
        # wedge attribute
        #-------------------------------------------------------------------------------------            
        if wedge is False or everything is False:
            self.wedge = None
        else:
            if isinstance(wedge,(int,float)):
                self.wedge = wedge
            elif wedge is None:
                self.wedge = 0
            else:
                raise Exception("Wedge angle must be an int or a float")
                
        # axis attribute        
        #-------------------------------------------------------------------------------------
        if axis is False or everything is False:
            self.axis = None
        else:
            if axis in ['x','y','z']:
                self.axis = axis
            elif axis is None:
                self.axis = 'x'
            else:
                raise Exception("Axis must be specified as one of 'x', 'y', or 'z'")
        
        
        
        # norm1 and norm2 attributes
        #-------------------------------------------------------------------------------------
        if wedge is False or axis is False or everything is False:
            self.norm1 = None
            self.norm2 = None
        else:
            base = numpy.array([[0.,0.,1.]]).T
            if wedge is None:
                self.norm1 = base
                self.norm2 = base
            elif isinstance(wedge,(int,float)):
                self.norm1 = rotvec(base,wedge/2,self.axis)
                self.norm2 = rotvec(base,-wedge/2,self.axis)
            else:
                raise Exception("Wedge angle must be an int or a float")
            
            
                
            
        