# -*- coding: utf-8 -*-
"""
Created on Mon May 20 23:14:02 2019

@author: Ryan Goetz, ryan.m.goetz@gmail.com
"""
import wavecalc
#from wavecalc import functions
import numpy
from numpy import random
from wavecalc.functions import aux_rotate_copy as rotate_copy
from wavecalc.functions import transmit as transmit
from wavecalc.functions import reflect as reflect
#import wavecalc.functions as fun   # for when we define wave + surface


class wave:
    ''' A class for EM waves '''
    def __init__(self,kvec=None,efield=None,medium=None):
        
        if kvec is None:
            self.kvec = numpy.array([[0.,0.,1.]]).T
       
        elif str(type(kvec)) == "<class 'numpy.ndarray'>" and numpy.shape(kvec) == (3,1):
            self.kvec = kvec
        
        elif kvec == 'random':
            self.kvec = random.rand(3,1)
        
        elif kvec is False:
            self.kvec = None
        
        else:
            raise Exception('Must specify kvec as (3,1) numpy.ndarray')
        
        
        
        
        if efield is None:
            self.efield = numpy.array([[1.,0.,0.]]).T
        
        elif str(type(efield)) == "<class 'numpy.ndarray'>" and numpy.shape(efield) == (3,1):
            self.efield = efield
        
        elif efield == 'random':
            self.efield = random.rand(3,1)
        
        elif efield is False:
            self.efield = None
        
        else:
            raise Exception('Must specify efield as (3,1) numpy.ndarray')        
        
        
        
        
        if medium is None:
            self.medium = numpy.array([[1.,0.,0.],
                                       [0.,1.,0.],
                                       [0.,0.,1.]])
    
        elif str(type(medium)) == "<class 'numpy.ndarray'>" and numpy.shape(medium) == (3,3):
            self.medium = medium
            
        elif medium == 'random':
            self.medium = random.rand(3,3)
            
        elif medium is False:
            self.medium = None
            
        else:
            raise Exception('Must specify medium as a (3,3) numpy.ndarray')
            
            
            
            
    def __add__(self,other):
        if isinstance(other,wavecalc.classes.medium):
            self.medium = other.epsilon
        elif isinstance(other,wavecalc.classes.surface):
            return reflect(self,other)
        else:
            raise Exception("Waves can only add with media or surfaces")
    
        
        
    def __sub__(self,other):
        if isinstance(other,wavecalc.classes.surface):
            return transmit(self,other)
        else:
            raise Exception("Waves can only subtract surfaces")
            
            
    def __neg__(self):
        if (isinstance(self.kvec,numpy.ndarray) 
            and numpy.shape(self.kvec) == (3,1)):
            return wave(kvec=-self.kvec,efield=self.efield,medium=self.medium)
        else:
            raise Exception("Wave must have a properly formed kvec to be negated")

        
        
    def poynting(self):
        return print('Need to add a Poynting attribute')



    def rotate(self,ang,axis=None,medmove=None,verbose=None):
        rotate_copy(self,ang,axis,medmove,verbose)



class surface:
    ''' A class for planar interfaces '''
    def __init__(self,normal=None,into=None,out=None):
        
        if normal is None:
            self.normal = numpy.array([[0.,0.,1.]]).T
        
        elif str(type(normal)) == "<class 'numpy.ndarray'>" and numpy.shape(normal) == (3,1):
            self.normal = normal
            
        elif normal == 'random':
            self.normal = random.rand(3,1)
       
        elif normal is False:
            self.normal = None
        
        else:
            raise Exception('Must specify normal as (3,1) numpy.ndarray')
            
        
        
        
        if into is None:
            self.into = numpy.array([[1.,0.,0.],
                                     [0.,1.,0],
                                     [0.,0.,1.]])
        
        elif str(type(into)) == "<class 'numpy.ndarray'>" and numpy.shape(into) == (3,3):
            self.into = into
            
        elif into == 'random':
            self.into = random.rand(3,3)
        
        elif into is False:
            self.into = None
        
        else:
            raise Exception('Must specify into as a (3,3) numpy.ndarray')
            
        
        
        
        if out is None:
            self.out = numpy.array([[1.,0.,0.],
                                    [0.,1.,0.],
                                    [0.,0.,1.]])
        
        elif str(type(out)) == "<class 'numpy.ndarray'>" and numpy.shape(out) == (3,3):
            self.out = out
        
        elif out == 'random':
            self.out = random.rand(3,3)
        
        elif out is False:
            self.out = None
        
        else:
            raise Exception('Must specify out as a (3,3) numpy.ndarray')
            
            
    def __add__(self,other):
        if isinstance(other,wavecalc.classes.medium):
            if (isinstance(other.epsilon,numpy.ndarray)
                and numpy.shape(other.epsilon) == (3,3)):
                return surface(normal=self.normal,into=other.epsilon,out=self.out)
            else:
                raise Exception('Medium must have a properly defined dieletric tensor')
        elif isinstance(other,wavecalc.classes.wave):
            return transmit(other,self)
        else:
            raise Exception('Surfaces can only add media or waves')
       
        
        
    def __sub__(self,other):
        if isinstance(other,wavecalc.classes.medium):
            if (isinstance(other.epsilon,numpy.ndarray)
                and numpy.shape(other.epsilon) == (3,3)):
                return surface(normal=self.normal,into=self.into,out=other.epsilon)
            else:
                raise Exception('Medium must have a properly defined dieletric tensor')
        elif isinstance(other,wavecalc.classes.wave):
            return reflect(other,self)
        else:
            raise Exception('Surfaces can only subtract waves or media')
            
    
    def __neg__(self):
        if (isinstance(self.normal,numpy.ndarray) 
            and numpy.shape(self.normal) == (3,1)):
            return surface(normal=-self.normal,into=self.into,out=self.out)
        else:
            raise Exception("Surface must have a ")
            
    
    
    def __invert__(self):
        return surface(normal=self.normal,into=self.out,out=self.into)
    
    
    
    def rotate(self,ang,axis=None,medmove=None,verbose=None):
        rotate_copy(self,ang,axis,medmove,verbose)
        
            
        
        
        
        
        
        
        
        
class medium:
    ''' A class for dielectric media '''
    def __init__(self,epx=None,epy=None,epz=None,epsilon=None):
        if (isinstance(epx,complex) 
            or isinstance(epy,complex) 
            or isinstance(epz,complex)):
            dat = complex
        else:
            dat = float
        if epsilon is None:
            self.epsilon = numpy.array([[1.,0.,0.],
                                        [0.,1.,0.],
                                        [0.,0.,1.]],dtype=dat)
            if (isinstance(epx,int) 
                or isinstance(epx,float) 
                or isinstance(epx,complex)):
                self.epsilon[0,0] = epx
            if (isinstance(epy,int) 
                or isinstance(epy,float) 
                or isinstance(epy,complex)):
                self.epsilon[1,1] = epy
            if (isinstance(epz,int) 
                or isinstance(epz,float) 
                or isinstance(epz,complex)):
                self.epsilon[2,2] = epz
        elif str(type(epsilon)) == "<class 'numpy.ndarray'>" and numpy.shape(epsilon) == (3,3):
            self.epsilon = epsilon
        elif epsilon == 'random':
            self.epsilon = random.rand(3,3)
        elif epsilon is False:
            self.epsilon = None
        else:
            raise Exception('Must specify epsilon as a (3,3) numpy.ndarray')
    

    def epx(self,epsilon_xx):
        if (not isinstance(self.epsilon,numpy.ndarray) 
            or not numpy.shape(self.epsilon)==(3,3)):
            raise Exception("Medium must have proper epsilon attribute")
        if (isinstance(epsilon_xx,int) 
            or isinstance(epsilon_xx,float) 
            or isinstance(epsilon_xx,complex)):
            self.epsilon[0,0] = epsilon_xx
        else:
            raise Exception("Epsilon_xx must be an int, float, or complex")



    def epy(self,epsilon_yy):
        if (not isinstance(self.epsilon,numpy.ndarray) 
            or not numpy.shape(self.epsilon)==(3,3)):
            raise Exception("Medium must have proper epsilon attribute")
        if (isinstance(epsilon_yy,int) 
            or isinstance(epsilon_yy,float) 
            or isinstance(epsilon_yy,complex)):
            self.epsilon[1,1] = epsilon_yy
        else:
            raise Exception("Epsilon_yy must be an int, float, or complex")



    def epz(self,epsilon_zz):
        if (not isinstance(self.epsilon,numpy.ndarray) 
            or not numpy.shape(self.epsilon)==(3,3)):
            raise Exception("Medium must have proper epsilon attribute")
        if (isinstance(epsilon_zz,int) 
            or isinstance(epsilon_zz,float) 
            or isinstance(epsilon_zz,complex)):
            self.epsilon[2,2] = epsilon_zz
        else:
            raise Exception("Epsilon_zz must be an int, float, or complex")








        
    def __add__(self,other):
        if isinstance(other,wavecalc.classes.medium):
            if (isinstance(self.epsilon,numpy.ndarray)
                and isinstance(other.epsilon,numpy.ndarray) 
                and numpy.shape(other.epsilon) == (3,3)
                and numpy.shape(self.epsilon) == (3,3)):
                return surface(out=self.epsilon,into=other.epsilon)
            else:
                raise Exception('Both media must have properly defined dielectric tensors')
        elif (isinstance(other,wavecalc.classes.surface) 
              and isinstance(self.epsilon,numpy.ndarray) 
              and numpy.shape(self.epsilon) == (3,3)):
            return surface(normal=other.normal,into=other.into,out=self.epsilon)
        else:
            raise Exception('Cannot add a medium to a non-medium')
            
    def __sub__(self,other):
        if isinstance(other,wavecalc.classes.medium):
            if (isinstance(self.epsilon,numpy.ndarray)
                and isinstance(other.epsilon,numpy.ndarray) 
                and numpy.shape(other.epsilon) == (3,3)
                and numpy.shape(self.epsilon) == (3,3)):
                return surface(out=other.epsilon,into=self.epsilon)
            else:
                raise Exception('Both media must have defined dielectric tensors')
        elif (isinstance(other,wavecalc.classes.surface) 
                and isinstance(self.epsilon,numpy.ndarray) 
                and numpy.shape(self.epsilon) == (3,3)):
            return surface(normal=other.normal,into=self.epsilon,out=other.out)
        else:
            raise Exception('Cannot add a medium to a non-medium')
        
                
    def rotate(self,ang,axis=None,medmove=None,verbose=None):
        rotate_copy(self,ang,axis,medmove,verbose)
        

    

            
        
        