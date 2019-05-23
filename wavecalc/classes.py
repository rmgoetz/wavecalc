# -*- coding: utf-8 -*-
"""
Created on Mon May 20 23:14:02 2019

@author: Ryan
"""
import wavecalc
import numpy
from numpy import random
#import wavecalc.functions as fun   # for when we define wave + surface


class wave:
    ''' A class for EM waves '''
    def __init__(self,kvec=None,efield=None,medium=None):
        
        if kvec is None:
            self.kvec = numpy.array([[0,0,1]]).T
       
        elif str(type(kvec)) == "<class 'numpy.ndarray'>" and numpy.shape(kvec) == (3,1):
            self.kvec = kvec
        
        elif kvec == 'random':
            self.kvec = random.rand(3,1)
        
        elif kvec is False:
            self.kvec = None
        
        else:
            raise Exception('Must specify kvec as (3,1) numpy.ndarray')
        
        
        
        
        if efield is None:
            self.efield = numpy.array([[1,0,0]]).T
        
        elif str(type(efield)) == "<class 'numpy.ndarray'>" and numpy.shape(efield) == (3,1):
            self.efield = efield
        
        elif efield == 'random':
            self.efield = random.rand(3,1)
        
        elif efield is False:
            self.efield = None
        
        else:
            raise Exception('Must specify efield as (3,1) numpy.ndarray')        
        
        
        
        
        if medium is None:
            self.medium = numpy.array([[1,0,0],
                                       [0,1,0],
                                       [0,0,1]])
    
        elif str(type(medium)) == "<class 'numpy.ndarray'>" and numpy.shape(medium) == (3,3):
            self.medium = medium
            
        elif medium == 'random':
            self.medium = random.rand(3,3)
            
        elif medium is False:
            self.medium = None
            
        else:
            raise Exception('Must specify medium as a (3,3) numpy.ndarray')
            
            
    def __neg__(self):
        if (isinstance(self.kvec,numpy.ndarray) 
            and numpy.shape(self.kvec) == (3,1)):
            return wave(kvec=-self.kvec,efield=self.efield,medium=self.medium)
        else:
            raise Exception("Wave must have a properly formed kvec to be negated")

        
        
    def poynting(self):
        return print('Need to add a Poynting attribute')


class surface:
    ''' A class for planar interfaces '''
    def __init__(self,normal=None,into=None,out=None):
        
        if normal is None:
            self.normal = numpy.array([[0,0,1]]).T
        
        elif str(type(normal)) == "<class 'numpy.ndarray'>" and numpy.shape(normal) == (3,1):
            self.normal = normal
            
        elif normal == 'random':
            self.normal = random.rand(3,1)
       
        elif normal is False:
            self.normal = None
        
        else:
            raise Exception('Must specify normal as (3,1) numpy.ndarray')
            
        
        
        
        if into is None:
            self.into = numpy.array([[1,0,0],
                                     [0,1,0],
                                     [0,0,1]])
        
        elif str(type(into)) == "<class 'numpy.ndarray'>" and numpy.shape(into) == (3,3):
            self.into = into
            
        elif into == 'random':
            self.into = random.rand(3,3)
        
        elif into is False:
            self.into = None
        
        else:
            raise Exception('Must specify into as a (3,3) numpy.ndarray')
            
        
        
        
        if out is None:
            self.out = numpy.array([[1,0,0],
                                    [0,1,0],
                                    [0,0,1]])
        
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
                return wavecalc.classes.surface(normal=self.normal,into=other.epsilon,out=self.out)
            else:
                raise Exception('Medium must have a properly defined dieletric tensor')
        elif isinstance(other,wavecalc.classes.wave):
            ''' Do the interaction '''
        else:
            raise Exception('Surfaces can only be added to media')
       
        
        
    def __sub__(self,other):
        if isinstance(other,wavecalc.classes.medium):
            if (isinstance(other.epsilon,numpy.ndarray)
                and numpy.shape(other.epsilon) == (3,3)):
                return wavecalc.classes.surface(normal=self.normal,into=self.into,out=other.epsilon)
            else:
                raise Exception('Medium must have a properly defined dieletric tensor')
        elif isinstance(other,wavecalc.classes.wave):
            ''' Do the interaction '''
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
            
class medium:
    ''' A class for dielectric media '''
    def __init__(self,epsilon=None):
        if epsilon is None:
            self.epsilon = numpy.array([[1,0,0],
                                        [0,1,0],
                                        [0,0,1]])
        elif str(type(epsilon)) == "<class 'numpy.ndarray'>" and numpy.shape(epsilon) == (3,3):
            self.epsilon = epsilon
        elif epsilon == 'random':
            self.epsilon = random.rand(3,3)
        elif epsilon is False:
            self.epsilon = None
        else:
            raise Exception('Must specify epsilon as a (3,3) numpy.ndarray')
            
    def __add__(self,other):
        if isinstance(other,wavecalc.classes.medium):
            if (isinstance(self.epsilon,numpy.ndarray)
                and isinstance(other.epsilon,numpy.ndarray) 
                and numpy.shape(other.epsilon) == (3,3)
                and numpy.shape(self.epsilon) == (3,3)):
                return wavecalc.classes.surface(out=self.epsilon,into=other.epsilon)
            else:
                raise Exception('Both media must have properly defined dielectric tensors')
        elif (isinstance(other,wavecalc.classes.surface) 
              and isinstance(self.epsilon,numpy.ndarray) 
              and numpy.shape(self.epsilon) == (3,3)):
            return wavecalc.classes.surface(normal=other.normal,into=other.into,out=self.epsilon)
        else:
            raise Exception('Cannot add a medium to a non-medium')
            
    def __sub__(self,other):
        if isinstance(other,wavecalc.classes.medium):
            if (isinstance(self.epsilon,numpy.ndarray)
                and isinstance(other.epsilon,numpy.ndarray) 
                and numpy.shape(other.epsilon) == (3,3)
                and numpy.shape(self.epsilon) == (3,3)):
                return wavecalc.classes.surface(out=other.epsilon,into=self.epsilon)
            else:
                raise Exception('Both media must have defined dielectric tensors')
        elif (isinstance(other,wavecalc.classes.surface) 
                and isinstance(self.epsilon,numpy.ndarray) 
                and numpy.shape(self.epsilon) == (3,3)):
            return wavecalc.classes.surface(normal=other.normal,into=self.epsilon,out=other.out)
        else:
            raise Exception('Cannot add a medium to a non-medium')
        
                
            
        
        