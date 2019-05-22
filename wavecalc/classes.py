# -*- coding: utf-8 -*-
"""
Created on Mon May 20 23:14:02 2019

@author: Ryan
"""
import wavecalc
import numpy
from numpy import random


class wave:
    ''' A class for EM waves '''
    def __init__(self,kvec=None,efield=None,medium=None):
        
        if kvec is None:
            self.kvec = numpy.array([0,0,1])
       
        elif str(type(kvec)) == "<class 'numpy.ndarray'>" and numpy.shape(kvec) == (3,):
            self.kvec = kvec
        
        elif kvec == 'random':
            self.kvec = random.rand(3,)
        
        elif kvec is False:
            self.kvec = None
        
        else:
            raise Exception('Must specify kvec as (3,) numpy.ndarray')
        
        
        
        
        if efield is None:
            self.efield = numpy.array([1,0,0])
        
        elif str(type(efield)) == "<class 'numpy.ndarray'>" and numpy.shape(efield) == (3,):
            self.efield = efield
        
        elif efield == 'random':
            self.efield = random.rand(3,)
        
        elif efield is False:
            self.efield = None
        
        else:
            raise Exception('Must specify efield as (3,) numpy.ndarray')        
        
        
        
        
        if medium is None:
            self.medium = numpy.matrix([[1,0,0],
                                     [0,1,0],
                                     [0,0,1]])
    
        elif str(type(medium)) == "<class 'numpy.matrix'>" and numpy.shape(medium) == (3,3):
            self.medium = medium
            
        elif medium == 'random':
            self.medium = numpy.asmatrix(random.rand(3,3))
            
        elif medium is False:
            self.medium = None
            
        else:
            raise Exception('Must specify medium as a (3,3) numpy.matrix')

        
        
    def poynting(self):
        return print('Need to add a Poynting attribute')


class surface:
    ''' A class for planar interfaces '''
    def __init__(self,normal=None,into=None,out=None):
        
        if normal is None:
            self.normal = numpy.array([0,0,1])
        
        elif str(type(normal)) == "<class 'numpy.ndarray'>" and numpy.shape(normal) == (3,):
            self.normal = normal
            
        elif normal == 'random':
            self.normal = random.rand(3,)
       
        elif normal is False:
            self.normal = None
        
        else:
            raise Exception('Must specify normal as (3,) numpy.ndarray')
            
        
        
        
        if into is None:
            self.into = numpy.matrix([[1,0,0],
                                     [0,1,0],
                                     [0,0,1]])
        
        elif str(type(into)) == "<class 'numpy.matrix'>" and numpy.shape(into) == (3,3):
            self.into = into
            
        elif into == 'random':
            self.into = numpy.asmatrix(random.rand(3,3))
        
        elif into is False:
            self.into = None
        
        else:
            raise Exception('Must specify into as a (3,3) numpy.matrix')
            
        
        
        
        if out is None:
            self.out = numpy.matrix([[1,0,0],
                                     [0,1,0],
                                     [0,0,1]])
        
        elif str(type(out)) == "<class 'numpy.matrix'>" and numpy.shape(out) == (3,3):
            self.out = out
        
        elif out == 'random':
            self.out = numpy.asmatrix(random.rand(3,3))
        
        elif out is False:
            self.out = None
        
        else:
            raise Exception('Must specify out as a (3,3) numpy.matrix')
            
            
class medium:
    ''' A class for dielectric media '''
    def __init__(self,epsilon=None):
        if epsilon is None:
            self.epsilon = numpy.matrix([[1,0,0],
                                     [0,1,0],
                                     [0,0,1]])
        elif str(type(epsilon)) == "<class 'numpy.matrix'>" and numpy.shape(epsilon) == (3,3):
            self.epsilon = epsilon
        elif epsilon == 'random':
            self.epsilon = numpy.asmatrix(random.rand(3,3))
        elif epsilon is False:
            self.epsilon = None
        else:
            raise Exception('Must specify epsilon as a (3,3) numpy.matrix')
            
    def __add__(self,other):
        if isinstance(other,wavecalc.classes.medium):
            if (isinstance(self.epsilon,numpy.matrix)
                and isinstance(other.epsilon,numpy.matrix) 
                and numpy.shape(other.epsilon) == (3,3)
                and numpy.shape(self.epsilon) == (3,3)):
                return wavecalc.classes.surface(out=self.epsilon,into=other.epsilon)
            else:
                raise Exception('Both media must have defined dielectric tensors')
        else:
            raise Exception('Cannot add a medium to a non-medium')
            
    def __sub__(self,other):
        if isinstance(other,wavecalc.classes.medium):
            if (isinstance(self.epsilon,numpy.matrix)
                and isinstance(other.epsilon,numpy.matrix) 
                and numpy.shape(other.epsilon) == (3,3)
                and numpy.shape(self.epsilon) == (3,3)):
                return wavecalc.classes.surface(out=other.epsilon,into=self.epsilon)
            else:
                raise Exception('Both media must have defined dielectric tensors')
        else:
            raise Exception('Cannot add a medium to a non-medium')
        
                
            
        
        