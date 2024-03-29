# -*- coding: utf-8 -*-
import wavecalc
import numpy
import copy
from numpy import random
from collections import UserList
from wavecalc.aux_funcs import aux_rotate_copy as rotate_copy
from wavecalc.main_funcs import transmit as transmit
from wavecalc.main_funcs import reflect as reflect
from wavecalc.main_funcs import crash as crash
from wavecalc.aux_funcs import aux_goodtest as goodtest
from wavecalc.aux_funcs import aux_fixmode, aux_clean, aux_check_ab
"""
Created on Mon May 20 23:14:02 2019

@author: Ryan Goetz, ryan.m.goetz@gmail.com
last update: April 10, 2020
"""
'''
Table of Contents

    wave Class -------------------------------- Line 116
       
        __init__ ------------------------------ Line 118
        pol ----------------------------------- Line 228
        amp ----------------------------------- Line 251
        poynting ------------------------------ Line
        rotate -------------------------------- Line
        fixmode ------------------------------- Line
        index --------------------------------- Line
        clean --------------------------------- Line
        k ------------------------------------- Line
        propagate ----------------------------- Line
        jones --------------------------------- Line
        __add__ ------------------------------- Line
        __sub__ ------------------------------- Line
        __pos__ ------------------------------- Line
        __neg__ ------------------------------- Line
        __invert__ ---------------------------- Line
        __matmul__ ---------------------------- Line
        __eq__ -------------------------------- Line
        
       
        
    surface Class ----------------------------- Line 508
    
        __init__ ------------------------------ Line
        rotate -------------------------------- Line
        clean --------------------------------- Line
        __add__ ------------------------------- Line
        __sub__ ------------------------------- Line
        __pos__ ------------------------------- Line
        __neg__ ------------------------------- Line
        __matmul__ ---------------------------- Line
        __invert__ ---------------------------- Line
        __eq__ -------------------------------- Line
        
        
        
    medium Class ------------------------------ Line 724

        __init__ ------------------------------ Line
        epx ----------------------------------- Line
        epy ----------------------------------- Line
        epz ----------------------------------- Line
        ep_all -------------------------------- Line
        rotate -------------------------------- Line
        clean --------------------------------- Line
        __add__ ------------------------------- Line
        __pos__ ------------------------------- Line
        __sub__ ------------------------------- Line
        __eq__ -------------------------------- Line



    bundle Class ------------------------------ Line 941
    
        __init__ ------------------------------ Line
        append -------------------------------- Line
        __setitem__ --------------------------- Line
        
        
    
    chain Class ------------------------------- Line 1010
    
        __init__ ------------------------------ Line
        append -------------------------------- Line
        __setitem__ --------------------------- Line
        
        
        
    example_class Class ----------------------- Line 1080
    
        __init__ ------------------------------ Line
        


Last line check June 8, 2019

'''

####################################################################################
#----------------------------------------------------------------------------------#
#--w---------------------------w------w------w-----------------w------wwwwwwwwwww--#
#---w-------------------------w------w-w------w---------------w------w-------------#
#----w-----------------------w------w---w------w-------------w------w--------------#
#-----w---------------------w------w-----w------w-----------w------w---------------#
#------w---------w---------w------w-------w------w---------w------wwwwwwwwwww------#
#-------w-------w-w-------w------wwwwwwwwwww------w-------w------w-----------------#
#--------w-----w---w-----w------w-----------w------w-----w------w------------------#
#---------w---w-----w---w------w-------------w------w---w------w-------------------#
#----------w-w-------w-w------w---------------w------w-w------w--------------------#
#-----------w---------w------w-----------------w------w------wwwwwwwwwww-----------#
#----------------------------------------------------------------------------------#
####################################################################################
class wave:
    ''' A class for EM waves '''
    def __init__(self,kvec=None,efield=None,phase=None,**kwargs):
        
        medium = kwargs.pop('medium',None)
        pol = kwargs.pop('pol',None)
        amp = kwargs.pop('amp',None)
        ab = kwargs.pop('ab',None)
        fix = kwargs.pop('fix',None)
        everything = kwargs.pop('everything',None)
        
        if len(kwargs) != 0:
            print("Undefined arguments passed to wave.__init__")
        
        # kvec attribute
        #-------------------------------------------------------------------------------------        
        if kvec is False or everything is False:
            self.kvec = None
        
        elif kvec is None:
            self.kvec = numpy.array([[0.,0.,1.]]).T
       
        elif type(kvec) is numpy.ndarray and numpy.shape(kvec) == (3,1):
            self.kvec = kvec
            fix = False
        
        elif kvec == 'random':
            self.kvec = random.rand(3,1)
        
        else:
            raise Exception("Must specify kvec as (3,1) numpy.ndarray, 'random', False, or None")
        
        
        # efield attribute
        #-------------------------------------------------------------------------------------
        if efield is False or everything is False:
            self.efield = None
            
        elif efield is None:
            if type(pol) is str:
                if isinstance(amp,(int,float)):
                    multiplier = amp
                else:
                    multiplier = 1
                if pol == 'x':
                    self.efield = multiplier*numpy.array([[1.,0.,0.]]).T
                elif pol == 'y':
                    self.efield = multiplier*numpy.array([[0.,1.,0.]]).T
                elif pol == 'z':
                    self.efield = multiplier*numpy.array([[0.,0.,1.]]).T
            elif type(pol) is numpy.ndarray and numpy.shape(pol) == (3,1):
                sqr_mod = numpy.conj(pol).T @ pol
                normpol = pol/numpy.sqrt(sqr_mod)
                if isinstance(amp,(int,float)):
                    self.efield = amp*normpol
                else:
                    self.efield = normpol
            elif kvec is None and isinstance(amp,(int,float,complex)):
                self.efield = numpy.array([[amp,0.,0.]]).T
            elif kvec is None:
                self.efield = numpy.array([[1.,0.,0.]]).T
            else:
                self.efield = None
                
        elif type(efield) is numpy.ndarray and numpy.shape(efield) == (3,1):
            self.efield = efield
        
        elif efield == 'random':
            self.efield = random.rand(3,1)
        
        else:
            raise Exception("Must specify efield as (3,1) numpy.ndarray, 'random', False, or None")        
        
        
        # phase attribute
        #-------------------------------------------------------------------------------------
        if phase is False or everything is False:
            self.phase = None
            
        elif phase is None:
            self.phase = 0.
       
        elif isinstance(phase,(int,float,complex)):
            self.phase = phase
        
        elif phase == 'random':
            self.phase = 2*numpy.pi*random.rand()
        
        else:
            raise Exception("Must specify phase as an int, float, complex, 'random', False, or None")
        
        
        # medium attribute
        #-------------------------------------------------------------------------------------
        if medium is False or everything is False:
            self.medium = None
        
        if medium is None:
            self.medium = numpy.array([[1.,0.,0.],
                                       [0.,1.,0.],
                                       [0.,0.,1.]])
    
        elif type(medium) is numpy.ndarray and numpy.shape(medium) == (3,3):
            self.medium = medium
            if fix is not False:
                self.fixmode(ab=ab)
            
        elif medium == 'random':
            self.medium = random.rand(3,3)
            if fix is not False:
                self.fixmode(ab=ab)
            
        else:
            raise Exception("Must specify medium as a (3,3) numpy.ndarray, 'random', False, or None")
            
                    
    
    # Custom methods            
    #-----------------------------------------------------------------------------------------
    def pol(self,polar=None):
        if polar is None:
            if type(self.efield) is numpy.ndarray and numpy.shape(self.efield) == (3,1):
                sqr_mod = (numpy.conj(self.efield).T @ self.efield)[0,0].real
                norm = numpy.sqrt(sqr_mod)
                if norm == 0:
                    normpol = 0*self.efield
                else:
                    normpol = self.efield/norm
                return normpol
            else:
                print("No electric field found")
        elif type(polar) is numpy.ndarray and numpy.shape(polar) == (3,1):
            sqr_mod = (numpy.conj(self.efield).T @ self.efield)[0,0].real
            norm = numpy.sqrt(sqr_mod)
            sqr_mod_pol = (numpy.conj(polar).T @ polar)[0,0].real
            norm_pol = numpy.sqrt(sqr_mod_pol)
            self.efield = (norm/norm_pol)*polar
        else:
            raise Exception("New polarization must be a (3, 1) numpy.ndarray")
        
        
        
    def amp(self,amplitude=None,rel=False):
        if amplitude is None:
            if type(self.efield) is numpy.ndarray and numpy.shape(self.efield) == (3,1):
                sqr_mod = (numpy.conj(self.efield).T @ self.efield)[0,0]
                norm = numpy.sqrt(sqr_mod)
                return norm.real    
            else:
                print("No appropriate electric field found")
        elif isinstance(amplitude,(int,float)):
            if rel:
                self.efield = amplitude*self.efield
            else:
                sqr_mod = (numpy.conj(self.efield).T @ self.efield)[0,0]
                norm = numpy.sqrt(sqr_mod)
                norm = norm.real
                self.efield = (amplitude/norm)*self.efield 
        else:
            raise Exception("New amplitude must be an int or a float")
            
     
        
    def poynting(self,**kwargs):
        
        scale = kwargs.pop('scale',None)
        norm = kwargs.pop('norm',None)
        
        if len(kwargs) != 0:
            print("Undefined argument passed to wave.poynting")
        
        if goodtest(self,test_type='wave'):
            if scale is None or not isinstance(scale,(int,float)):
                scale = 1
            ee = self.efield
            kk = self.kvec
            hh = numpy.cross(kk.T,ee.T).T 
            hhc = numpy.conj(hh)
            S = 0.5*scale*numpy.cross(ee.T,hhc.T).T
            if norm:
                Snorm = numpy.sqrt(numpy.conj(S.T) @ S)[0,0]
                S = S/Snorm    
            return S
        else:
            raise Exception("Wave object is not well-formed, check for improper attributes")



    def rotate(self,ang,axis,**kwargs):
        
        medmove = kwargs.pop('medmove',None)
        fix = kwargs.pop('fix',None)
        verbose = kwargs.pop('verbose',None)
        
        if len(kwargs) != 0:
            print("Undefined argument passed to wave.rotate")
        
        if goodtest(self,test_type='wave'):
            medQ = isinstance(self.medium,numpy.ndarray)
            if fix and medQ:
                back = aux_check_ab(self)
            rotate_copy(self,ang,axis,medmove,verbose)
            if fix and medQ:
                self.fixmode(ab=back,conserve=True)
        else:
            raise Exception("Wave object is not well-formed, check for improper attributes")
        
    
    
    def fixmode(self,**kwargs):
        
        ab = kwargs.pop('ab',None)
        k0 = kwargs.pop('k0',None)
        conserve = kwargs.pop('conserve',None)
        verbose = kwargs.pop('verbose',None)
        
        if len(kwargs) != 0:
            print("Undefined argument passed to wave.fixmode")
        
        if goodtest(self,test_type='wave'):
            back = aux_fixmode(wave=self,ab=ab,k0=k0,conserve=conserve,verbose=verbose)
            self.kvec = back[0]
            self.efield = back[1]
        else:
            raise Exception("Wave object is not well-formed, check for improper attributes")
            
            
            
    def index(self,k0=None):
        if goodtest(self,test_type='wave'):
            if k0 is None:
                k0 = 1
            if numpy.shape(self.kvec) == (3,1):
                kre = self.kvec.real
                kim = self.kvec.imag
                Nre = numpy.sqrt((kre.T @ kre)[0,0])/k0
                Nim = numpy.sqrt((kim.T @ kim)[0,0])/k0
                return [Nre,Nim]
            else:
                raise Exception("No wave vector found")
        else:
            raise Exception("Wave object is not well-formed, check for improper attributes")
          
            
    
    def clean(self,tol=None):
        if goodtest(self,test_type='wave'):
            self.kvec = aux_clean(self.kvec,tol)
            self.efield = aux_clean(self.efield,tol)
            self.medium = aux_clean(self.medium,tol)
        else:
            raise Exception("Wave object is not well-formed, check for improper attributes")
            
            
            
    def k(self,k0=None):
        if goodtest(self,test_type='wave'):
            if type(self.kvec) is numpy.ndarray:
                vec = copy.deepcopy(self.kvec)
                vecC = numpy.conj(vec)
                norm = numpy.sqrt(vecC.T @ vec)[0,0]
                if isinstance(k0,(int,float)) and k0>0:
                    return k0*norm.real
                else:
                    return norm.real
            else:
                raise Exception("No kvec found")
        else:
            raise Exception("Wave object is not well-formed, check for improper attributes")
            
            
            
    def propagate(self,vec,k0=None,wrap=None):
        if k0 is None:
            k0 = 1
        elif not isinstance(k0,(int,float)) or k0<=0:
            k0 = 1
            print("Invalid input for k0, setting equal to 1")
        if goodtest(self,test_type='wave'):
            if type(vec) is numpy.ndarray and numpy.shape(vec)==(3,1):
                currentphase = copy.deepcopy(self.phase)
                kvec = copy.deepcopy(self.kvec)
                newphase = currentphase + k0*(kvec.T @ vec)[0,0]
                if wrap:
                    self.phase = newphase % (2*numpy.pi)
                else:
                    self.phase = newphase
            else:
                raise Exception("Propagation direction must be specified as a (3,1) numpy array")
        else:
            raise Exception("Wave object is not well-formed, check for improper attributes")
            
            
            
    def jones(self):
        if goodtest(self,test_type='wave'):
            if type(self.efield) is numpy.ndarray:
                if isinstance(self.phase,(int,float,complex)):
                    return numpy.exp(1j*self.phase)*self.efield
                else:
                    return self.efield
            else:
                raise Exception("No electric field found")
        else:
            raise Exception("Wave object is not well-formed, check for improper attributes")

    
        
    
    # Overloaded methods            
    #-----------------------------------------------------------------------------------------
    def __add__(self,other):
        if isinstance(other,wavecalc.classes.medium):
            new = copy.deepcopy(self)
            new.medium = other.epsilon
            return new
        elif isinstance(other,wavecalc.classes.surface):
            return reflect(self,other,coat=other.coat)
        else:
            raise Exception("Waves can only add with media or surfaces")
     
        
        
    def __sub__(self,other):
        if isinstance(other,wavecalc.classes.surface):
            return transmit(self,other)
        else:
            raise Exception("Waves can only subtract surfaces")
            
            
            
    def __pos__(self):
        new = copy.deepcopy(self)
        new.fixmode()
        new.clean()
        return new
            
            
            
    def __neg__(self):
        if (isinstance(self.efield,numpy.ndarray) 
            and numpy.shape(self.efield)==(3,1)):
            new = copy.deepcopy(self)
            new.efield = -self.efield+0
            return new
        else:
            raise Exception("Wave must have a properly formed efield to be negated")
            
            
            
    def __invert__(self):
        if (isinstance(self.kvec,numpy.ndarray) 
            and numpy.shape(self.kvec)==(3,1)):
            new = copy.deepcopy(self)
            new.kvec = -self.kvec+0
            return new
        else:
            raise Exception("Wave must have properly formed kvec to invert")
            
    
    
    def __matmul__(self,other):
        if isinstance(other,wavecalc.classes.surface):
            return crash(self,other,coat=other.coat,combine_same=True)
        else:
            raise Exception("Waves can only crash onto surfaces")
            
    
    
    def __eq__(self,other):
        if isinstance(other,wavecalc.classes.wave):
            ks = numpy.prod(self.kvec == other.kvec)
            es = numpy.prod(self.efield == other.efield)
            ms = numpy.prod(self.medium == other.medium)
            ph = self.phase == other.phase
            if ks*es*ms*ph:
                return True
            else:
                return False
        else:
            return False

        
        
    
########################################################################################################################
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
########################################################################################################################
class surface:
    ''' A class for planar interfaces '''
    def __init__(self,normal=None,**kwargs):
        
        into = kwargs.pop('into',None)
        out = kwargs.pop('out',None)
        coat = kwargs.pop('coat',None)
        everything = kwargs.pop('everything',None)
        
        if len(kwargs) != 0:
            print("Undefined arguments passed to surface.__init__")
        
        
        
        # normal attribute
        #-------------------------------------------------------------------------------------
        if normal is False or everything is False:
            self.normal = None
        
        elif normal is None:
            self.normal = numpy.array([[0.,0.,1.]]).T
        
        elif type(normal) is numpy.ndarray and numpy.shape(normal) == (3,1):
            self.normal = normal
            
        elif normal == 'random':
            self.normal = random.rand(3,1)
        
        else:
            raise Exception('Must specify normal as (3,1) numpy.ndarray')
            
        
        
        # into attribute
        #-------------------------------------------------------------------------------------
        if into is False or everything is False:
            self.into = None
        
        elif into is None:
            self.into = numpy.array([[1.,0.,0.],
                                     [0.,1.,0],
                                     [0.,0.,1.]])
        
        elif type(into) is numpy.ndarray and numpy.shape(into) == (3,3):
            self.into = into
            
        elif into == 'random':
            self.into = random.rand(3,3)
     
        else:
            raise Exception('Must specify into as a (3,3) numpy.ndarray')
            
        
        
        # out attribute
        #-------------------------------------------------------------------------------------
        if out is False or everything is False:
            self.out = None
        
        elif out is None:
            self.out = numpy.array([[1.,0.,0.],
                                    [0.,1.,0.],
                                    [0.,0.,1.]])
        
        elif type(out) is numpy.ndarray and numpy.shape(out) == (3,3):
            self.out = out
        
        elif out == 'random':
            self.out = random.rand(3,3)
    
        
        else:
            raise Exception('Must specify out as a (3,3) numpy.ndarray')
        
        
        # coat attribute
        #-------------------------------------------------------------------------------------
        if coat is False or everything is False:
            self.coat = None
        
        if coat is not None:
            if (coat == 'HR' or coat == 'hr'):
                self.coat = 'HR'
            elif (coat == 'AR' or coat == 'ar'):
                self.coat = 'AR'
            else:
                raise Exception("Must specify coat as 'HR' or 'AR'")
        else:
            self.coat = None
            
    
    # Custom methods            
    #-----------------------------------------------------------------------------------------
    def rotate(self,ang,axis,**kwargs):
        
        medmove = kwargs.pop('medmove',None)
        fix = kwargs.pop('fix',None)
        verbose = kwargs.pop('verbose',None)
        
        if len(kwargs) != 0:
            print ("Undefined argument passed to surface.rotate")
        
        if fix is not None:
            print("The fix option has no meaning when rotating surfaces")
        rotate_copy(self,ang,axis,medmove,verbose)
        
        
    def clean(self,tol=None):
        if goodtest(self,test_type='surface'):
            self.normal = aux_clean(self.normal,tol)
            self.out = aux_clean(self.out,tol)
            self.into = aux_clean(self.into,tol)
        else:
            raise Exception("Surface object is not well-formed, check for improper attributes")



    # Overloaded methods            
    #-----------------------------------------------------------------------------------------        
    def __add__(self,other):
        if isinstance(other,wavecalc.classes.medium):
            new = copy.deepcopy(self)
            new.into = other.epsilon
            return new
        elif isinstance(other,wavecalc.classes.wave):
            return transmit(other,self,coat=self.coat)
        else:
            raise Exception('Surfaces can only add media or waves')
       
        
        
    def __sub__(self,other):
        if isinstance(other,wavecalc.classes.medium):
            new = copy.deepcopy(self)
            new.out = other.epsilon
            return new
        elif isinstance(other,wavecalc.classes.wave):
            return reflect(other,self,coat=self.coat)
        else:
            raise Exception('Surfaces can only subtract waves or media')



    def __pos__(self):
        new = copy.deepcopy(self)
        new.clean()
        return new
            
    
    
    def __neg__(self):
        if (isinstance(self.normal,numpy.ndarray) 
            and numpy.shape(self.normal) == (3,1)):
            new = copy.deepcopy(self)
            new.normal = -self.normal+0
            return new
        else:
            raise Exception("Surface must have a properly defined normal")
            
            
    
    def __matmul__(self,other):
        if isinstance(other,wavecalc.classes.wave):
            return crash(other,self,coat=self.coat,combine_same=True)
        else:
            raise Exception("Only waves can crash onto surfaces")
    
    
    
    def __invert__(self):
        new = copy.deepcopy(self)
        new.out = self.into
        new.into = self.out
        return new
    
    
    
    def __eq__(self,other):
        if isinstance(other,wavecalc.classes.surface):
            ns = numpy.prod(self.normal == other.normal)
            os = numpy.prod(self.out == other.out)
            Is = numpy.prod(self.into == other.into)
            cs = numpy.prod(self.coat == other.coat)
            if ns*os*Is*cs:
                return True
            else:
                return False
        else:
            return False
    
    
    

        
            
        
        
        
        
        
        
        
########################################################################################################################
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
########################################################################################################################        
class medium:
    ''' A class for dielectric media '''
    def __init__(self,epx=None,epy=None,epz=None,**kwargs):
        
        ep_all = kwargs.pop('ep_all',None)
        epsilon = kwargs.pop('epsilon',None)
        everything = kwargs.pop('everything',None)
        
        if len(kwargs) != 0:
            print("Undefined arguments passed to medium.__init__")
        
        # epsilon attribute
        #-------------------------------------------------------------------------------------
        if (type(epx) is complex 
            or type(epy) is complex 
            or type(epz) is complex 
            or type(ep_all) is complex):
            dat = complex
        else:
            dat = float
        
        if epsilon is False or everything is False:
            self.epsilon = None
        
        elif epsilon is None:
            self.epsilon = numpy.array([[1.,0.,0.],
                                        [0.,1.,0.],
                                        [0.,0.,1.]],dtype=dat)
            if isinstance(epx,(int,float,complex)):
                self.epsilon[0,0] = epx
            if isinstance(epy,(int,float,complex)):
                self.epsilon[1,1] = epy
            if isinstance(epz,(int,float,complex)):
                self.epsilon[2,2] = epz
            if isinstance(ep_all,(int,float,complex)):
                self.epsilon[0,0] = ep_all
                self.epsilon[1,1] = ep_all
                self.epsilon[2,2] = ep_all
            if epx == 'random':
                self.epsilon[0,0] = random.rand()
            if epy == 'random':
                self.epsilon[1,1] = random.rand()
            if epz == 'random':
                self.epsilon[2,2] = random.rand()
            if ep_all == 'random':
                self.epsilon[0,0] = random.rand()
                self.epsilon[1,1] = random.rand()
                self.epsilon[2,2] = random.rand()
        elif type(epsilon) is numpy.ndarray and numpy.shape(epsilon) == (3,3):
            self.epsilon = epsilon
        elif epsilon == 'random':
            self.epsilon = random.rand(3,3)
        else:
            raise Exception('Must specify epsilon as a (3,3) numpy.ndarray')
    
    
        
    # Custom methods            
    #-----------------------------------------------------------------------------------------
    def epx(self,epsilon_xx=None):
        if (type(self.epsilon) is not numpy.ndarray 
            or not numpy.shape(self.epsilon)==(3,3)):
            raise Exception("Medium must have proper epsilon attribute")
        if epsilon_xx is None:
            return self.epsilon[0,0]
        elif isinstance(epsilon_xx,(int,float)):
            self.epsilon[0,0] = epsilon_xx
        elif type(epsilon_xx) is complex:
            self.epsilon = numpy.asarray(self.epsilon,dtype=complex)
            self.epsilon[0,0] = epsilon_xx
        else:
            raise Exception("'epsilon_xx' must be an int, float, or complex")



    def epy(self,epsilon_yy=None):
        if (type(self.epsilon) is not numpy.ndarray 
            or not numpy.shape(self.epsilon)==(3,3)):
            raise Exception("Medium must have proper epsilon attribute")
        if epsilon_yy is None:
            return self.epsilon[1,1]
        elif isinstance(epsilon_yy,(int,float)):
            self.epsilon[1,1] = epsilon_yy
        elif isinstance(epsilon_yy,complex):
            self.epsilon = numpy.asarray(self.epsilon,dtype=complex)
            self.epsilon[1,1] = epsilon_yy
        else:
            raise Exception("'epsilon_yy' must be an int, float, or complex")



    def epz(self,epsilon_zz=None):
        if (type(self.epsilon) is not numpy.ndarray 
            or not numpy.shape(self.epsilon)==(3,3)):
            raise Exception("Medium must have proper epsilon attribute")
        if epsilon_zz is None:
            return self.epsilon[2,2]
        elif isinstance(epsilon_zz,(int,float)):
            self.epsilon[2,2] = epsilon_zz
        elif isinstance(epsilon_zz,complex):
            self.epsilon = numpy.asarray(self.epsilon,dtype=complex)
            self.epsilon[2,2] = epsilon_zz
        else:
            raise Exception("'epsilon_zz' must be an int, float, or complex")



    def ep_all(self,epsilon_all=None):
        if (type(self.epsilon) is not numpy.ndarray 
            or not numpy.shape(self.epsilon)==(3,3)):
            raise Exception("Medium must have proper epsilon attribute")
        if epsilon_all is None:
            return [self.epsilon[0,0],self.epsilon[1,1],self.epsilon[2,2]]
        elif isinstance(epsilon_all,(int,float)):
            self.epsilon[0,0] = epsilon_all
            self.epsilon[1,1] = epsilon_all
            self.epsilon[2,2] = epsilon_all
        elif isinstance(epsilon_all,complex):
            self.epsilon = numpy.asarray(self.epsilon,dtype=complex)
            self.epsilon[0,0] = epsilon_all
            self.epsilon[1,1] = epsilon_all
            self.epsilon[2,2] = epsilon_all
        else:
            raise Exception("'epsilon_all' must be an int, float, or complex")



    def rotate(self,ang,axis,**kwargs):
        
        medmove = kwargs.pop('medmove',None)
        fix = kwargs.pop('fix',None)
        verbose = kwargs.pop('verbose',None)
        
        if len(kwargs) != 0:
            print("Undefined argument passed to medium.rotate")
        
        if fix is not None:
            print("The fix option has no meaning when rotating media")
        rotate_copy(self,ang,axis,medmove,verbose)



    def clean(self,tol=None):
        if goodtest(self):
            self.epsilon = aux_clean(self.epsilon,tol)
        else:
            raise Exception("Medium object is not well-formed, check for improper attributes")



    # Overloaded methods            
    #-----------------------------------------------------------------------------------------    
    def __add__(self,other):
        if isinstance(other,wavecalc.classes.medium):
            return surface(out=self.epsilon,into=other.epsilon)
        elif isinstance(other,wavecalc.classes.surface):
            new = copy.deepcopy(other)
            new.out = self.epsilon
            return new
        elif isinstance(other,wavecalc.classes.wave):
            new = copy.deepcopy(other)
            new.medium = self.epsilon
            return new
        else:
            raise Exception('Cannot add a medium to a non-wavecalc objects')
    
    
        
    def __sub__(self,other):
        if isinstance(other,wavecalc.classes.medium):
                return surface(out=other.epsilon,into=self.epsilon)
        elif isinstance(other,wavecalc.classes.surface):
            new = copy.deepcopy(other)
            new.into = self.epsilon
            return new
        else:
            raise Exception('Cannot add a medium to a non-medium')
            
   
    
    def __pos__(self):
        new = copy.deepcopy(self)
        new.clean()
        return new
    
    
         
    def __eq__(self,other):
        if isinstance(other,wavecalc.classes.medium):
            es = numpy.prod(self.epsilon == other.epsilon)
            if es:
                return True
            else:
                return False
        else:
            return False
        
                





########################################################################################################################
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
########################################################################################################################
class bundle(UserList):
    ''' A class for a bundle of waves '''
    def __init__(self,initdata=None):
        if isinstance(initdata,wavecalc.classes.wave):
            initdata = [initdata]
        else:
            truth = isinstance(initdata,list)*all(isinstance(x,wavecalc.classes.wave) for x in initdata)
            if not truth:
                raise Exception("Bundles can only be made from waves and lists of waves")
        super().__init__(initdata)

    
    
    # Overloaded methods            
    #-----------------------------------------------------------------------------------------
    def append(self,item):
        if isinstance(item,wavecalc.classes.wave):
            self.data.append(item)
        elif isinstance(item,list):
            if all(isinstance(x,wavecalc.classes.wave) for x in item):
                for x in item:
                    self.data.append(x)
            else:
                raise Exception("Can only append waves and lists of waves to a bundle")
        else:
            raise Exception("Can only append waves and lists of waves to a bundle")


    
    def __setitem__(self,i,item):
        if isinstance(item,wavecalc.classes.wave):
            self.data[i] = item
        else:
            raise Exception("Bundle components must be waves")



    def insert(self,i,item):
        if isinstance(item,wavecalc.classes.wave):
            self.data.insert(i,item)
        else:
            raise Exception("Bundle comnponents must be waves")













########################################################################################################################
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
########################################################################################################################
class chain(UserList):
    ''' A class for a chain of surfaces '''
    def __init__(self,initdata=None):
        if isinstance(initdata,wavecalc.classes.surface):
            initdata = [initdata]
        else:
            truth = isinstance(initdata,list)*all(isinstance(x,wavecalc.classes.surface) for x in initdata)
            if not truth:
                raise Exception("Chains can only be made from surfaces and lists of surfaces")
        super().__init__(initdata)



    # Overloaded methods            
    #-----------------------------------------------------------------------------------------
    def append(self,item):
        if isinstance(item,wavecalc.classes.surface):
            self.data.append(item)
        elif isinstance(item,list):
            if all(isinstance(x,wavecalc.classes.surface) for x in item):
                for x in item:
                    self.data.append(x)
            else:
                raise Exception("Can only append surfaces and lists of surfaces to a chain")
        else:
            raise Exception("Can only append surfaces and lists of surfaces to a chain")


    
    def __setitem__(self,i,item):
        if isinstance(item,wavecalc.classes.surface):
            self.data[i] = item
        else:
            raise Exception("Bundle components must be surfaces")



    def insert(self,i,item):
        if isinstance(item,wavecalc.classes.surface):
            self.data.insert(i,item)
        else:
            raise Exception("Bundle comnponents must be surfaces")
        













########################################################################################################################
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
########################################################################################################################
class example_class:
    def __init__(self,example_attribute=None,everything=None):
        
        if everything is False:
            self.example_attribute = None
            
        elif example_attribute == 'random':
            self.example_attribute = random.rand()
            
        else:
            self.example_attribute = example_attribute
        

            
        
        