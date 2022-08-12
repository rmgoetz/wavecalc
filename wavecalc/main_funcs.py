# -*- coding: utf-8 -*-
import numpy
import wavecalc
import copy
from wavecalc.aux_funcs import *
"""
Created on Fri Apr 10 18:48:47 2020
author: Ryan Goetz, ryan.m.goetz@gmail.com
last update: May 20, 2022
"""
'''
Table of Contents:

    Foreground Functions:
        
        crash --------------------------------- Line 33
        
        modes --------------------------------- Line 165
        
        rotate -------------------------------- Line 267
        
        reflect ------------------------------- Line 420
        
        transmit ------------------------------ Line 515


Last line check: April 10, 2020

'''



def crash(wave,surface,**kwargs):
    ''' Outputs the reflected and transmitted waves from 'wave' incident on 'surface' '''
    
     
    ####################################################################################################
    #                                                                                                  #
    # The crash function.                                                                              #
    #                                                                                                  #
    # INPUTS:                                                                                          #
    #         wave - The input wave, given as a wavecalc wave object or (3,1) numpy array.             #
    #      surface - The medium of the transmission, given as a wavecalc medium object or a (3,3)      #
    #                numpy ndarray.                                                                    #
    #           k0 - The magnitude of the wave vector in vacuum, useful for normalizing the results,   #
    #                given as an int or a float.                                                       #
    #         coat - An option to make the surface AR or HR coated with strings 'AR' and 'HR'.         #
    # combine_same - An option that if set to True will combine waves that have very similar wave      #
    #                vectors.                                                                          #                                                                 
    #      verbose - If set to True, prints more information about the calculation.                    #
    #                                                                                                  # 
    #                                                                                                  #
    # Outputs a list of up to four waves resulting from wave incident on surface.                      #
    #                                                                                                  # 
    #                                                                                                  #
    # Last Updated: August 1, 2019                                                                     #
    #                                                                                                  #
    ####################################################################################################
    
    k0 = kwargs.pop('k0',None)
    coat = kwargs.pop('coat',None)
    combine_same = kwargs.pop('combine_same',None)
    verbose = kwargs.pop('verbose',None)
    
    if len(kwargs) != 0:
        print("Undefined arguments passed to wavecalc.functions.crash")
    
    
    # Handle bad surface objects - necessary for the coating handling
    #---------------------------------------------------------------------------------------------------
    if not aux_goodtest(surface,'surface'):
        raise Exception("'surface' argument is not a properly formed wavecalc surface")
    
    
    # Handle the coating choice
    #---------------------------------------------------------------------------------------------------
    if coat is False:
        coat = None
    elif coat not in ['AR','ar','HR','hr']:
        coat = surface.coat
        
    
    # Handle the combine_same option
    #---------------------------------------------------------------------------------------------------
    if not (combine_same is True or combine_same is False):
        combine_same = False
    
    
    # Handle (3,1) arrays as wave input
    #---------------------------------------------------------------------------------------------------
    if isinstance(wave,numpy.ndarray):
        if not numpy.shape(wave)==(3,1):
            raise Exception("If 'wave' parameter is given as a numpy.ndarray it must have shape (3,1)")
        if verbose:
            print("Input parameter 'wave' is assumed to be the wave vector")
        if not isinstance(surface,wavecalc.classes.surface):
            raise Exception("'surface' parameter must be given as a wavecalc surface object")
        if (not isinstance(surface.out,numpy.ndarray) 
            or not numpy.shape(surface.out)==(3,3)):
            raise Exception("No out medium found: surface.out must be a (3,3) numpy ndarray")
        if (not isinstance(surface.into,numpy.ndarray) 
            or not numpy.shape(surface.into)==(3,3)):
            raise Exception("No into medium found: surface.into must be a (3,3) numpy ndarray")
        
        k = wave
        ef = None
        ph = 0.
        s = surface.normal
        ep1 = surface.out
        ep2 = surface.into
        if k0 is None:
            k0 = 1
            if verbose:
                print("Assuming k0 = 1")
        
        sol = aux_waveinterf(k,ef,ph,s,ep1,ep2,k0,coating=coat,verbose=verbose)
        
        return sol
        
                
    # Handle wavecalc waves as wave input
    #---------------------------------------------------------------------------------------------------
    if isinstance(wave,wavecalc.classes.wave):
        if not numpy.shape(wave.kvec)==(3,1):
            raise Exception("If 'wave' parameter is given as a wavecalc wave, kvec attribute must be a numpy.ndarray with shape (3,1)")
        if not isinstance(surface,wavecalc.classes.surface):
            raise Exception("'surface' parameter must be given as a wavecalc surface object")
        if (not isinstance(surface.out,numpy.ndarray) 
            or not numpy.shape(surface.out)==(3,3)):
            raise Exception("No out medium found: surface.out must be a (3,3) numpy ndarray")
        if (not isinstance(surface.into,numpy.ndarray) 
            or not numpy.shape(surface.into)==(3,3)):
            raise Exception("No into medium found: surface.into must be a (3,3) numpy ndarray")
        
        k = wave.kvec
        ef = wave.efield
        ph = wave.phase
        s = surface.normal
        ep1 = surface.out
        ep2 = surface.into
        if k0 is None:
            k0 = 1
            if verbose:
                print("Assuming k0 = 1")
        
        sol = aux_waveinterf(k,ef,ph,s,ep1,ep2,k0,coating=coat,same=combine_same,verbose=verbose)
        
        return sol
    
    
    # Handle unsupported objects
    #---------------------------------------------------------------------------------------------------            
    else:
        raise Exception("'wave' parameter must be given as a wavecalc wave or (3,1) numpy.ndarray")
#
#
#
#
#
#         
#
#
#
#
def modes(ob,med,**kwargs):
    ''' Outputs allowed modes of a medium parallel to input object '''
    
     
    ####################################################################################################
    #                                                                                                  #
    # The modes function.                                                                              #
    #                                                                                                  #
    # INPUTS:                                                                                          #
    #      ob - The input wave, given as a wavecalc wave object or (3,1) numpy array.                  #
    #     med - The medium of the transmission, given as a wavecalc medium object or a (3,3) numpy     #
    #           array.                                                                                 #
    #      k0 - The magnitude of the wave vector in vacuum, useful for normalizing the results, given  #
    #           as an int or a float.                                                                  #                                                                 
    # verbose - If set to True, prints more information about the calculation.                         #
    #                                                                                                  # 
    #                                                                                                  #
    # Outputs a list of wave modes in the medium parallel to ob.                                       #
    #                                                                                                  # 
    #                                                                                                  #
    # Last Updated: April 29, 2020                                                                     #
    #                                                                                                  #
    ####################################################################################################


    k0 = kwargs.pop('k0',None)
    verbose = kwargs.pop('verbose',None)
    
    if len(kwargs) != 0:
        print("Undefined arguments passed to wavecalc.functions.modes")

    # Handle (3,1) arrays for ob input
    #---------------------------------------------------------------------------------------------------
    if isinstance(ob,numpy.ndarray):
        if not numpy.shape(ob)==(3,1):
            raise Exception("If ob is given as a numpy.ndarray it must have shape (3,1)")
        elif isinstance(med,wavecalc.classes.medium):
            if (not numpy.shape(med.epsilon)==(3,3) 
                or not isinstance(med.epsilon,numpy.ndarray)):
                raise Exception("med.epsilon is not properly formed")
            else:
                vec = ob
                medi = med.epsilon
        elif isinstance(med,numpy.ndarray):
            if not numpy.shape(med)==(3,3):
                raise Exception("If medium is given as a numpy.ndarray it must have shape (3,3)")
            else:
                vec = ob
                medi = med
        else:
            raise Exception("med must be given as a wavecalc medium or (3,3) numpy.ndarray")
                
            
    # Handle wavecalc waves for ob input
    #---------------------------------------------------------------------------------------------------
    elif isinstance(ob,wavecalc.classes.wave):
        if not numpy.shape(ob.kvec)==(3,1):
            raise Exception("If ob is given as a wavecalc wave, ob.kvec must be a numpy.ndarray with shape (3,1")
        elif isinstance(med,wavecalc.classes.medium):
            if (not numpy.shape(med.epsilon)==(3,3) 
                or not isinstance(med.epsilon,numpy.ndarray)):
                raise Exception("med.epsilon is not properly formed")
            else:
                vec = ob.kvec
                medi = med.epsilon
        elif isinstance(med,numpy.ndarray):
            if not numpy.shape(med)==(3,3):
                raise Exception("If medium is given as a numpy.ndarray it must have shape (3,3)")
            else:
                vec = ob.kvec
                medi = med        
        else:
            raise Exception("med must be given as a wavecalc medium or (3,3) numpy.ndarray")
    
    
    # Handle unsupported ob types
    #---------------------------------------------------------------------------------------------------            
    else:
        raise Exception("ob must be given as a wavecalc wave or (3,1) numpy.ndarray")
    
    
    # Verify the wave vector is real
    #---------------------------------------------------------------------------------------------------
    realvec = vec.real
    if realvec is not vec:
        raise Exception("Wave vector must be real")
    
    
    # Get the solution
    #---------------------------------------------------------------------------------------------------    
    else:
        return aux_modecalc(vec,medi,k0,verbose)
#
#
#
#
#
#         
#
#
#
#
def rotate(ob,ang,axis,**kwargs):
    ''' Rotates wavecalc objects around a specified axis: 'x', 'y', or 'z' '''
    
    
    ####################################################################################################
    #                                                                                                  #
    # The rotation function for wavecalc object classes.                                               #
    #                                                                                                  #
    # INPUTS:                                                                                          #
    #      ob - The input object to be rotated, either a wavecalc object, (3,1) numpy ndarray or (3,3) #
    #           numpy ndarray.                                                                         #
    #     ang - The angle of rotation in degrees, given as an int or a float.                          #
    #    axis - The axis about which the rotation is to be performed, either 'x', 'y', or 'z'.         #
    # medmove - If set to 'with', the media associated with waves and surfaces rotates as well. If set #
    #           to 'only', the medium (media) is (are) rotated but not the other object attributes.    #
    #           For surfaces: if set to 'into' only the into medium is rotated with the surface, and   #
    #           similarly for 'out'; if set to 'onlyinto' or 'onlyout', these media are rotated while  #
    #           the surface normal remains fixed.                                                      #
    # verbose - If set to True, prints some information about the calculation.                         #
    #                                                                                                  #
    #                                                                                                  #
    # Returns the rotated object in accordance with its transformation properties.                     #
    #                                                                                                  # 
    #                                                                                                  #
    # Last Updated: May 20, 2022                                                                     #
    #                                                                                                  #
    ####################################################################################################
    
    medmove = kwargs.pop('medmove',None)
    verbose = kwargs.pop('verbose',None)
    
    if len(kwargs) != 0:
        print("Undefined arguments passed to wavecalc.functions.rotate")
    
    
    # Reset 'medmove' to None for invalid 'medmove' inputs
    #---------------------------------------------------------------------------------------------------
    medmove_opts = set([None,'with','only','onlyinto','onlyout','into','out'])
    if medmove not in medmove_opts:
        str1 = "If specified, 'medmove' must be set to one of the following: \n"
        str2 = "'with', 'only', 'into', 'out', 'onlyinto', 'onlyout' \n "
        str3 = " Variable 'medmove' will be set to None for following calculation. "
        print(str1+str2+str3)
        medmove = None
        
        
    # Handle (3,1) numpy ndarrays as ob input
    #---------------------------------------------------------------------------------------------------
    if isinstance(ob,numpy.ndarray) and numpy.shape(ob) == (3,1):
        if medmove is not None:
            str1 = "When 'ob' is a numpy ndarray, 'medmove' option has no meaning \n"
            str2 = "and will thus be ignored"
            print(str1+str2)
        return aux_rotvec(ob,ang,axis)
        
    
    # Handle (3,3) numpy ndarrays as ob input
    #---------------------------------------------------------------------------------------------------
    if isinstance(ob,numpy.ndarray) and numpy.shape(ob) == (3,3):
        if medmove is not None:
            str1 = "When 'ob' is a numpy ndarray, 'medmove' option has no meaning \n"
            str2 = "and will thus be ignored"
            print(str1+str2)
        return aux_rottens(ob,ang,axis)
        
        
    # Handle wavecalc wave instances as ob input
    #---------------------------------------------------------------------------------------------------
    if isinstance(ob,wavecalc.classes.wave):
        if not aux_goodtest(ob):
            raise Exception('Your wavecalc wave has improper attributes')
        if medmove in medmove_opts-set([None,'with','only']):
            str1 ="For wavecalc waves, if specified, 'medmove' must be set to one of the following: 'with', 'only'.\n"
            str2 = "Variable 'medmove' will be set to None for following calculation."
            print(str1+str2)
            medmove = None
        if medmove == 'only':
            kk = ob.kvec
            ee = ob.efield
            mm = aux_rottens(ob.medium,ang,axis)
        else:
            kk = aux_rotvec(ob.kvec,ang,axis)
            ee = aux_rotvec(ob.efield,ang,axis)
            mm = ob.medium
            if medmove == 'with':
                mm = aux_rottens(ob.medium,ang,axis)
        if verbose:
            print('New kvec :',kk)
            print('New efield :',ee)
            print('New medium :',mm)
        return wavecalc.classes.wave(kvec=kk,efield=ee,medium=mm)
    
    
    # Handle wavecalc surface instances as ob input
    #---------------------------------------------------------------------------------------------------
    elif isinstance(ob,wavecalc.classes.surface):
        if not aux_goodtest(ob):
            raise Exception('Your wavecalc surface has improper attributes')
        nn = ob.normal
        oo = ob.out
        ii = ob.into
        if medmove == 'only':
            oo = aux_rottens(ob.out,ang,axis)
            ii = aux_rottens(ob.into,ang,axis)
        elif medmove == 'onlyout':
            oo = aux_rottens(ob.out,ang,axis)
        elif medmove == 'onlyinto':
            ii = aux_rottens(ob.into,ang,axis)
        else:
            nn = aux_rotvec(ob.normal,ang,axis)
            if medmove == 'with':
                oo = aux_rottens(ob.out,ang,axis)
                ii = aux_rottens(ob.into,ang,axis)
            elif medmove == 'out':
                oo = aux_rottens(ob.out,ang,axis)
            elif medmove == 'into':
                ii = aux_rottens(ob.into,ang,axis)
        if verbose:
            print('New normal :',nn)
            print('New out :',ii)
            print('New into :',oo)
        return wavecalc.classes.surface(normal=nn,into=ii,out=oo)
    
    
    # Handle wavecalc medium instances as ob input
    #---------------------------------------------------------------------------------------------------        
    elif isinstance(ob,wavecalc.classes.medium):
        if not aux_goodtest(ob):
            raise Exception('Your wavecalc medium has improper attributes')
        if medmove is not None:
            str1 ="For wavecalc media, 'medmove' option has no meaning. \n"
            str2 = "Variable 'medmove' will be set to None for following calculation."
            print(str1+str2)
        ee = aux_rottens(ob.epsilon,ang,axis)
        if verbose:
            print('epsilon :',ee)
        return wavecalc.classes.medium(epsilon=ee)
    
    
    # Handle unsupported object inputs
    #---------------------------------------------------------------------------------------------------
    else:
        raise Exception("Argument 'ob' must be a (3,1) numpy array, (3,3) numpy array, or a wavecalc wave, surface, or medium")
#
#
#
#
#
#         
#
#
#
#
def reflect(wav,surf,**kwargs):
    ''' Outputs reflection waves '''
    
     
    ####################################################################################################
    #                                                                                                  #
    # The wave reflection function.                                                                    #
    #                                                                                                  #
    # INPUTS:                                                                                          #
    #          wav - The input wave, given as a wavecalc wave object.                                  #
    #         surf - The surface normal, given as a wavecalc surface object.                           #
    #           k0 - The magnitude of the wave vector in vacuum, useful for normalizing the results,   #
    #                given as an int or a float.                                                       #
    #         coat - An option to make the surface AR or HR coated with strings 'AR' and 'HR'.         #    
    # combine_same - An option that if set to True will combine waves that have very similar wave      #
    #                vectors.                                                                          #                                                                                  
    #      verbose - If set to True, prints more information about the calculation.                    #
    #                                                                                                  # 
    #                                                                                                  #
    # Outputs a list of up to two reflection wavecalc wave objects.                                    #
    #                                                                                                  # 
    #                                                                                                  #
    # Last Updated: July 31, 2019                                                                      #
    #                                                                                                  #
    ####################################################################################################
    
    
    k0 = kwargs.pop('k0',None)
    coat = kwargs.pop('coat',None)
    combine_same = kwargs.pop('combine_same',None)
    verbose = kwargs.pop('verbose',None)
    
    if len(kwargs) != 0:
        print("Undefined arguments passed to wavecalc.functions.reflect")
    
    
    # Handle the coating choice
    #---------------------------------------------------------------------------------------------------
    coats=['AR','ar','HR','hr']
    if coat is not None and coat not in coats:
        coat = None
        
    if not (combine_same is True or combine_same is False):
        combine_same = False
    
    
    # Handle k0
    #---------------------------------------------------------------------------------------------------
    if k0 is None:
        k0 = 1
    elif not isinstance(k0,(int,float)):
        raise Exception("If specified, 'k0' must be an int or float")
    
    
    # Handle wav
    #---------------------------------------------------------------------------------------------------
    if (isinstance(wav,wavecalc.classes.wave) 
            and aux_goodtest(wav) 
            and isinstance(wav.kvec,numpy.ndarray)):
        k = wav.kvec
        ph = wav.phase
        ef = wav.efield

    else:
        raise Exception("The 'wav' input must be a wavecalc wave with proper attributes")
        
    
    # Handle surf
    #---------------------------------------------------------------------------------------------------
    if (isinstance(surf,wavecalc.classes.surface) 
                and aux_goodtest(surf) 
                and isinstance(surf.normal,numpy.ndarray)
                and isinstance(surf.out,numpy.ndarray)
                and isinstance(surf.into,numpy.ndarray)):
        s = surf.normal
        ep1 = surf.out
        ep2 = surf.into
    else:
        raise Exception("The 'surf' input must be a wavecalc surface with proper normal, out, and into attributes")
    
    
    # Get the output
    #---------------------------------------------------------------------------------------------------
    sol = aux_waveinterf(k,ef,ph,s,ep1,ep2,k0,act='refl',coating=coat,same=combine_same,verbose=verbose)     
    return sol
#
#
#
#
#
#         
#
#
#
#
def transmit(wav,surf,**kwargs):
    ''' Outputs transmission waves '''
    
     
    ####################################################################################################
    #                                                                                                  #
    # The wave transmission function.                                                                  #
    #                                                                                                  #
    # INPUTS:                                                                                          #
    #          wav - The input wave, given as a wavecalc wave object.                                  #
    #         surf - The surface normal, given as a wavecalc surface object.                           #
    #           k0 - The magnitude of the wave vector in vacuum, useful for normalizing the results,   #
    #                given as an int or a float.                                                       #
    #         coat - An option to make the surface AR or HR coated with strings 'AR' and 'HR'.         #
    # combine_same - An option that if set to True will combine waves that have very similar wave      #
    #                vectors.                                                                          #                                                                                                                                                   
    #      verbose - If set to True, prints more information about the calculation.                    #
    #                                                                                                  # 
    #                                                                                                  #
    # Outputs a list of up to two transmitted wavecalc wave objects.                                   #
    #                                                                                                  # 
    #                                                                                                  #
    # Last Updated: July 31, 2019                                                                      #
    #                                                                                                  #
    ####################################################################################################
    
    k0 = kwargs.pop('k0',None)
    coat = kwargs.pop('coat',None)
    combine_same = kwargs.pop('combine_same',None)
    verbose = kwargs.pop('verbose',None)
    
    if len(kwargs) != 0:
        print("Undefined arguments passed to wavecalc.functions.transmit")
    
    
    # Handle the coating choice
    #---------------------------------------------------------------------------------------------------
    coats=['AR','ar','HR','hr']    
    if coat is not None and coat not in coats:
        coat = None
    
    
    if not (combine_same is True or combine_same is False):
        combine_same = False
    
    
    # Handle k0
    #---------------------------------------------------------------------------------------------------
    if k0 is None:
        k0 = 1
    elif not isinstance(k0,(int,float)):
        raise Exception("If specified, 'k0' must be an int or float")
    
    
    # Handle wav
    #---------------------------------------------------------------------------------------------------
    if (isinstance(wav,wavecalc.classes.wave) 
            and aux_goodtest(wav,type_test='wave') 
            and isinstance(wav.kvec,numpy.ndarray)):
        k = wav.kvec
        ph = wav.phase
        ef = wav.efield

    else:
        raise Exception("The 'wav' input must be a wavecalc wave with proper attributes")
        
    
    # Handle surf
    #---------------------------------------------------------------------------------------------------
    if (isinstance(surf,wavecalc.classes.surface) 
                and aux_goodtest(surf,type_test='surface') 
                and isinstance(surf.normal,numpy.ndarray)
                and isinstance(surf.out,numpy.ndarray)
                and isinstance(surf.into,numpy.ndarray)):
        s = surf.normal
        ep1 = surf.out
        ep2 = surf.into
    else:
        raise Exception("The 'surf' input must be a wavecalc surface attributes")
    
    
    # Get the output
    #---------------------------------------------------------------------------------------------------
    sol = aux_waveinterf(k,ef,ph,s,ep1,ep2,k0,act='trans',coating=coat,same=combine_same,verbose=verbose)        
    return sol