# -*- coding: utf-8 -*-
import numpy
import wavecalc
'''
Created May 20, 2019
author: Ryan Goetz, ryan.m.goetz@gmail.com
last update: June 9, 2019 21:39 EST
'''
'''
Table of Contents:

    Foreground Functions:
        
        crash --------------------------------- Line 100
        
        modes --------------------------------- Line 240
        
        rotate -------------------------------- Line 354
        
        reflect ------------------------------- Line 519
        
        transmit ------------------------------ Line 624
        
    
    Background Functions:
        
        aux_booker_interf --------------------- Line 730
        
        aux_check_same ------------------------ Line 883
        
        aux_clean ----------------------------- Line 1004
        
        aux_coat_handle ----------------------- Line 1109
        
        aux_coord_transform ------------------- Line 1205
        
        aux_field_match ----------------------- Line 1337
        
        aux_fixmode --------------------------- Line 1446
        
        aux_goodtest -------------------------- Line 1585
        
        aux_goodtest_wav ---------------------- Line 1678
        
        aux_goodtest_surf --------------------- Line 1746
        
        aux_goodtest_med ---------------------- Line 1820
        
        aux_maxwell_eigenvec ------------------ Line 1874
        
        aux_modecalc -------------------------- Line 1979
        
        aux_quarttest ------------------------- Line 2086
        
        aux_realtest -------------------------- Line 2146
        
        aux_rotate_copy ----------------------- Line 2260
        
        aux_rotmatrix ------------------------- Line 2392
        
        aux_rottens --------------------------- Line 2463

        aux_rotvec ---------------------------- Line 2515
        
        aux_waveinterf ------------------------ Line 2567
           
        
    Under Construction:
        
        aux_complex_killer -------------------- Line 2852
        
        aux_ferrari --------------------------- Line 2909
        
        aux_root_order ------------------------ Line 2958


Last line check: June 9, 2019

'''




















def crash(wave,surface,k0=None,coat=None,combine_same=None,verbose=None):
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
    # Last Updated: June 2, 2019                                                                       #
    #                                                                                                  #
    ####################################################################################################
    
    
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
        if verbose == True:
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
        s = surface.normal
        ep1 = surface.out
        ep2 = surface.into
        if k0 is None:
            k0 = 1
            if verbose == True:
                print("Assuming k0 = 1")
        
        sol = aux_waveinterf(k,ef,s,ep1,ep2,k0,coating=coat,verbose=verbose)
        
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
        s = surface.normal
        ep1 = surface.out
        ep2 = surface.into
        if k0 is None:
            k0 = 1
            if verbose == True:
                print("Assuming k0 = 1")
        
        sol = aux_waveinterf(k,ef,s,ep1,ep2,k0,coating=coat,same=combine_same,verbose=verbose)
        
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
#
#
#
#
#
#
#         
#
def modes(ob,med,k0=None,verbose=None):
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
    # Last Updated: May 25, 2019                                                                       #
    #                                                                                                  #
    ####################################################################################################


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
    if isinstance(ob,wavecalc.classes.wave):
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
    if not realvec is vec:
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
#
#
#
#
#
#
#         
#
def rotate(ob,ang,axis,medmove=None,verbose=None):
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
    # Last Updated: May 25, 2019                                                                       #
    #                                                                                                  #
    ####################################################################################################
    
    
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
        if not medmove is None:
            str1 = "When 'ob' is a numpy ndarray, 'medmove' option has no meaning \n"
            str2 = "and will thus be ignored"
            print(str1+str2)
        return aux_rotvec(ob,ang,axis)
        
    
    # Handle (3,3) numpy ndarrays as ob input
    #---------------------------------------------------------------------------------------------------
    if isinstance(ob,numpy.ndarray) and numpy.shape(ob) == (3,3):
        if not medmove is None:
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
        if verbose == True:
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
        if verbose == True:
            print('New normal :',nn)
            print('New out :',ii)
            print('New int :',oo)
        return wavecalc.classes.surface(normal=nn,into=ii,out=oo)
    
    
    # Handle wavecalc medium instances as ob input
    #---------------------------------------------------------------------------------------------------        
    elif isinstance(ob,wavecalc.classes.medium):
        if not aux_goodtest(ob):
            raise Exception('Your wavecalc medium has improper attributes')
        if not medmove is None:
            str1 ="For wavecalc media, 'medmove' option has no meaning. \n"
            str2 = "Variable 'medmove' will be set to None for following calculation."
            print(str1+str2)
        ee = aux_rottens(ob.epsilon,ang,axis)
        if verbose == True:
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
#
#
#
#
#
#
#         
#
def reflect(wav,surf,k0=None,coat=None,combine_same=None,verbose=None):
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
    # Last Updated: June 2, 2019                                                                       #
    #                                                                                                  #
    ####################################################################################################
    
    
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
        if wav.efield is not None:
            ef = wav.efield
        else:
            ef = None
    else:
        raise Exception("The 'wav' input must be a wavecalc wave with proper kvec attribute")
        
    
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
    sol = aux_waveinterf(k,ef,s,ep1,ep2,k0,act='refl',coating=coat,same=combine_same,verbose=verbose)     
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
#
#
#
#
#
#
#         
#       
def transmit(wav,surf,k0=None,coat=None,combine_same=None,verbose=None):
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
    # Last Updated: June 1, 2019                                                                       #
    #                                                                                                  #
    ####################################################################################################
    
    
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
        if wav.efield is not None:
            ef = wav.efield
        else:
            ef = None
    else:
        raise Exception("The 'wav' input must be a wavecalc wave with proper kvec attribute")
        
    
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
    sol = aux_waveinterf(k,ef,s,ep1,ep2,k0,act='trans',coating=coat,same=combine_same,verbose=verbose)        
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
#
#
#
#
#
#
#         
#
def aux_booker_interf(kx,med1,med2,verbose=None):
    ''' A behind-the-scenes function for solving the Booker equation at a boundary '''
    
   
    ####################################################################################################
    #                                                                                                  #
    # The Booker quartic solution at an interface                                                      #
    #                                                                                                  #
    # INPUTS:                                                                                          #
    #      kx - The the x component of normalized input wave vector in solver coordinates, given as a  #
    #           float or int.                                                                          #
    #    med1 - The incident/reflection medium given as a (3,3) numpy ndarray. Must be given in solver #
    #           coordinates.                                                                           #
    #    med2 - The transmission medium given as a (3,3) numpy ndarray. Must be given in solver        #
    #           coordinates.                                                                           #
    # verbose - If set to True, prints more information about the calculation                          #
    #                                                                                                  # 
    #                                                                                                  #
    # Outputs a list of the two reflection wave vectors and two transmission wave vectors as (3,1)     #
    # numpy ndarrays.                                                                                  # 
    #                                                                                                  # 
    #                                                                                                  #
    # Last Updated: May 30, 2019                                                                       #
    #                                                                                                  #
    ####################################################################################################
    

    
    # Compute the minors of epsilon 1 in the solution coordinates
    #---------------------------------------------------------------------------------------------------
    M1 = numpy.zeros([3,3],dtype=complex)
    for i in range(3):
        for j in range(3):
            minora = numpy.delete(med1,i,axis=0)
            minor = numpy.delete(minora,j,axis=1)
            M1[i,j] = numpy.linalg.det(minor)
            
    
    # Compute the minors of epsilon 2 in the solution coordinates
    #---------------------------------------------------------------------------------------------------
    M2 = numpy.zeros([3,3],dtype=complex)
    for i in range(3):
        for j in range(3):
            minora = numpy.delete(med2,i,axis=0)
            minor = numpy.delete(minora,j,axis=1)
            M2[i,j] = numpy.linalg.det(minor)
   
    
    # Some convenient quantities to calculate
    #---------------------------------------------------------------------------------------------------
    ID = numpy.array([[1,0,0],[0,1,0],[0,0,1]])
    
    detep1 = med1[0,0]*M1[0,0]-med1[0,1]*M1[0,1]+med1[0,2]*M1[0,2]
    sig1 = med1+med1.T
    delt1 = (med1[0,0]+med1[1,1]+med1[2,2])*ID-med1
    
    detep2 = med2[0,0]*M2[0,0]-med2[0,1]*M2[0,1]+med2[0,2]*M2[0,2]
    sig2 = med2+med2.T
    delt2 = (med2[0,0]+med2[1,1]+med2[2,2])*ID-med2
    
    
    # The coefficients of the two (refl and trans) quartics in kz/k0
    #---------------------------------------------------------------------------------------------------
    DELTA1 = kx*sig1[0,2]
    SIGMA1 = (kx**2)*delt1[1,1]-(M1[0,0]+M1[1,1])
    PSI1 = (kx**3)*sig1[0,2]+kx*(M1[0,2]+M1[2,0])
    GAMMA1 = (kx**4)*med1[0,0]-(kx**2)*(M1[1,1]+M1[2,2])+detep1
    
    DELTA2 = kx*sig2[0,2]
    SIGMA2 = (kx**2)*delt2[1,1]-(M2[0,0]+M2[1,1])
    PSI2 = (kx**3)*sig2[0,2]+kx*(M2[0,2]+M2[2,0])
    GAMMA2 = (kx**4)*med2[0,0]-(kx**2)*(M2[1,1]+M2[2,2])+detep2
   
    
    # Verbosity
    #---------------------------------------------------------------------------------------------------
    if verbose==True:
        print('DELTA_r =',DELTA1)
        print("SIGMA_r =",SIGMA1)
        print("PSI_r =",PSI1)
        print("GAMMA_r =",GAMMA1)
        print('DELTA_t =',DELTA2)
        print("SIGMA_t =",SIGMA2)
        print("PSI_t =",PSI2)
        print("GAMMA_t =",GAMMA2)
    
    
    
    # Roots of the quartics
    #---------------------------------------------------------------------------------------------------
    coeffs_low_to_high_1 = [GAMMA1,PSI1,SIGMA1,DELTA1,med1[2,2]]
    coeffs_low_to_high_2 = [GAMMA2,PSI2,SIGMA2,DELTA2,med2[2,2]]
    kzs1 = numpy.polynomial.polynomial.polyroots(coeffs_low_to_high_1)
    kzs2 = numpy.polynomial.polynomial.polyroots(coeffs_low_to_high_2)
    kzs1 = numpy.sort(kzs1)
    kzs2 = numpy.sort(kzs2)
    
    kzout = [kzs1[0],kzs1[1],kzs2[2],kzs2[3]]
    
    
    # Eliminate small imaginary numerical junk
    #---------------------------------------------------------------------------------------------------
    for i in range(4):
        if abs(kzout[i].imag) < 1e-10:
            kzout[i] = kzout[i].real
    
    
    # Verbosity
    #---------------------------------------------------------------------------------------------------
    if verbose==True:
        print("Quartic roots are approximated as: ",kzout)
        
    
    # Verify that the roots have been properly solved for
    #---------------------------------------------------------------------------------------------------
    if not aux_quarttest(coeffs_low_to_high_1,kzs1):
        raise Exception("Quartic root solver failed to sufficiently approximate reflection roots")
        
    if not aux_quarttest(coeffs_low_to_high_2,kzs2):
        raise Exception("Quartic root solver failed to sufficiently approximate transmission roots")
    
    
    # Return the roots
    #---------------------------------------------------------------------------------------------------
    return kzout
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
#
#
#
#
#
#
#         
#    
def aux_check_same(lis,switch=None):
    ''' A behind-the-scenes function for combining waves with the same wave vector '''
    
   
    ####################################################################################################
    #                                                                                                  #
    # The sameness checker                                                                             #
    #                                                                                                  #
    # INPUTS:                                                                                          #
    #    lis - A list of two or four wavecalc waves to be tested for sameness.                         #
    # switch - An option for running the sameness checker. If set to False, the function returns its   #
    #         input.                                                                                   #
    #                                                                                                  # 
    #                                                                                                  #
    # Outputs the input list with like wave vector waves combined. Only compares lis[0] to lis[1], and #
    # lis[2] to lis[3], and vice versa. Does not compare lis[0] to lis[3], for example.                #
    #                                                                                                  # 
    #                                                                                                  #
    # Last Updated: June 2, 2019                                                                       #
    #                                                                                                  #
    ####################################################################################################
    
    
    # Handle te switch options
    #---------------------------------------------------------------------------------------------------
    if switch is False:
        return lis
    
    else:
        ell = len(lis)
        LIS = []
        
        # Handle when lis is two waves
        #-----------------------------------------------------------------------------------------------
        if ell == 2:
            E1 = lis[0]
            E2 = lis[1]
            E1c = E1 
            E2c = E2
            E1c.clean()
            E2c.clean()
            delta = E1c.kvec-E2c.kvec
            diff = (numpy.conj(delta).T @ delta)[0,0].real
            if diff < 1e-14:
                efsum = E1.efield +E2.efield
                k = E1.kvec
                med = E1.medium
                wave = wavecalc.classes.wave(kvec=k,efield=efsum,medium=med)
                LIS.append(wave)
            else:
                LIS.append(E1)
                LIS.append(E2)
            return LIS
        
        
        # Handle when lis is four waves
        #-----------------------------------------------------------------------------------------------
        else:
            E1 = lis[0]
            E2 = lis[1]
            E3 = lis[2]
            E4 = lis[3]
            E1c = E1 
            E2c = E2
            E3c = E3
            E4c = E4
            E1c.clean()
            E2c.clean()
            E3c.clean()
            E4c.clean()
            delta1 = E1c.kvec-E2c.kvec
            delta2 = E3c.kvec-E4c.kvec
            diff1 = (numpy.conj(delta1).T @ delta1)[0,0].real
            diff2 = (numpy.conj(delta2).T @ delta2)[0,0].real
            if diff1 < 1e-14:
                efsum = E1.efield +E2.efield
                k = E1.kvec
                med = E1.medium
                wave = wavecalc.classes.wave(kvec=k,efield=efsum,medium=med)
                LIS.append(wave)
            else:
                LIS.append(E1)
                LIS.append(E2)
            if diff2 < 1e-14:
                efsum = E3.efield +E4.efield
                k = E3.kvec
                med = E3.medium
                wave = wavecalc.classes.wave(kvec=k,efield=efsum,medium=med)
                LIS.append(wave)
            else:
                LIS.append(E3)
                LIS.append(E4)
            return LIS
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
#
#
#
#
#
#
#         
#    
def aux_clean(ob,tolerance=None):
    ''' A behind-the-scenes function for cleaning vectors and matrices of numerical dirt '''
    
   
    ####################################################################################################
    #                                                                                                  #
    # The cleaning function                                                                            #
    #                                                                                                  #
    # INPUTS:                                                                                          #
    #        ob - The object to be cleaned, either a (3,1) or (3,3) numpy ndarray                      #
    # tolerance - The tolerance to consider when cleaning, as a float or int. Defaults to 1e-8.        #
    #                                                                                                  # 
    #                                                                                                  #
    # Outputs the cleaned version of the object, for a given tolerance                                 #
    #                                                                                                  # 
    #                                                                                                  #
    # Last Updated: June 2, 2019                                                                       #
    #                                                                                                  #
    ####################################################################################################
    
    
    # Clone the input object
    #---------------------------------------------------------------------------------------------------
    dummy = ob
    
    
    # Handle the tolerance input
    #---------------------------------------------------------------------------------------------------
    if tolerance is None:
        tol = 1e-8
    elif isinstance(tolerance,(int,float)):
        tol = float(tolerance)
    else:
        raise Exception("tolerance can only be specified as a float or an int")
        
        
    # Determine the digits
    #---------------------------------------------------------------------------------------------------
    dig = -round(numpy.log(tol)/numpy.log(10)) 
    dig = int(dig)
    
    
    # Handle Nonetype input for ob
    #---------------------------------------------------------------------------------------------------
    if ob is None:
        return None
    
    
    # Handle (3,1) numpy ndarray as ob input
    #---------------------------------------------------------------------------------------------------
    elif numpy.shape(ob)==(3,1):
        for i in range(3):
            digg = dig
            remainder = numpy.abs(round(dummy[i,0],digg)-dummy[i,0])
            while remainder < tol and remainder < numpy.abs(dummy[i,0]):
                digg = digg-1
                remainder = numpy.abs(round(dummy[i,0],digg)-dummy[i,0])
            dummy[i,0] = round(dummy[i,0],digg+1)
        if numpy.isreal(dummy).all():
            dummy = numpy.asarray(dummy,dtype=float)
        return dummy
    
    
    # Handle (3,3) numpy ndarray as ob input
    #---------------------------------------------------------------------------------------------------
    else:
        for i in range(3):
            for j in range(3):
                digg = dig
                remainder = numpy.abs(round(dummy[i,j],digg)-dummy[i,j])
                while remainder < tol and remainder < numpy.abs(dummy[i,j]):
                    digg = digg-1
                    remainder = numpy.abs(round(dummy[i,j],digg)-dummy[i,j])
                dummy[i,j] = round(dummy[i,j],digg+1)
        if numpy.isreal(dummy).all():
            dummy = numpy.asarray(dummy,dtype=float)
        return dummy
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
#
#
#
#
#
#
#         
#    
def aux_coat_handle(Ein,Ea,Eb,Eg,En,ep1,ep2,coat):
    ''' A behind-the-scenes function for handling AR and HR surface coatings '''
    
   
    ####################################################################################################
    #                                                                                                  #
    # The coating function                                                                             #
    #                                                                                                  #
    # INPUTS:                                                                                          #
    #  Ein - The incident wave, as a (3,1) numpy ndarray.                                              #
    #   Ea - The alpha wave, as a (3,1) numpy ndarray.                                                 #
    #   Eb - The beta wave, as a (3,1) numpy ndarray.                                                  #
    #   Eg - The gamma wave, as a (3,1) numpy ndarray.                                                 #
    #   En - The nu wave, as a (3,1) numpy ndarray.                                                    #
    #  ep1 - The incident/reflection material dielectric tensor, as a (3,3) numpy ndarray.             #
    #  ep2 - The refraction/transmission material dielectric tensor, as a (3,3) numpy ndarray.         #
    # coat - An option to make the surface AR or HR coated with strings 'AR' and 'HR'.                 #
    #                                                                                                  # 
    #                                                                                                  #
    # Outputs a list of the reflection and transmission fields after adjusting amplitudes to           #
    # correspond with coating properties.                                                              #
    #                                                                                                  # 
    #                                                                                                  #
    # Last Updated: June 2, 2019                                                                       #
    #                                                                                                  #
    ####################################################################################################
    
    
    # The energy density of the incident wave
    #---------------------------------------------------------------------------------------------------
    dfield_inc = ep1 @ Ein
    u_inc = (numpy.conj(Ein).T @ dfield_inc)[0,0].real
    
    
    # Handle the case of no coating
    #---------------------------------------------------------------------------------------------------
    if coat is None:
        return [Ea,Eb,Eg,En]
    
    
    # Handle the AR case
    #---------------------------------------------------------------------------------------------------
    elif coat == 'ar' or coat == 'AR':
        Ea = numpy.array([[0.,0.,0.]]).T
        Eb = numpy.array([[0.,0.,0.]]).T
        dfield_g = ep2 @ Eg
        dfield_n = ep2 @ En
        u_tran_g  = (numpy.conj(Eg).T @ dfield_g)[0,0].real
        u_tran_n  = (numpy.conj(En).T @ dfield_n)[0,0].real
        u_tran = u_tran_g + u_tran_n
        Eg = numpy.sqrt(u_inc/u_tran)*Eg
        En = numpy.sqrt(u_inc/u_tran)*En
        return [Ea,Eb,Eg,En] 
        
        
    # Handle the HR case    
    #---------------------------------------------------------------------------------------------------
    else:
        Eg = numpy.array([[0.,0.,0.]]).T
        En = numpy.array([[0.,0.,0.]]).T
        dfield_a = ep1 @ Ea
        dfield_b = ep1 @ Eb
        u_refl_a  = (numpy.conj(Ea).T @ dfield_a)[0,0].real
        u_refl_b  = (numpy.conj(Eb).T @ dfield_b)[0,0].real
        u_refl = u_refl_a + u_refl_b
        Ea = numpy.sqrt(u_inc/u_refl)*Ea
        Eb = numpy.sqrt(u_inc/u_refl)*Eb
        return [Ea,Eb,Eg,En] 
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
#
#
#
#
#
#
#         
#
def aux_coord_transform(k,s,verbose=None):
    ''' A behind-the-scenes function calculating the solver coordinate transform matrix and inverse '''
    
   
    ####################################################################################################
    #                                                                                                  #
    # The solver coordinate transform calculator                                                       #
    #                                                                                                  #
    # INPUTS:                                                                                          #
    #       k - The input wave vector given as a (3,1) numpy ndarray.                                  #
    #       s - The surface normal vector given as a (3,1) numpy ndarray.                              #
    # verbose - If set to True, prints more information about the calculation.                         #
    #                                                                                                  # 
    #                                                                                                  #
    # Outputs a list of the transformation matrix U and its inverse as (3,3) numpy ndarrays            #
    #                                                                                                  # 
    #                                                                                                  #
    # Last Updated: June 2, 2019                                                                       #
    #                                                                                                  #
    ####################################################################################################   
    
    
    # Normalize the surface vector, set zp to s, and find the magnitude of the wave vector
    #---------------------------------------------------------------------------------------------------
    sstar = numpy.conj(s)
    snorm = numpy.sqrt((sstar.T @ s)[0,0]).real
    S = s/snorm
    zp = S
    kstar = numpy.conj(k)
    knorm = numpy.sqrt((kstar.T @ k)[0,0]).real
    
    
    # Check if the wave is incident on the surface 
    #---------------------------------------------------------------------------------------------------
    kdotS = (k.T @ S)[0,0] 
    if kdotS.real<=0:
        str1 = 'Wave vector is not incident on interface'
        raise Exception(str1)
    
    
    # Naively calculate xp
    #---------------------------------------------------------------------------------------------------
    xpa = k-kdotS*S
    xpastar = numpy.conj(xpa)
    xnorm = numpy.sqrt((xpastar.T @ xpa)[0,0]).real
    
    
    # Handle when k and s are parallel 
    #---------------------------------------------------------------------------------------------------
    if xnorm<1e-14:
        
        xpa = numpy.array([[0.,0.,0.]]).T
        full_index = numpy.array([0,1,2])
        where_zeros = numpy.where(zp<1e-14)[0]
        where_not = numpy.setdiff1d(full_index,where_zeros)
        num_zeros = len(where_zeros)
        
        # Case where only 1 component of zp is non-zero
        #-----------------------------------------------------------------------------------------------
        if num_zeros==2:
            j_ind = where_zeros[0]
            xpa[j_ind,0] = 1  
            xp = xpa
            
        # Case where at least two components of zp are non-zero
        #-----------------------------------------------------------------------------------------------
        else:
            ij = where_not[:2]
            i_ind = ij[0]
            j_ind = ij[1]
            xpa[i_ind,0] = 1
            xpa[j_ind,0] = -zp[i_ind,0]/zp[j_ind,0]
            xpastar = numpy.conj(xpa)
            newnorm = numpy.sqrt(1+(zp[i_ind,0]/zp[j_ind,0])**2)  
            xp = xpa/newnorm
                
    # Else when k and s are not parallel 
    #---------------------------------------------------------------------------------------------------
    else:
        xp=xpa/xnorm
    
    yp = numpy.cross(zp.T,xp.T).T
    
    
    # Determine the transformation matrix and its indices
    #---------------------------------------------------------------------------------------------------    
    U = numpy.array([[xp[0,0],yp[0,0],zp[0,0]],[xp[1,0],yp[1,0],zp[1,0]],[xp[2,0],yp[2,0],zp[2,0]]])
    Uinv = numpy.linalg.inv(U)
    
    
    # Verbosity
    #---------------------------------------------------------------------------------------------------
    if verbose==True:
        print("x' =",xp.T)
        print("y' =",yp.T)
        print("z' =",zp.T)
        print('k_hat . s_hat = ',kdotS/knorm)
        print('U = ',U)
        print('Uinv = ',Uinv)
    

    # Return the transform matrices
    #---------------------------------------------------------------------------------------------------
    return [U,Uinv]
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
#
#
#
#
#
#
#         
#
def aux_field_match(k_alpha,k_beta,k_gamma,k_nu,alpha,beta,gamma,nu,kin,Ein,verbose=None):
    ''' A behind-the-scenes function matching electric fields with boundary conditions '''
    
   
    ####################################################################################################
    #                                                                                                  #
    # The field matching function                                                                      #
    #                                                                                                  #
    # INPUTS:                                                                                          #
    # k_alpha - The alpha wave vector in solver coordinates given as a (3,1) numpy ndarray.            #
    #  k_beta - The beta wave vector in solver coordinates given as a (3,1) numpy ndarray.             #
    # k_gamma - The gamma wave vector in solver coordinates given as a (3,1) numpy ndarray.            #
    #    k_nu - The nu wave vector in solver coordinates given as a (3,1) numpy ndarray.               #
    #   alpha - The normed alpha field vector in solver coordinates given as a (3,1) numpy ndarray.    #
    #    beta - The normed beta field vector in solver coordinates given as a (3,1) numpy ndarray.     #
    #   gamma - The normed gamma field vector in solver coordinates given as a (3,1) numpy ndarray.    #
    #      nu - The normed nu field vector in solver coordinates given as a (3,1) numpy ndarray.       #
    #    k_in - The incident wave vector in solver coordinates given as a (3,1) numpy ndarray.         #
    #     Ein - The incident field vector in solver coordinates given as a (3,1) numpy ndarray.        #
    # verbose - If set to True, prints more information about the calculation                          #
    #                                                                                                  # 
    #                                                                                                  #
    # Outputs the four reflection and transmission coefficients as a (4,) numpy ndarray.               # 
    #                                                                                                  # 
    #                                                                                                  #
    # Last Updated: May 30, 2019                                                                       #
    #                                                                                                  #
    ####################################################################################################
    
    
    # Initialize the system matrix and result vector 
    #---------------------------------------------------------------------------------------------------
    M = numpy.zeros([4,4],dtype=complex)
    b = numpy.zeros([4],dtype=complex)
    
    
    # The first row of the system matrix
    #---------------------------------------------------------------------------------------------------
    M[0,0] = alpha[0,0]
    M[0,1] = beta[0,0]
    M[0,2] = -gamma[0,0]
    M[0,3] = -nu[0,0]
    
    # The second row of the system matrix
    #---------------------------------------------------------------------------------------------------
    M[1,0] = alpha[1,0]
    M[1,1] = beta[1,0]
    M[1,2] = -gamma[1,0]
    M[1,3] = -nu[1,0]
    
    # The third row of the system matrix
    #---------------------------------------------------------------------------------------------------
    M[2,0] = k_alpha[2,0]*alpha[1,0]
    M[2,1] = k_beta[2,0]*beta[1,0]
    M[2,2] = -k_gamma[2,0]*gamma[1,0]
    M[2,3] = -k_nu[2,0]*nu[1,0]
    
    # The fourth row of the system matrix
    #---------------------------------------------------------------------------------------------------
    M[3,0] = k_alpha[2,0]*alpha[0,0] - k_alpha[0,0]*alpha[2,0]
    M[3,1] = k_beta[2,0]*beta[0,0] - k_beta[0,0]*beta[2,0]
    M[3,2] = k_gamma[0,0]*gamma[2,0] - k_gamma[2,0]*gamma[0,0]
    M[3,3] = k_nu[0,0]*nu[2,0] - k_nu[2,0]*nu[0,0]
    
    
    # The result vector
    #---------------------------------------------------------------------------------------------------
    b[0] = -Ein[0,0]
    b[1] = -Ein[1,0]
    b[2] = -kin[2,0]*Ein[1,0]
    b[3] = kin[0,0]*Ein[2,0]-kin[2,0]*Ein[0,0]
    
    
    # Get the solution vector
    #---------------------------------------------------------------------------------------------------
    sol = numpy.linalg.solve(M,b)
    
    
    # Return the solution vector
    #---------------------------------------------------------------------------------------------------
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
#
#
#
#
#
#
#         
#    
def aux_fixmode(wave=None,ab=None,k0=None,verbose=None):
    ''' A behind-the-scenes function for rectifying waves with their media '''
    
   
    ####################################################################################################
    #                                                                                                  #
    # The fixmode function                                                                             #
    #                                                                                                  #
    # INPUTS:                                                                                          #
    #   kvec - The wave vector of the wave to be fixed, as a (3,1) numpy ndarray.                      #
    # efield - The wave vector of the wave to be fixed, as a (3,1) numpy ndarray.                      #
    #    med - The medium of the wave as a (3,3) numpy ndarray.                                        #
    #     ab - Whether the a- or b-wave are to be chosen. If set to 0, the a-wave, if 1 the b-wave. If #
    #          left unspecified, defaults to 0.                                                        #
    #     k0 - The wave number in vacuum, given as an int, float, or complex. Useful for normalization #
    #          when needed. If unspecified, defaults to 1.                                             #
    #                                                                                                  # 
    #                                                                                                  #
    # Outputs a list of the new wave vector and electric field.                                        #
    #                                                                                                  # 
    #                                                                                                  #
    # Last Updated: June 2, 2019                                                                       #
    #                                                                                                  #
    ####################################################################################################
    
    if wave is None:
        return
    
    reals = aux_realtest(wave)
    
    kvec = wave.kvec
    efield = wave.efield
    med = wave.medium
    
    
    if ab is None:
        ab = 0
        
    if ab != 0 and ab != 1:
        ab = 0
    
    if k0 is None:
        k0 = 1
        
    if kvec is None:
        raise Exception("Must have a wave vector to rectify the mode")
        
    if med is None:
        med = numpy.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
        print("No medium found, assuming free space")
    
    ab = ab+2
    
    kvec_star = numpy.conj(kvec)
    kvec_sqr_mod = (kvec_star.T @ kvec)[0,0]
    kvec_norm = numpy.sqrt(kvec_sqr_mod.real)    
    kvec_hat = kvec/kvec_norm
    
    # get the new wave vector
    #---------------------------------------------------------------------------------------------------
    out = aux_modecalc(kvec,med,k0=k0)
    new_k_value = out[ab]
    new_kvec = new_k_value*kvec_hat
    
    
    # get the new electric field
    #---------------------------------------------------------------------------------------------------
    if efield is None:
        new_efield = aux_maxwell_eigenvec(new_kvec,med,k0,verbose=verbose)
        if reals[2]:
            new_efield = new_efield.real
    
    else:
        
        # change problem to solver coordinates
        #-----------------------------------------------------------------------------------------------
        mats = aux_coord_transform(new_kvec+efield,new_kvec,verbose=verbose)
        U = mats[0]
        Uinv = mats[1]
        new_kvec = Uinv @ new_kvec
        ef = Uinv @ efield
        new_med = Uinv @ med @ U
        

        out = aux_maxwell_eigenvec(new_kvec,new_med,k0,verbose=verbose,switch=2)
        
        if isinstance(out,list):
            candidate_1 = numpy.conj(out[0]).T
            candidate_2 = numpy.conj(out[1]).T
            proj_1 = (candidate_1 @ ef)[0,0]
            proj_2 = (candidate_2 @ ef)[0,0]
            best = numpy.argmax([numpy.abs(proj_1),numpy.abs(proj_2)])
            new_efield = out[best]
            
        else:
            new_efield = out
        
        # transform back to lab coordinates
        #-----------------------------------------------------------------------------------------------
        new_kvec = U @ new_kvec
        new_kvec = new_kvec.real
        new_efield = U @ new_efield
        if reals[1] and reals[2]:
            new_efield = new_efield.real
        
    
    
    # return the new kvec and efield
    #---------------------------------------------------------------------------------------------------
    sol = [new_kvec, new_efield]
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
#
#
#
#
#
#
#         
#    
def aux_goodtest(ob,test_type=None):
    ''' A behind-the-scenes function for testing whether wavecalc objects have proper attributes '''
    
   
    ####################################################################################################
    #                                                                                                  #
    # The goodness test                                                                                #
    #                                                                                                  #
    # INPUTS:                                                                                          #
    #        ob - The wavecalc object to be tested.                                                    #       
    # test_type - An option to specify what type the object should be, either 'wave', 'surface', or    #
    #            'medium'.                                                                             #
    #                                                                                                  # 
    #                                                                                                  #
    # Outputs True if the object is well formed, and False otherwise                                   #
    #                                                                                                  # 
    #                                                                                                  #
    # Last Updated: June 5, 2019                                                                       #
    #                                                                                                  #
    ####################################################################################################
    
    
    # Handle improper arguments
    #---------------------------------------------------------------------------------------------------
    if (not isinstance(ob,(wavecalc.classes.wave,wavecalc.classes.surface,wavecalc.classes.medium)) 
        and ob is not None):
        raise Exception("Non-wavecalc object passed as argument")
    
    
    # Handle if the test type is specified
    #---------------------------------------------------------------------------------------------------    
    types_to_test = ['wave','surface','medium']
    if test_type in types_to_test:
        if ob is None:
            return False
        elif test_type == 'wave':
            return aux_goodtest_wav(ob)
        elif test_type == 'surface':
            return aux_goodtest_surf(ob)
        else:
            return aux_goodtest_med(ob)
   

    # Handle if test type is unspecified and object is a wave
    #---------------------------------------------------------------------------------------------------     
    elif isinstance(ob,wavecalc.classes.wave):
        return aux_goodtest_wav(ob)
    
    
    # Handle if test type is unspecified and object is a surface
    #---------------------------------------------------------------------------------------------------
    elif isinstance(ob,wavecalc.classes.surface):
        return aux_goodtest_surf(ob)
    
    
    # Handle if test type is unspecified and object is a medium
    #---------------------------------------------------------------------------------------------------
    elif isinstance(ob,wavecalc.classes.medium):
            return aux_goodtest_med(ob)
    
    
    # Handle if test type is unspecified and object is None
    #---------------------------------------------------------------------------------------------------
    else:
        return True
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
#
#
#
#
#
#
#         
#    
def aux_goodtest_wav(wav):
    ''' A behind-the-scenes function for testing whether wavecalc wave objects have proper attributes '''
    
   
    ####################################################################################################
    #                                                                                                  #
    # The goodness test for waves                                                                      #
    #                                                                                                  #
    # INPUTS:                                                                                          #
    # wav - The wavecalc wave object to be tested.                                                     #       
    #                                                                                                  # 
    #                                                                                                  #
    # Outputs True if the wave object is well formed, and False otherwise.                             #
    #                                                                                                  # 
    #                                                                                                  #
    # Last Updated: June 5, 2019                                                                       #
    #                                                                                                  #
    ####################################################################################################
    
    if wav.kvec is None or (type(wav.kvec) is numpy.ndarray 
                               and numpy.shape(wav.kvec) == (3,1)):
        ktest = True
    else:
        ktest = False
        
        
    if wav.efield is None or (type(wav.efield) is numpy.ndarray 
                                 and numpy.shape(wav.efield) == (3,1)):
        etest = True
    else:
        etest = False
        
    if wav.medium is None or (type(wav.medium) is numpy.ndarray 
                                and numpy.shape(wav.medium) == (3,3)):
        mtest = True
    else:
        mtest = False
        
    good = ktest*etest*mtest
    return good
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
#
#
#
#
#
#
#         
#    
def aux_goodtest_surf(surf):
    ''' A behind-the-scenes function for testing whether wavecalc surface objects have proper attributes '''
    
   
    ####################################################################################################
    #                                                                                                  #
    # The goodness test for surfaces                                                                   #
    #                                                                                                  #
    # INPUTS:                                                                                          #
    # surf - The wavecalc surface object to be tested.                                                 #       
    #                                                                                                  # 
    #                                                                                                  #
    # Outputs True if the surface object is well formed, and False otherwise.                          #
    #                                                                                                  # 
    #                                                                                                  #
    # Last Updated: June 5, 2019                                                                       #
    #                                                                                                  #
    ####################################################################################################
    
    if surf.normal is None or (type(surf.normal) is numpy.ndarray 
                                 and numpy.shape(surf.normal) == (3,1)):
        ntest = True
    else:
        ntest = False
        
        
    if surf.out is None or (type(surf.out) is numpy.ndarray 
                              and numpy.shape(surf.out) == (3,3)):
        otest = True
    else:
        otest = False
            
        
    if surf.into is None or (type(surf.into) is numpy.ndarray 
                               and numpy.shape(surf.into) == (3,3)):
        itest = True
    else:
        itest = False
        
    if surf.coat is None or surf.coat in ['hr','HR','ar','AR']:
        ctest = True
    else:
        ctest = False
            
    good = ntest*otest*itest*ctest
    return good
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
#
#
#
#
#
#
#         
#    
def aux_goodtest_med(med):
    ''' A behind-the-scenes function for testing whether wavecalc medium objects have proper attributes '''
    
   
    ####################################################################################################
    #                                                                                                  #
    # The goodness test for media                                                                      #
    #                                                                                                  #
    # INPUTS:                                                                                          #
    # med - The wavecalc medium object to be tested.                                                   #       
    #                                                                                                  # 
    #                                                                                                  #
    # Outputs True if the medium object is well formed, and False otherwise.                           #
    #                                                                                                  # 
    #                                                                                                  #
    # Last Updated: June 5, 2019                                                                       #
    #                                                                                                  #
    ####################################################################################################
    
    if med.epsilon is None or (type(med.epsilon) is numpy.ndarray 
                                  and numpy.shape(med.epsilon) == (3,3)):
        eptest = True
    else:
        eptest = False
        
    return eptest
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
#
#
#
#
#
#
#         
#
def aux_maxwell_eigenvec(k,ep,k0,verbose=None,switch=None):
    ''' A behind-the-scenes function for finding specific eigenvectors of the Maxwell operator '''
    
   
    ####################################################################################################
    #                                                                                                  #
    # The Maxwell eigenvector finder                                                                   #
    #                                                                                                  #
    # INPUTS:                                                                                          #
    #       k - The wave vector in solver coordinates, given as a (3,1) numpy ndarray.                 #
    #      ep - The medium given as a (3,3) numpy ndarray in solver coordinates.                       #
    # verbose - If set to True, prints more information about the calculation                          #
    #  switch - An option to handle cases of eigenvalue degeneracy. If set to 0, will choose the first #
    #           eigenvector, if set to 1 will choose the second. If set to 2 will return both.         #
    #           Defaults to 0.                                                                         #
    #                                                                                                  # 
    #                                                                                                  #
    # Outputs the eigenvector corresponding to the Maxwell operator with eigenvalue k^2 as a (3,1)     # 
    # numpy ndarray.                                                                                   # 
    #                                                                                                  # 
    #                                                                                                  #
    # Last Updated: June 2, 2019                                                                       #
    #                                                                                                  #
    ####################################################################################################
    
    if (switch is None or not isinstance(switch,int) or switch > 2 or switch < 0):
        switch = 0
    
    MAXOP = k @ k.T + (k0**2)*ep
    val = (k.T @ k)[0,0]
    err = numpy.imag(val)
    val = numpy.real(val)
    
    if err > 1e-6 and verbose==True:
        print("Imaginary part of eigenvalue is very large")
    
    tol = 1e-6
    
    e_vals, e_vecs = numpy.linalg.eig(MAXOP)
    
    e_vals = numpy.real(e_vals)
    
    ell = len(e_vals)
    
    over = numpy.zeros(ell)
    under = numpy.zeros(ell)
    
    for i in range(0,ell):
        if e_vals[i] > val-tol:
            over[i] = 1
        if e_vals[i] < val+tol:
            under[i] = 1
    
    if verbose == True:
        print('eigenvalue should be:',val)
        print('eigenvalues approximated as:',e_vals)
        print('over*under = ',over*under)
    
    if numpy.sum(over*under) == 1:    
        good = numpy.where(over*under==1)[0][0]
        sol = numpy.array([e_vecs[:,good]]).T
        return sol
    elif numpy.sum(over*under) == 2:
        if switch == 2:
            goods = numpy.where(over*under==1)[0]
            good_1 = goods[0]
            good_2 = goods[1]
            sol_1 = numpy.array([e_vecs[:,good_1]]).T
            sol_2 = numpy.array([e_vecs[:,good_2]]).T
            sol = [sol_1, sol_2]
            return sol
        else:
            good = numpy.where(over*under==1)[0][switch]
            sol = numpy.array([e_vecs[:,good]]).T
            return sol
    else:
        raise Exception('numpy has failed to appropriately approximate eigenvalues')
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
#
#
#
#
#
#
#         
#    
def aux_modecalc(vector,medium,k0=None,verbose=None):
    ''' Solves the Booker quartic for unbounded media '''
    
    
    ####################################################################################################
    #                                                                                                  #
    # The modes auxiliary function.                                                                    #
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
    # Last Updated: May 25, 2019                                                                       #
    #                                                                                                  #
    ####################################################################################################
 
    
    vector_star = numpy.conj(vector)
    vec_sqr_mod = (vector_star.T @ vector)[0,0]
    vecnorm = numpy.sqrt(vec_sqr_mod.real)
    k = vector/vecnorm
    
    if k0 is None:
        print("Assuming k0 = 1")
        k0 = 1
    elif (not isinstance(k0,int) and not isinstance(k0,float)):
        raise Exception("If specified, k0 must be given as an int or a float")
    
    ID = numpy.array([[1,0,0],[0,1,0],[0,0,1]])
    

    
    ### Build the cofactor matrix
    Cof = numpy.zeros([3,3],dtype=complex)
    for i in range(3):
        for j in range(3):
            minora = numpy.delete(medium,i,axis=0)
            minor = numpy.delete(minora,j,axis=1)
            Cof[i,j] = ((-1)**(i+j))*numpy.linalg.det(minor)
            
    adj = Cof.T
    trace = adj[0,0]+adj[1,1]+adj[2,2]
    
    A = k.T @ medium @ k
    A = A[0,0]
    
    B = (k0**2)*k.T @ (adj - trace*ID) @ k
    B = B[0,0]
    
    C = (k0**4)*numpy.linalg.det(medium)
    
    qplus = (-B + numpy.sqrt((B**2)-4*A*C))/(2*A)
    qminus = (-B - numpy.sqrt((B**2)-4*A*C))/(2*A)
    
    if verbose == True:
        print("adjugate matrix : ",adj)
        print("A : ",A)
        print("B : ",B)
        print("C : ",C)
    
    kau = numpy.sqrt(qminus)
    kbu = numpy.sqrt(qplus)
    kad = -numpy.sqrt(qminus)
    kbd = -numpy.sqrt(qplus)
    
    #kau = [kau,kau*k]
    #kbu = [kbu,kbu*k]
    #kad = [kad,kad*k]
    #kbd = [kbd,kbd*k]
    
    return [kbd,kad,kau,kbu]
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
#
#
#
#
#
#
#         
# 
def aux_quarttest(coeffs,solves):
    ''' A behind-the-scenes function for testing a quartic has been properly solved '''
    
   
    ####################################################################################################
    #                                                                                                  #
    # The quartic solutions test                                                                       #
    #                                                                                                  #
    # INPUTS:                                                                                          #
    # coeffs - The list of coefficients of the quartic, in order of increasing power                   #
    # solves - The list of solutions provided by the solving method                                    # 
    #                                                                                                  # 
    #                                                                                                  #
    # Outputs True if all of the solutions are roots of the quartic, to within 1e-6.                   #
    #                                                                                                  # 
    #                                                                                                  #
    # Last Updated: May 23, 2019                                                                       #
    #                                                                                                  #
    ####################################################################################################
    
    passes = 1
    
    for s in solves:
        fo = coeffs[4]*(s)**4
        th = coeffs[3]*(s)**3
        tw = coeffs[2]*(s)**2
        on = coeffs[1]*(s)
        ze = coeffs[0]
        trial = (numpy.abs(fo+th+tw+on+ze) < 1e-6)
        passes = trial*passes
        
    return passes
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
#
#
#
#
#
#
#         
#    
def aux_realtest(ob):
    ''' A behind-the-scenes function for testing whether wavecalc objects have all real attributes '''
    
   
    ####################################################################################################
    #                                                                                                  #
    # The realness test                                                                                #
    #                                                                                                  #
    # INPUTS:                                                                                          #
    #  ob - The wavecalc object to be tested                                                           #
    #                                                                                                  # 
    #                                                                                                  #
    # Outputs True if the object is well formed, and False otherwise                                   #
    #                                                                                                  # 
    #                                                                                                  #
    # Last Updated: June 2, 2019                                                                       #
    #                                                                                                  #
    ####################################################################################################
    
    if not aux_goodtest(ob):
        raise Exception("Bad wavecalc object or non-wavecalc object passed as argument")
    
    # Handle wave instances
    #---------------------------------------------------------------------------------------------------
    elif isinstance(ob,wavecalc.classes.wave):
        if ob.kvec is None:
            ktest = True
        else:
            X = numpy.isreal(ob.kvec)
            ktest = numpy.prod(X)
        
        if ob.efield is None:
            etest = True
        else:
            X = numpy.isreal(ob.efield)
            etest = numpy.prod(X)
        
        if ob.medium is None:
            mtest = True
        else:
            X = numpy.isreal(ob.medium)
            mtest = numpy.prod(X)
        
        
        rtest = numpy.array([ktest,etest,mtest])
        return rtest
    
    # Handle surface instances
    #---------------------------------------------------------------------------------------------------
    elif isinstance(ob,wavecalc.classes.surface):
        if ob.normal is None:
            ntest = True
        else:
            X = numpy.isreal(ob.normal)
            ntest = numpy.prod(X)
        
        if ob.out is None:
            otest = True            
        else:
            X = numpy.isreal(ob.out)
            otest = numpy.prod(X)
        
        if ob.into is None:
            itest = True
        else:
            X = numpy.isreal(ob.into)
            itest = numpy.prod(X)
            
        rtest = numpy.array([ntest,otest,itest])
        return rtest
            
    # Handle media instances
    #---------------------------------------------------------------------------------------------------    
    elif isinstance(ob,wavecalc.classes.medium):
        if ob.epsilon is None:
            rtest = True
        else:
            X = numpy.isreal(ob.epsilon)
            rtest = numpy.prod(X)
        
        return numpy.array([rtest])
    
    # Handle None instances
    #---------------------------------------------------------------------------------------------------
    else:
        return True
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
#
#
#
#
#
#
#         
#
def aux_rotate_copy(ob,ang,axis,medmove=None,verbose=None):
    ''' Rotates wavecalc objects around a specified axis: 'x', 'y', or 'z' \n
        For use with rotate class methods'''
    
    
    ####################################################################################################
    #                                                                                                  #
    # The rotation function for wavecalc object methods.                                               #
    #                                                                                                  #
    # INPUTS:                                                                                          #
    #      ob - The input object to be rotated, either a wave, surface, or medium.                     #
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
    # Rotates the object in accordance with its transformation properties.                             #
    #                                                                                                  # 
    #                                                                                                  #
    # Last Updated: May 25, 2019                                                                       #
    #                                                                                                  #
    ####################################################################################################
    
    # Raise exception for invalid 'medmove' settings
    #---------------------------------------------------------------------------------------------------
    medmove_opts = set([None,'with','only','onlyinto','onlyout','into','out'])
    if medmove not in medmove_opts:
        str1 = "If specified, 'medmove' must be set to one of the following: \n"
        str2 = "'with', 'only', 'into', 'out', 'onlyinto', 'onlyout' \n "
        str3 = " Variable 'medmove' will be set to None for following calculation. "
        print(str1+str2+str3)
        medmove = None
        
        
    # Handle wave instances
    #---------------------------------------------------------------------------------------------------
    if isinstance(ob,wavecalc.classes.wave):
        if not aux_goodtest_wav(ob):
            raise Exception('Your wavecalc wave has improper attributes')
        if medmove in medmove_opts-set([None,'with','only']):
            str1 ="For wavecalc waves, if specified, 'medmove' must be set to one of the following: 'with', 'only'.\n"
            str2 = "Variable 'medmove' will be set to None for following calculation."
            print(str1+str2)
            medmove = None
        if medmove == 'only':
            ob.medium = aux_rottens(ob.medium,ang,axis)
        else:
            ob.kvec = aux_rotvec(ob.kvec,ang,axis)
            ob.efield = aux_rotvec(ob.efield,ang,axis)
            if medmove == 'with':
                ob.medium = aux_rottens(ob.medium,ang,axis)
        if verbose == True:
            print('New kvec :',ob.kvec)
            print('New efield :',ob.efield)
            print('New medium :',ob.medium)
    
    # Handle surface instances
    #---------------------------------------------------------------------------------------------------
    elif isinstance(ob,wavecalc.classes.surface):
        if not aux_goodtest_surf(ob):
            raise Exception('Your wavecalc surface has improper attributes')
        if medmove == 'only':
            ob.out = aux_rottens(ob.out,ang,axis)
            ob.into = aux_rottens(ob.into,ang,axis)
        elif medmove == 'onlyout':
            ob.out = aux_rottens(ob.out,ang,axis)
        elif medmove == 'onlyinto':
            ob.into = aux_rottens(ob.into,ang,axis)
        else:
            ob.normal = aux_rotvec(ob.normal,ang,axis)
            if medmove == 'with':
                ob.out = aux_rottens(ob.out,ang,axis)
                ob.into = aux_rottens(ob.into,ang,axis)
            elif medmove == 'out':
                ob.out = aux_rottens(ob.out,ang,axis)
            elif medmove == 'into':
                ob.into = aux_rottens(ob.into,ang,axis)
        if verbose == True:
            print('New normal :',ob.normal)
            print('New out :',ob.into)
            print('New int :',ob.out)
       
    # Handle medium instances
    #---------------------------------------------------------------------------------------------------        
    elif isinstance(ob,wavecalc.classes.medium):
        if not aux_goodtest_med(ob):
            raise Exception('Your wavecalc medium has improper attributes')
        if not medmove is None:
            str1 ="For wavecalc media, 'medmove' option has no meaning. \n"
            str2 = "Variable 'medmove' will be set to None for following calculation."
            print(str1+str2)
        ob.epsilon = aux_rottens(ob.epsilon,ang,axis)
        if verbose == True:
            print('epsilon :',ob.epsilon)
    
    # Handle unsupported object instances
    #---------------------------------------------------------------------------------------------------
    else:
        raise Exception("Argument 'ob' must be a wavecalc wave, surface, or medium")
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
#
#
#
#
#
#
#         
#    
def aux_rotmatrix(ang,axis):
    ''' A behind-the-scenes function for building rotation matrices '''
    
    
    ####################################################################################################
    #                                                                                                  #
    # The rotation matrix                                                                              #
    #                                                                                                  #
    # INPUTS:                                                                                          #
    #  ang - The angle of rotation in degrees, given as an int or a float.                             #
    # axis - The axis about which the rotation is to be performed, either 'x', 'y', or 'z'.            #
    #                                                                                                  # 
    #                                                                                                  #
    # Outputs the corresponding rotation matrix as a (3,3) numpy matrix.                               #
    #                                                                                                  # 
    #                                                                                                  #
    # Last Updated: May 21, 2019                                                                       #
    #                                                                                                  #
    ####################################################################################################
    
    if not isinstance(ang,int) and not isinstance(ang,float):
        raise Exception("Variable 'ang' must be specified as either a float or an int")
    
    if not ang is int(ang) and not ang is float(ang):
        raise Exception("Variable 'ang' must be specified as either a float or an int")
    
    pi = 3.141592653589793115997963468544185161590576171875
    angle = ang*pi/180
    
    if axis == 'x':
        return numpy.array([[1,0,0],
                           [0,numpy.cos(angle),-numpy.sin(angle)],
                           [0,numpy.sin(angle),numpy.cos(angle)]])
    elif axis == 'y':
        return numpy.array([[numpy.cos(angle),0,numpy.sin(angle)],
                            [0,1,0],
                            [-numpy.sin(angle),0,numpy.cos(angle)]])
    elif axis == 'z':
        return numpy.array([[numpy.cos(angle),-numpy.sin(angle),0],
                            [numpy.sin(angle),numpy.cos(angle),0],
                            [0,0,1]])
    else:
        raise Exception("Axis variable must be set to one of the following strings: 'x', 'y', 'z' ")
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
#
#
#
#
#
#
#         
#    
def aux_rottens(tens,ang,axis):
    ''' A behind-the-scenes function to rotate tensors around a specified axis: 'x', 'y', or 'z' '''
    
    
    ####################################################################################################
    #                                                                                                  #
    # The tensor (matrix) rotation function                                                            #
    #                                                                                                  #
    # INPUTS:                                                                                          #
    # tens - The input tensor to be rotated, given as a (3,3) numpy array.                             #
    #  ang - The angle of rotation in degrees, given as an int or a float.                             #
    # axis - The axis about which the rotation is to be performed, either 'x', 'y', or 'z'.            #
    #                                                                                                  # 
    #                                                                                                  #
    # Outputs the rotated tensor as a (3,3) numpy array.                                               #
    #                                                                                                  # 
    #                                                                                                  #
    # Last Updated: May 21, 2019                                                                       #
    #                                                                                                  #
    ####################################################################################################
    Rinv = aux_rotmatrix(-ang,axis)
    Ro = aux_rotmatrix(ang,axis)
    sol = Ro @ tens @ Rinv
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
#
#
#
#
#
#
#         
#
def aux_rotvec(vec,ang,axis):
    ''' A behind-the-scenes function to rotate vectors around a specified axis: 'x', 'y', or 'z' '''
    
    
    ####################################################################################################
    #                                                                                                  #
    # The vector rotation function                                                                     #
    #                                                                                                  #
    # INPUTS:                                                                                          #
    #  vec - The input vector to be rotated, given as a (3,1) array.                                   #
    #  ang - The angle of rotation in degrees, given as an int or a float.                             #
    # axis - The axis about which the rotation is to be performed, either 'x', 'y', or 'z'.            #
    #                                                                                                  # 
    #                                                                                                  #
    # Outputs the rotated vector as a (3,1) numpy array.                                               #
    #                                                                                                  # 
    #                                                                                                  #
    # Last Updated: May 21, 2019                                                                       #
    #                                                                                                  #
    ####################################################################################################
    
    Ro = aux_rotmatrix(ang,axis)
    sol = Ro @ vec
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
#
#
#
#
#
#
#         
#
def aux_waveinterf(k,ef,s,ep1,ep2,k0,act=None,coating=None,same=None,verbose=None):
    ''' Outputs reflection and/or transmission waves as a list of (3,1) arrays or wavecalc waves '''
    
    
    ####################################################################################################
    #                                                                                                  #
    # The wave interface function                                                                      #
    #                                                                                                  #
    # INPUTS:                                                                                          #
    #       k - The input wave vector, given as a (3,1) array.                                         #
    #      ef - The input electric field vector, given as a (3,1) array.                               #
    #       s - The surface normal vector which defines interface surface, given as a (3,1) array.     #
    #     ep1 - The dielectric tensor of the reflection medium, given as a (3,3) numpy ndarray.        #
    #     ep2 - The dielectric tensor of the transmission medium, given as a (3,3) numpy array.        #
    #      k0 - The magnitude of the wave vector in vacuum, useful for normalizing the results, given  #
    #           as an int or a float.                                                                  #
    #     act - The action of interest, either a reflection denoted by setting the variable to 'refl', #
    #           or a transmission denoted by setting the variable to 'trans'. Leaving the variable     #
    #           unspecified results in both outputs.                                                   #
    # verbose - If set to True, prints more information about the calculation.                         #
    #                                                                                                  # 
    #                                                                                                  #
    # Outputs a (4,) list of (3,1) arrays or wavecalc waves, corresponding to the output wave vectors  #
    # or output waves respectively.                                                                    #
    #                                                                                                  # 
    #                                                                                                  #
    # Last Updated: May 30, 2019                                                                       #
    #                                                                                                  #
    ####################################################################################################
    

    # Get the coordinate transform matrices
    #---------------------------------------------------------------------------------------------------    
    coord = aux_coord_transform(k,s,verbose=verbose)
    U = coord[0]
    Uinv = coord[1]
    
    
    # Transform the inputs into the solver basis coordinates
    #---------------------------------------------------------------------------------------------------
    k_p = Uinv @ k
    
    if ef is not None:
        ef_p = Uinv @ ef
    
    ep1_p = Uinv @ ep1 @ U
    ep2_p = Uinv @ ep2 @ U
    
    
    # Normalize the kx component
    #---------------------------------------------------------------------------------------------------
    kx = k_p[0,0]/k0
    
    
    # Get the solutions to the quartic
    #---------------------------------------------------------------------------------------------------
    kz_vals = aux_booker_interf(kx,ep1_p,ep2_p,verbose=verbose)


    # The four (normalized and transposed) solution wave vectors in solver coordinates
    #---------------------------------------------------------------------------------------------------
    k_alpha_p = numpy.array([[kx,0,kz_vals[1]]])
    k_beta_p = numpy.array([[kx,0,kz_vals[0]]])
    k_gamma_p = numpy.array([[kx,0,kz_vals[2]]])
    k_nu_p = numpy.array([[kx,0,kz_vals[3]]])
    
    
    # Un-normalize and transpose the wave vectors
    #---------------------------------------------------------------------------------------------------
    k_alpha_p = k0*k_alpha_p.T
    k_beta_p = k0*k_beta_p.T
    k_gamma_p = k0*k_gamma_p.T
    k_nu_p = k0*k_nu_p.T
    
    
    # Transform the solutions back to lab coordinates
    #---------------------------------------------------------------------------------------------------
    k_alpha = U @ k_alpha_p
    k_beta = U @ k_beta_p
    k_gamma = U @ k_gamma_p
    k_nu = U @ k_nu_p
    
    
    # Delete any unecessary complex datatypes
    #---------------------------------------------------------------------------------------------------
    whole_list = [k_alpha,k_beta,k_gamma,k_nu]
    for i, wvec in enumerate(whole_list):
        really = numpy.isreal(wvec).all()
        if really:
            whole_list[i] = numpy.asarray(wvec,dtype=float)
    
    
    # Verbosity
    #---------------------------------------------------------------------------------------------------
    if verbose == True:
        if act == 'refl':
            print("k_alpha = ",k_alpha.T)
            print("k_ beta = ",k_beta.T)
        elif act == 'trans':
            print("k_gamma = ",k_gamma.T)
            print("k_ nu = ",k_nu.T)
        else:
            print("k_alpha = ",k_alpha.T)
            print("k_ beta = ",k_beta.T)
            print("k_gamma = ",k_gamma.T)
            print("k_ nu = ",k_nu.T)
    
    
    # Build the arrays of solutions
    #---------------------------------------------------------------------------------------------------
    both = numpy.array(whole_list)
    refls = numpy.array(whole_list[0:2])  #[k_alpha,k_beta]
    transs = numpy.array(whole_list[2:4])  # [k_gamma,k_nu]


    # Calculate effective indices of refraction
    #---------------------------------------------------------------------------------------------------
    k_alpha_re = k_alpha.real
    k_beta_re = k_beta.real
    k_gamma_re = k_gamma.real
    k_nu_re = k_nu.real
    
    alpha_norm = numpy.sqrt((k_alpha_re.T @ k_alpha_re)[0,0])    
    beta_norm = numpy.sqrt((k_beta_re.T @ k_beta_re)[0,0])
    gamma_norm = numpy.sqrt((k_gamma_re.T @ k_gamma_re)[0,0])
    nu_norm = numpy.sqrt((k_nu_re.T @ k_nu_re)[0,0])
    
    n_alpha = alpha_norm/k0
    n_beta = beta_norm/k0
    n_gamma = gamma_norm/k0
    n_nu = nu_norm/k0
    
    
    # Get the field eigenvectors in the solver frame
    #---------------------------------------------------------------------------------------------------
    alpha_p = aux_maxwell_eigenvec(k_alpha_p,ep1_p,k0,verbose=verbose,switch=0)
    beta_p = aux_maxwell_eigenvec(k_beta_p,ep1_p,k0,verbose=verbose,switch=1)
    gamma_p = aux_maxwell_eigenvec(k_gamma_p,ep2_p,k0,verbose=verbose,switch=0)
    nu_p = aux_maxwell_eigenvec(k_nu_p,ep2_p,k0,verbose=verbose,switch=1)
    
    
    # Transform the normalized eigenvectors to the lab frame
    #---------------------------------------------------------------------------------------------------
    alpha_hat = U @ alpha_p
    beta_hat = U @ beta_p
    gamma_hat = U @ gamma_p
    nu_hat = U @ nu_p
    
    
    # Handle the case of no efield
    #---------------------------------------------------------------------------------------------------
    if ef is None:
        if act=='refl':
            if verbose==True:
                print('alpha-wave effective index=',n_alpha)
                print('beta-wave effective index=',n_beta)
            alpha_wave = wavecalc.classes.wave(kvec=k_alpha,medium=ep1,pol=alpha_hat)
            beta_wave = wavecalc.classes.wave(kvec=k_beta,medium=ep1,pol=beta_hat)
            sol = aux_check_same([alpha_wave, beta_wave],switch=same)
            return sol
        elif act=='trans':
            if verbose==True:
                print('gamma-wave effective index=',n_gamma)
                print('nu-wave effective index=',n_nu)
            gamma_wave = wavecalc.classes.wave(kvec=k_gamma,medium=ep2,pol=gamma_hat)
            nu_wave = wavecalc.classes.wave(kvec=k_nu,medium=ep2,pol=nu_hat)
            sol = aux_check_same([gamma_wave, nu_wave],switch=same)
            return sol
        else:
            if verbose==True:
                print('alpha-wave (refl) effective index=',n_alpha)
                print('beta-wave (refl) effective index=',n_beta)
                print('gamma-wave (trans) effective index=',n_gamma)
                print('nu-wave (trans) effective index=',n_nu)
            alpha_wave = wavecalc.classes.wave(kvec=k_alpha,medium=ep1,pol=alpha_hat)
            beta_wave = wavecalc.classes.wave(kvec=k_beta,medium=ep1,pol=beta_hat)
            gamma_wave = wavecalc.classes.wave(kvec=k_gamma,medium=ep2,pol=gamma_hat)
            nu_wave = wavecalc.classes.wave(kvec=k_nu,medium=ep2,pol=nu_hat)
            sol = aux_check_same([alpha_wave, beta_wave, gamma_wave, nu_wave],switch=same)
            return sol
        
    
    
    # Find the field reflection and transmission coefficients
    #---------------------------------------------------------------------------------------------------
    coeffs = aux_field_match(k_alpha_p,k_beta_p,k_gamma_p,k_nu_p,
                             alpha_p,beta_p,gamma_p,nu_p,
                             k_p,ef_p,
                             verbose=verbose)
    
    R_alpha = coeffs[0]
    R_beta = coeffs[1]
    T_gamma = coeffs[2]
    T_nu = coeffs[3]
    
    
    
    # Construct the four field vectors in the lab frame
    #---------------------------------------------------------------------------------------------------
    alpha_ef = R_alpha*alpha_hat
    beta_ef = R_beta*beta_hat
    gamma_ef = T_gamma*gamma_hat
    nu_ef = T_nu*nu_hat
    
    
    # Adjust the field vectors to account for the surface coating
    #---------------------------------------------------------------------------------------------------  
    out_fields = aux_coat_handle(Ein=ef,Ea=alpha_ef,Eb=beta_ef,Eg=gamma_ef,En=nu_ef,ep1=ep1,ep2=ep2,coat=coating)


    # Delete any unecessary complex datatypes in the field vectors
    #---------------------------------------------------------------------------------------------------
    for i, evec in enumerate(whole_list):
        really = numpy.isreal(evec).all()
        if really:
            out_fields[i] = numpy.asarray(evec,dtype=float)


    # Get the four field vectors
    #---------------------------------------------------------------------------------------------------    
    alpha_ef = out_fields[0]
    beta_ef = out_fields[1]
    gamma_ef = out_fields[2]
    nu_ef = out_fields[3]
    
    
    # Create the four new waves as wave objects
    #---------------------------------------------------------------------------------------------------
    alpha_wave = wavecalc.classes.wave(kvec=k_alpha,efield=alpha_ef,medium=ep1) 
    beta_wave = wavecalc.classes.wave(kvec=k_beta,efield=beta_ef,medium=ep1)
    gamma_wave = wavecalc.classes.wave(kvec=k_gamma,efield=gamma_ef,medium=ep2)
    nu_wave = wavecalc.classes.wave(kvec=k_nu,efield=nu_ef,medium=ep2)
    
    
    # Create the solution sets and eliminate redunancies if necessary
    #---------------------------------------------------------------------------------------------------
    both = aux_check_same([alpha_wave,beta_wave,gamma_wave,nu_wave],switch=same)
    refls = aux_check_same([alpha_wave,beta_wave],switch=same)
    transs = aux_check_same([gamma_wave,nu_wave],switch=same)
    
    if act=='refl':
        if verbose==True:
            print('alpha-wave effective index=',n_alpha)
            print('beta-wave effective index=',n_beta)
        return refls
    elif act=='trans':
        if verbose==True:
            print('gamma-wave effective index=',n_gamma)
            print('nu-wave effective index=',n_nu)
        return transs
    else:
        if verbose==True:
            print('alpha-wave (refl) effective index=',n_alpha)
            print('beta-wave (refl) effective index=',n_beta)
            print('gamma-wave (trans) effective index=',n_gamma)
            print('nu-wave (trans) effective index=',n_nu)
        return both
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
#
#
#
#
#
#
#         
#    
def aux_complex_killer(ob):
    ''' A behind-the-scenes function for recasting complexes with 0 imaginary as floats '''
    
   
    ####################################################################################################
    #                                                                                                  #
    # The coating function                                                                             #
    #                                                                                                  #
    # INPUTS:                                                                                          #
    # ob - The object to be recast, either a WaveCalc wave, surface, or medium.                        #
    #                                                                                                  # 
    #                                                                                                  #
    # Returns the object recast as a float type if all imaginary parts are zero.                       #
    #                                                                                                  # 
    #                                                                                                  #
    # Last Updated: June 7, 2019                                                                       #
    #                                                                                                  #
    #################################################################################################### 
    
    if not isinstance(ob,(wavecalc.classes.wave,wavecalc.classes.surface,wavecalc.classes.medium)):
        return ob
    
    elif isinstance(ob,wavecalc.classes.wave):
        return ob
    elif isinstance(ob,wavecalc.classes.surface):
        return ob
    else:
        return ob
       
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
#
#
#
#
#
#
#         
#
def aux_ferrari(coeffs):
    ''' A behind-the-scenes function that uses the Ferrari method to solve quartics '''
    
   
    ####################################################################################################
    #                                                                                                  #
    # The quartic solutions test                                                                       #
    #                                                                                                  #
    # INPUTS:                                                                                          #
    # coeffs - The list of coefficients of the quartic, in order of increasing power                   #
    # solves - The list of solutions provided by the solving method                                    # 
    #                                                                                                  # 
    #                                                                                                  #
    # Outputs True if all of the solutions are roots of the quartic, to within 1e-6.                   #
    #                                                                                                  # 
    #                                                                                                  #
    # Last Updated: May 23, 2019                                                                       #
    #                                                                                                  #
    ####################################################################################################
    
    return print("Not sure if this is necessary yet")
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
#
#
#
#
#
#
#         
#
def aux_root_order(lis):
    ''' Sorts roots of the Booker Quartic into [kbd, kad, kau, kbu, number_of_complex_roots] '''
    
   
    ####################################################################################################
    #                                                                                                  #
    # The root ordering function.                                                                      #
    #                                                                                                  #
    # INPUTS:                                                                                          #
    # lis - The list of roots of the Booker quartic.                                                   # 
    #                                                                                                  # 
    #                                                                                                  #
    # Outputs the sorted list of roots ordered as [kbd,kad,kau,kbu,number_of_complex_roots]            #
    #                                                                                                  # 
    #                                                                                                  #
    # Last Updated: May 23, 2019                                                                       #
    #                                                                                                  #
    ####################################################################################################

    complexes = numpy.where(numpy.abs(numpy.imag(lis))>1e-9)[0]
    if len(complexes)==0:
        sol = numpy.sort(lis)
        sol = numpy.append(sol,0)
        return sol
    elif len(complexes)==2:
        return 


