# -*- coding: utf-8 -*-
import numpy
import wavecalc
'''
Created May 20, 2019
author: Ryan Goetz, ryan.m.goetz@gmail.com
last update: May 31, 2019 16:59 EST
'''
'''
Table of Contents:

    Foreground Functions:
        
        crash - Line 65
        
        modes - Line 168
        
        rotate - Line 269
        
        reflect - Line 412
        
        transmit - Line 518
        
    
    Background Functions:
        
        aux_booker_interf - Line 618
        
        aux_coord_trans - Line 747
        
        aux_field_match - Line 849
        
        aux_goodtest - Line 938
        
        aux_maxwell_eigenvec - Line 1052
        
        aux_modecalc - Line 1126
        
        aux_quarttest - Line 1231
        
        aux_rotate_copy - Line 1291 
        
        aux_rotmatrix - Line 1418
        
        aux_rottens - Line 1489

        aux_rotvec - Line 1541
        
        aux_waveinterf - Line 1593
           
        
    Under Construction:
        
        aux_ferrari - Line 1797
        
        aux_root_order - Line 1846


Last line check: May 31, 2019

'''



def crash(wave,surface,k0=None,verbose=None):
    ''' Outputs the reflected and transmitted waves from 'wave' incident on 'surface' '''
    
    ### Stauts: fully functional ### 
    ####################################################################################################
    #                                                                                                  #
    # The crash function.                                                                              #
    #                                                                                                  #
    # INPUTS:                                                                                          #
    #    wave - The input wave, given as a wavecalc wave object or (3,1) numpy array.                  #
    # surface - The medium of the transmission, given as a wavecalc medium object or a (3,3) numpy     #
    #           array.                                                                                 #
    #      k0 - The magnitude of the wave vector in vacuum, useful for normalizing the results, given  #
    #           as an int or a float.                                                                  #                                                                 
    # verbose - If set to True, prints more information about the calculation.                         #
    #                                                                                                  # 
    #                                                                                                  #
    # Outputs a list of four waves resulting from wave incident on surface.                            #
    #                                                                                                  # 
    #                                                                                                  #
    # Last Updated: May 31, 2019                                                                       #
    #                                                                                                  #
    ####################################################################################################

    ### Handle (3,1) arrays as wave input
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
        
        sol = aux_waveinterf(k,ef,s,ep1,ep2,k0,verbose=verbose)
        
        return sol
        
                
    ## Handle wavecalc waves as wave input
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
        
        sol = aux_waveinterf(k,ef,s,ep1,ep2,k0,verbose=verbose)
        
        return sol
                
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

    ### Handle (3,1) arrays as vector input
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
                
    ## Handle wavecalc waves as vector input
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
                
    else:
        raise Exception("ob must be given as a wavecalc wave or (3,1) numpy.ndarray")
    
    realvec = vec.real
    if not realvec is vec:
        raise Exception("Wave vector must be real")
        
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
def rotate(ob,ang,axis=None,medmove=None,verbose=None):
    ''' Rotates wavecalc objects around a specified axis: 'x', 'y', or 'z' '''
    
    
    ####################################################################################################
    #                                                                                                  #
    # The rotation function for wavecalc object classes.                                               #
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
    # Returns the rotated object in accordance with its transformation properties.                     #
    #                                                                                                  # 
    #                                                                                                  #
    # Last Updated: May 25, 2019                                                                       #
    #                                                                                                  #
    ####################################################################################################
    
    ### Raise exception for invalid 'medmove' settings
    medmove_opts = set([None,'with','only','onlyinto','onlyout','into','out'])
    if medmove not in medmove_opts:
        str1 = "If specified, 'medmove' must be set to one of the following: \n"
        str2 = "'with', 'only', 'into', 'out', 'onlyinto', 'onlyout' \n "
        str3 = " Variable 'medmove' will be set to None for following calculation. "
        print(str1+str2+str3)
        medmove = None
        
    ### Handle (3,) numpy arrays
    if isinstance(ob,numpy.ndarray) and numpy.shape(ob) == (3,1):
        return aux_rotvec(ob,ang,axis)
        
    ### Handle (3,3) numpy arrays
    if isinstance(ob,numpy.ndarray) and numpy.shape(ob) == (3,3):
        return aux_rottens(ob,ang,axis)
        
        
    ### Handle wave instances
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
    
    ### Handle surface instances
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
       
    ### Handle medium instances        
    elif isinstance(ob,wavecalc.classes.medium):
        if not aux_goodtest(ob):
            raise Exception('Your wavecalc medium has improper attributes')
        if medmove not in set([None]):
            str1 ="For wavecalc media, 'medmove' option has no meaning. \n"
            str2 = "Variable 'medmove' will be set to None for following calculation."
            print(str1+str2)
        ee = aux_rottens(ob.epsilon,ang,axis)
        if verbose == True:
            print('epsilon :',ee)
        return wavecalc.classes.medium(epsilon=ee)
    
    ### Handle unsupported object instances
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
def reflect(wav,surf,med=None,verbose=None,k0=None,HR=None):
    ''' Outputs reflection waves '''
    
     
    ####################################################################################################
    #                                                                                                  #
    # The wave reflection function.                                                                    #
    #                                                                                                  #
    # INPUTS:                                                                                          #
    #     wav - The input wave, given as a wavecalc wave object.                                       #
    #    surf - The surface normal, given as a wavecalc surface object.                                #
    #     med - The medium of the reflection, given as a wavecalc medium object.                       #
    #      k0 - The magnitude of the wave vector in vacuum, useful for normalizing the results, given  #
    #           as an int or a float.                                                                  #
    #      HR - If set to True, the interface is made to be fully reflective.                          #                                                                  
    # verbose - If set to True, prints more information about the calculation.                         #
    #                                                                                                  # 
    #                                                                                                  #
    # Outputs a list of reflection wavecalc wave objects.                                              #
    #                                                                                                  # 
    #                                                                                                  #
    # Last Updated: May 22, 2019                                                                       #
    #                                                                                                  #
    ####################################################################################################
    if not isinstance(wav,wavecalc.classes.wave):
        raise Exception("The input 'wav' must be a wavecalc wave object")
    if not isinstance(surf,wavecalc.classes.surface): 
        raise Exception("The input 'surf' must be a wavecalc surface object")
    if not aux_goodtest(wav):
        raise Exception('Your wavecalc wave has improper attributes')
    if not aux_goodtest(surf):
        raise Exception('Your wavecalc surface has improper attributes')
    if not aux_goodtest(med):
        raise Exception('Your wavecalc medium has improper attributes')
    if med is None:
        if surf.out is None:
            if wav.medium is None:
                raise Exception("No medium found within which to perform the calculation")
            else:
                reflmed = wav.medium
                if verbose is not False:
                    print('Using the wave medium as the reflection medium')
        else:
            reflmed = surf.out
            if verbose is not False:
                print("Using the surface 'out' medium as the reflection medium")
    else:
        reflmed = med.epsilon
    
    if surf.normal is None:
        raise Exception("Surface must have a normal")
        
    if wav.kvec is None:
        raise Exception("Wave must have a wave vector")
        
    if wav.efield is None:
        raise Exception("Wave must have an electric field")
    
    
    refla = wavecalc.classes.wave(efield=False,medium=reflmed)
    reflb = wavecalc.classes.wave(efield=False,medium=reflmed)
    
    if k0 is None:
        k0 = 1
        if verbose is not False:
            print("Assuming k0 = 1")
    
    kout = aux_waveinterf(wav.kvec,surf.normal,reflmed,transmed,k0,act='refl',verbose=verbose)
    
    refla.kvec = kout[0]
    reflb.kvec = kout[1]
    
    ''' Now handle the E-Field '''
    
    
    sol = [refla,reflb]
    
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
def transmit(wav,surf,med=None,verbose=None,k0=None,AR=None):
    ''' Outputs reflection waves '''
    
     
    ####################################################################################################
    #                                                                                                  #
    # The wave transmission function.                                                                  #
    #                                                                                                  #
    # INPUTS:                                                                                          #
    #     wav - The input wave, given as a wavecalc wave object.                                       #
    #    surf - The surface normal, given as a wavecalc surface object.                                #
    #     med - The medium of the transmission, given as a wavecalc medium object.                     #
    #      k0 - The magnitude of the wave vector in vacuum, useful for normalizing the results, given  #
    #           as an int or a float.                                                                  #
    #      AR - If set to True, the interface is made to be fully transmissive.                        #                                                                  
    # verbose - If set to True, prints more information about the calculation.                         #
    #                                                                                                  # 
    #                                                                                                  #
    # Outputs a list of transmission wavecalc wave objects.                                            #
    #                                                                                                  # 
    #                                                                                                  #
    # Last Updated: May 23, 2019                                                                       #
    #                                                                                                  #
    ####################################################################################################
    if not isinstance(wav,wavecalc.classes.wave):
        raise Exception("The input 'wav' must be a wavecalc wave object")
    if not isinstance(surf,wavecalc.classes.surface): 
        raise Exception("The input 'surf' must be a wavecalc surface object")
    if not aux_goodtest(wav):
        raise Exception('Your wavecalc wave has improper attributes')
    if not aux_goodtest(surf):
        raise Exception('Your wavecalc surface has improper attributes')
    if not aux_goodtest(med):
        raise Exception('Your wavecalc medium has improper attributes')
    if med is None:
        if surf.into is None:
            raise Exception("No medium found within which to perform the calculation")
        else:
            transmed = surf.into
            if verbose is not False:
                print("Using the surface 'into' medium as the transmission medium")
    else:
        transmed = med.epsilon
    
    if surf.normal is None:
        raise Exception("Surface must have a normal")
        
    if wav.kvec is None:
        raise Exception("Wave must have a wave vector")
        
    if wav.efield is None:
        raise Exception("Wave must have an electric field")
    
    transa = wavecalc.classes.wave(efield=False,medium=transmed)
    transb = wavecalc.classes.wave(efield=False,medium=transmed)
    
    if k0 is None:
        k0 = 1
        if verbose is not False:
            print("Assuming k0 = 1")
    
    kout = aux_waveinterf(wav.kvec,surf.normal,transmed,k0,act='trans',verbose=verbose)
    
    transa.kvec = kout[0]
    transb.kvec = kout[1]
    
    ''' Now handle the E-Field '''
    
    
    sol = [transa,transb]
    
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
    
    
    ### Compute the minors of epsilon 1 in the solution coordinates
    M1 = numpy.zeros([3,3],dtype=complex)
    for i in range(3):
        for j in range(3):
            minora = numpy.delete(med1,i,axis=0)
            minor = numpy.delete(minora,j,axis=1)
            M1[i,j] = numpy.linalg.det(minor)
    ###
    
    ### Compute the minors of epsilon 2 in the solution coordinates
    M2 = numpy.zeros([3,3],dtype=complex)
    for i in range(3):
        for j in range(3):
            minora = numpy.delete(med2,i,axis=0)
            minor = numpy.delete(minora,j,axis=1)
            M2[i,j] = numpy.linalg.det(minor)
    ###
    
    ### Some convenient quantities to calculate
    ID = numpy.array([[1,0,0],[0,1,0],[0,0,1]])
    
    detep1 = med1[0,0]*M1[0,0]-med1[0,1]*M1[0,1]+med1[0,2]*M1[0,2]
    sig1 = med1+med1.T
    delt1 = (med1[0,0]+med1[1,1]+med1[2,2])*ID-med1
    
    detep2 = med2[0,0]*M2[0,0]-med2[0,1]*M2[0,1]+med2[0,2]*M2[0,2]
    sig2 = med2+med2.T
    delt2 = (med2[0,0]+med2[1,1]+med2[2,2])*ID-med2
    ###
    
    ### The coefficients of the quartic in kz/k0
    DELTA1 = kx*sig1[0,2]
    SIGMA1 = (kx**2)*delt1[1,1]-(M1[0,0]+M1[1,1])
    PSI1 = (kx**3)*sig1[0,2]+kx*(M1[0,2]+M1[2,0])
    GAMMA1 = (kx**4)*med1[0,0]-(kx**2)*(M1[1,1]+M1[2,2])+detep1
    
    DELTA2 = kx*sig2[0,2]
    SIGMA2 = (kx**2)*delt2[1,1]-(M2[0,0]+M2[1,1])
    PSI2 = (kx**3)*sig2[0,2]+kx*(M2[0,2]+M2[2,0])
    GAMMA2 = (kx**4)*med2[0,0]-(kx**2)*(M2[1,1]+M2[2,2])+detep2
    ###
    
    if verbose==True:
        print('DELTA_r =',DELTA1)
        print("SIGMA_r =",SIGMA1)
        print("PSI_r =",PSI1)
        print("GAMMA_r =",GAMMA1)
        print('DELTA_t =',DELTA2)
        print("SIGMA_t =",SIGMA2)
        print("PSI_t =",PSI2)
        print("GAMMA_t =",GAMMA2)
    
    
    
    ### Roots of the quartic
    coeffs_low_to_high_1 = [GAMMA1,PSI1,SIGMA1,DELTA1,med1[2,2]]
    coeffs_low_to_high_2 = [GAMMA2,PSI2,SIGMA2,DELTA2,med2[2,2]]
    kzs1 = numpy.polynomial.polynomial.polyroots(coeffs_low_to_high_1)
    kzs2 = numpy.polynomial.polynomial.polyroots(coeffs_low_to_high_2)
    kzs1 = numpy.sort(kzs1)
    kzs2 = numpy.sort(kzs2)
    
    kzout = [kzs1[0],kzs1[1],kzs2[2],kzs2[3]]
    
    if verbose==True:
        print("Quartic roots are approximated as: ",kzout)
        
    if not aux_quarttest(coeffs_low_to_high_1,kzs1):
        raise Exception("Quartic root solver failed to sufficiently approximate reflection roots")
        
    if not aux_quarttest(coeffs_low_to_high_2,kzs2):
        raise Exception("Quartic root solver failed to sufficiently approximate transmission roots")
        
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
    # Last Updated: May 30, 2019                                                                       #
    #                                                                                                  #
    ####################################################################################################   
    
    ### Check if the wave is incident on the surface #############
    kdotS = (k.T @ s)[0,0] 
    if kdotS<=0:
        str1 = 'Wave vector is not incident on interface'
        raise Exception(str1)
    ##############################################################
    
    snorm = numpy.sqrt((s.T @ s)[0,0])
    S = s/snorm
    zp = S
    kstar = numpy.conj(k)
    knorm = numpy.sqrt((kstar.T @ k)[0,0].real)
    
    xpa = k-kdotS*S
    xpastar = numpy.conj(xpa)
    xnorm = numpy.sqrt((xpastar.T @ xpa)[0,0].real)
    
    ### Handle when k and s are parallel ############################
    if xnorm<1e-14:
        ### Consider when a component of zp is zero:
        if len(numpy.where(zp<1e-14)[0])>0:
            ### Case where 1 is zero:
            if len(numpy.where(zp<1e-14)[0])==0:
                goods = numpy.where(zp>1e-14)[0]
                xpa = numpy.zeros([3,1])
                xpa[goods[0]] = -zp[goods[1]]/zp[goods[0]]
                xpa[goods[1]] = 1
                xpastar = numpy.conj(xpa)
                newnorm = numpy.sqrt((xpastar.T @ xpa)[0,0].real)   
                xp = xpa/newnorm
            ### Other case is 2 are zero (all three can't be zero):
            else:
                good = numpy.where(zp<1e-14)[0]
                xpa = numpy.zeros([3,1])
                xpa[good[0]] = 1
                xp = xpa
                
    #################################################################
    
    else:
        xp=xpa/xnorm
        
        
    yp = numpy.cross(zp.T,xp.T).T
    U = numpy.array([[xp[0,0],yp[0,0],zp[0,0]],[xp[1,0],yp[1,0],zp[1,0]],[xp[2,0],yp[2,0],zp[2,0]]])
    Uinv = numpy.linalg.inv(U)
    
    if verbose==True:
        print('k_hat . s_hat = ',kdotS/knorm)
        print('U = ',U)
        print('Uinv = ',Uinv)
        
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
    
    
    M = numpy.zeros([4,4],dtype=complex)
    b = numpy.zeros([4],dtype=complex)
    
    M[0,0] = alpha[0,0]
    M[0,1] = beta[0,0]
    M[0,2] = -gamma[0,0]
    M[0,3] = -nu[0,0]
    
    M[1,0] = alpha[1,0]
    M[1,1] = beta[1,0]
    M[1,2] = -gamma[1,0]
    M[1,3] = -nu[1,0]
    
    M[2,0] = k_alpha[2,0]*alpha[1,0]
    M[2,1] = k_beta[2,0]*beta[1,0]
    M[2,2] = -k_gamma[2,0]*gamma[1,0]
    M[2,3] = -k_nu[2,0]*nu[1,0]
    
    M[3,0] = k_alpha[2,0]*alpha[0,0] - k_alpha[0,0]*alpha[2,0]
    M[3,1] = k_beta[2,0]*beta[0,0] - k_beta[0,0]*beta[2,0]
    M[3,2] = k_gamma[0,0]*gamma[2,0] - k_gamma[2,0]*gamma[0,0]
    M[3,3] = k_nu[0,0]*nu[2,0] - k_nu[2,0]*nu[0,0]
    
    b[0] = -Ein[0,0]
    b[1] = -Ein[1,0]
    b[2] = -kin[2,0]*Ein[1,0]
    b[3] = kin[0,0]*Ein[2,0]-kin[2,0]*Ein[0,0]
    
    sol = numpy.linalg.solve(M,b)
    
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
def aux_goodtest(ob):
    ''' A behind-the-scenes function for testing whether wavecalc objects have proper attributes '''
    
   
    ####################################################################################################
    #                                                                                                  #
    # The goodness test                                                                                #
    #                                                                                                  #
    # INPUTS:                                                                                          #
    #  ob - The wavecalc object to be tested                                                           #
    #                                                                                                  # 
    #                                                                                                  #
    # Outputs True if the object is well formed, and False otherwise                                   #
    #                                                                                                  # 
    #                                                                                                  #
    # Last Updated: May 22, 2019                                                                       #
    #                                                                                                  #
    ####################################################################################################
    if (not isinstance(ob,wavecalc.classes.wave) 
        and not isinstance(ob,wavecalc.classes.surface) 
        and not isinstance(ob,wavecalc.classes.medium)
        and ob is not None):
        raise Exception("Non-wavecalc object passed as argument")
    
    ### Handle wave instances
    elif isinstance(ob,wavecalc.classes.wave):
        if ob.kvec is None or (str(type(ob.kvec)) == "<class 'numpy.ndarray'>" 
                               and numpy.shape(ob.kvec) == (3,1)):
            ktest = True
        else:
            ktest = False
        
        
        if ob.efield is None or (str(type(ob.efield)) == "<class 'numpy.ndarray'>" 
                                 and numpy.shape(ob.efield) == (3,1)):
            etest = True
        else:
            etest = False
        
        if ob.medium is None or (str(type(ob.medium)) == "<class 'numpy.ndarray'>" 
                                and numpy.shape(ob.medium) == (3,3)):
            mtest = True
        else:
            mtest = False
        
        good = ktest*etest*mtest
        return good
    
    ### Handle surface instances
    elif isinstance(ob,wavecalc.classes.surface):
        if ob.normal is None or (str(type(ob.normal)) == "<class 'numpy.ndarray'>" 
                                 and numpy.shape(ob.normal) == (3,1)):
            ntest = True
        else:
            ntest = False
        
        
        if ob.out is None or (str(type(ob.out)) == "<class 'numpy.ndarray'>" 
                              and numpy.shape(ob.out) == (3,3)):
            otest = True
        else:
            otest = False
            
        
        if ob.into is None or (str(type(ob.into)) == "<class 'numpy.ndarray'>" 
                               and numpy.shape(ob.into) == (3,3)):
            itest = True
        else:
            itest = False
            
        good = ntest*otest*itest
        return good
            
    ### Handle media instances    
    elif isinstance(ob,wavecalc.classes.medium):
        if ob.epsilon is None or (str(type(ob.epsilon)) == "<class 'numpy.ndarray'>" 
                                  and numpy.shape(ob.epsilon) == (3,3)):
            eptest = True
        else:
            eptest = False
        
        return eptest
    
    ### Handle None instances
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
    #                                                                                                  # 
    #                                                                                                  #
    # Outputs the eigenvector corresponding to the Maxwell operator with eigenvalue k^2 as a (3,1)     # 
    # numpy ndarray.                                                                                   # 
    #                                                                                                  # 
    #                                                                                                  #
    # Last Updated: May 30, 2019                                                                       #
    #                                                                                                  #
    ####################################################################################################
    
    MAXOP = k @ k.T + (k0**2)*ep
    val = (k.T @ k)[0,0]
    tol = 1e-8
    
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
        sol = numpy.array([e_vecs[good]]).T
        return sol
    elif numpy.sum(over*under) == 2:
        good = numpy.where(over*under==1)[0][switch]
        sol = numpy.array([e_vecs[good]]).T
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
 
    
    vecnorm = numpy.sqrt(vector.T @ vector)
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
def aux_rotate_copy(ob,ang,axis=None,medmove=None,verbose=None):
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
    
    ### Raise exception for invalid 'medmove' settings
    medmove_opts = set([None,'with','only','onlyinto','onlyout','into','out'])
    if medmove not in medmove_opts:
        str1 = "If specified, 'medmove' must be set to one of the following: \n"
        str2 = "'with', 'only', 'into', 'out', 'onlyinto', 'onlyout' \n "
        str3 = " Variable 'medmove' will be set to None for following calculation. "
        print(str1+str2+str3)
        medmove = None
        
        
    ### Handle wave instances
    if isinstance(ob,wavecalc.classes.wave):
        if not aux_goodtest(ob):
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
    
    ### Handle surface instances
    elif isinstance(ob,wavecalc.classes.surface):
        if not aux_goodtest(ob):
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
       
    ### Handle medium instances        
    elif isinstance(ob,wavecalc.classes.medium):
        if not aux_goodtest(ob):
            raise Exception('Your wavecalc medium has improper attributes')
        if medmove not in set([None]):
            str1 ="For wavecalc media, 'medmove' option has no meaning. \n"
            str2 = "Variable 'medmove' will be set to None for following calculation."
            print(str1+str2)
        ob.epsilon = aux_rottens(ob.epsilon,ang,axis)
        if verbose == True:
            print('epsilon :',ob.epsilon)
    
    ### Handle unsupported object instances
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
def aux_rotmatrix(ang,axis=None):
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
def aux_rottens(tens,ang,axis=None):
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
def aux_rotvec(vec,ang,axis=None):
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
def aux_waveinterf(k,ef,s,ep1,ep2,k0,act=None,verbose=None):
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
    
    

    
    ####################################################################################################
    ####################################################################################################
    ##################################### CHANGE OF COORDINATES ########################################
    ####################################################################################################
    ####################################################################################################
    
    coord = aux_coord_transform(k,s,verbose=verbose)
    U = coord[0]
    Uinv = coord[1]
    
    k_p = Uinv @ k
    
    if ef is not None:
        ef_p = Uinv @ ef
    
    ep1_p = Uinv @ ep1 @ U
    ep2_p = Uinv @ ep2 @ U
    
    
    ####################################################################################################
    ####################################################################################################
    ###################################### SOLVE THE QUARTIC ###########################################
    ####################################################################################################
    ####################################################################################################
    
    ### Normalize the k components
    #norm = np.sqrt(kp[0]**2+kp[1]**2+kp[2]**2)
    kx = k_p[0,0]/k0
    ###
    
    kz_vals = aux_booker_interf(kx,ep1_p,ep2_p,verbose=verbose)


    ### The four (normalized and transposed) solution wave vectors in solver coordinates
    k_alpha_p = numpy.array([[kx,0,kz_vals[1]]])
    k_beta_p = numpy.array([[kx,0,kz_vals[0]]])
    k_gamma_p = numpy.array([[kx,0,kz_vals[2]]])
    k_nu_p = numpy.array([[kx,0,kz_vals[3]]])
    ###
    
    ### Un-normalize and transpose the wave vectors
    k_alpha_p = k0*k_alpha_p.T
    k_beta_p = k0*k_beta_p.T
    k_gamma_p = k0*k_gamma_p.T
    k_nu_p = k0*k_nu_p.T
    
    
    ### Transform the solutions back to lab coordinates
    k_alpha = U @ k_alpha_p
    k_beta = U @ k_beta_p
    k_gamma = U @ k_gamma_p
    k_nu = U @ k_nu_p
    
    ### Build the array of solutions and un-normalize
    both = numpy.array([k_alpha,k_beta,k_gamma,k_nu])*k0
    refls = numpy.array([k_alpha,k_beta])*k0
    transs = numpy.array([k_gamma,k_nu])*k0
    
    ### Calculate effective indices of refraction
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
    
    if ef is None:
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
        
    
    alpha_p = aux_maxwell_eigenvec(k_alpha_p,ep1_p,k0,verbose=verbose,switch=0)
    beta_p = aux_maxwell_eigenvec(k_beta_p,ep1_p,k0,verbose=verbose,switch=1)
    gamma_p = aux_maxwell_eigenvec(k_gamma_p,ep2_p,k0,verbose=verbose,switch=0)
    nu_p = aux_maxwell_eigenvec(k_nu_p,ep2_p,k0,verbose=verbose,switch=1)
    
    
    coeffs = aux_field_match(k_alpha_p,k_beta_p,k_gamma_p,k_nu_p,
                             alpha_p,beta_p,gamma_p,nu_p,
                             k_p,ef_p,
                             verbose=verbose)
    
    R_alpha = coeffs[0]
    R_beta = coeffs[1]
    T_gamma = coeffs[2]
    T_nu = coeffs[3]
    
    alpha_ef = R_alpha*(U @ alpha_p)
    beta_ef = R_beta*(U @ beta_p)
    gamma_ef = T_gamma*(U @ gamma_p)
    nu_ef = T_nu*(U @ nu_p)
    
    
    alpha_wave = wavecalc.classes.wave(k_alpha,alpha_ef,ep1) 
    beta_wave = wavecalc.classes.wave(k_beta,beta_ef,ep1)
    gamma_wave = wavecalc.classes.wave(k_gamma,gamma_ef,ep2)
    nu_wave = wavecalc.classes.wave(k_nu,nu_ef,ep2)
    
    both = [alpha_wave,beta_wave,gamma_wave,nu_wave]
    refls = [alpha_wave,beta_wave]
    transs = [gamma_wave,nu_wave]
    
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
        sol


