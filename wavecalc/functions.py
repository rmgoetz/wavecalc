# -*- coding: utf-8 -*-
import numpy
import wavecalc
#import warnings
#from wavecalc import classes


'''
Table of Contents:

    Foreground Functions:
        
        modes - Line 56
        
        rotate - Line 157
        
        reflect - Line 300
        
        transmit - Line 406
        
    
    Background Functions:
        
        __goodtest - Line 506
        
        __modecalc - Line 620
        
        __quarttest - Line 725
        
        rotate_copy - Line 785
        
        __rotmatrix - Line 912
        
        __rottens - Line 983

        __rotvec - Line 1035
        
        __waveinterf - Line 1087
        
        
        
        
    Under Construction:
        
        __ferrari - Line 1351
        
        __root_order - Line 1400


Last line check: May 26, 2019

'''



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
        return __modecalc(vec,medi,k0,verbose)
#
#
#
#
#
#         
#
#
#
#
#
#
#         
#
#
#
#
#
#
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
        return __rotvec(ob,ang,axis)
        
    ### Handle (3,3) numpy arrays
    if isinstance(ob,numpy.ndarray) and numpy.shape(ob) == (3,3):
        return __rottens(ob,ang,axis)
        
        
    ### Handle wave instances
    if isinstance(ob,wavecalc.classes.wave):
        if not __goodtest(ob):
            raise Exception('Your wavecalc wave has improper attributes')
        if medmove in medmove_opts-set([None,'with','only']):
            str1 ="For wavecalc waves, if specified, 'medmove' must be set to one of the following: 'with', 'only'.\n"
            str2 = "Variable 'medmove' will be set to None for following calculation."
            print(str1+str2)
            medmove = None
        if medmove == 'only':
            kk = ob.kvec
            ee = ob.efield
            mm = __rottens(ob.medium,ang,axis)
        else:
            kk = __rotvec(ob.kvec,ang,axis)
            ee = __rotvec(ob.efield,ang,axis)
            mm = ob.medium
            if medmove == 'with':
                mm = __rottens(ob.medium,ang,axis)
        if verbose == True:
            print('New kvec :',kk)
            print('New efield :',ee)
            print('New medium :',mm)
        return wavecalc.classes.wave(kvec=kk,efield=ee,medium=mm)
    
    ### Handle surface instances
    elif isinstance(ob,wavecalc.classes.surface):
        if not __goodtest(ob):
            raise Exception('Your wavecalc surface has improper attributes')
        nn = ob.normal
        oo = ob.out
        ii = ob.into
        if medmove == 'only':
            oo = __rottens(ob.out,ang,axis)
            ii = __rottens(ob.into,ang,axis)
        elif medmove == 'onlyout':
            oo = __rottens(ob.out,ang,axis)
        elif medmove == 'onlyinto':
            ii = __rottens(ob.into,ang,axis)
        else:
            nn = __rotvec(ob.normal,ang,axis)
            if medmove == 'with':
                oo = __rottens(ob.out,ang,axis)
                ii = __rottens(ob.into,ang,axis)
            elif medmove == 'out':
                oo = __rottens(ob.out,ang,axis)
            elif medmove == 'into':
                ii = __rottens(ob.into,ang,axis)
        if verbose == True:
            print('New normal :',nn)
            print('New out :',ii)
            print('New int :',oo)
        return wavecalc.classes.surface(normal=nn,into=ii,out=oo)
       
    ### Handle medium instances        
    elif isinstance(ob,wavecalc.classes.medium):
        if not __goodtest(ob):
            raise Exception('Your wavecalc medium has improper attributes')
        if medmove not in set([None]):
            str1 ="For wavecalc media, 'medmove' option has no meaning. \n"
            str2 = "Variable 'medmove' will be set to None for following calculation."
            print(str1+str2)
        ee = __rottens(ob.epsilon,ang,axis)
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
    if not __goodtest(wav):
        raise Exception('Your wavecalc wave has improper attributes')
    if not __goodtest(surf):
        raise Exception('Your wavecalc surface has improper attributes')
    if not __goodtest(med):
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
    
    kout = __waveinterf(wav.kvec,surf.normal,reflmed,k0,act='refl',verbose=verbose)
    
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
    if not __goodtest(wav):
        raise Exception('Your wavecalc wave has improper attributes')
    if not __goodtest(surf):
        raise Exception('Your wavecalc surface has improper attributes')
    if not __goodtest(med):
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
    
    kout = __waveinterf(wav.kvec,surf.normal,transmed,k0,act='trans',verbose=verbose)
    
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
def __goodtest(ob):
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
def __modecalc(vector,medium,k0=None,verbose=None):
    ''' Solves the Booker quartic for  '''
    
    
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
def __quarttest(coeffs,solves):
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
def rotate_copy(ob,ang,axis=None,medmove=None,verbose=None):
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
        if not __goodtest(ob):
            raise Exception('Your wavecalc wave has improper attributes')
        if medmove in medmove_opts-set([None,'with','only']):
            str1 ="For wavecalc waves, if specified, 'medmove' must be set to one of the following: 'with', 'only'.\n"
            str2 = "Variable 'medmove' will be set to None for following calculation."
            print(str1+str2)
            medmove = None
        if medmove == 'only':
            ob.medium = __rottens(ob.medium,ang,axis)
        else:
            ob.kvec = __rotvec(ob.kvec,ang,axis)
            ob.efield = __rotvec(ob.efield,ang,axis)
            if medmove == 'with':
                ob.medium = __rottens(ob.medium,ang,axis)
        if verbose == True:
            print('New kvec :',ob.kvec)
            print('New efield :',ob.efield)
            print('New medium :',ob.medium)
    
    ### Handle surface instances
    elif isinstance(ob,wavecalc.classes.surface):
        if not __goodtest(ob):
            raise Exception('Your wavecalc surface has improper attributes')
        if medmove == 'only':
            ob.out = __rottens(ob.out,ang,axis)
            ob.into = __rottens(ob.into,ang,axis)
        elif medmove == 'onlyout':
            ob.out = __rottens(ob.out,ang,axis)
        elif medmove == 'onlyinto':
            ob.into = __rottens(ob.into,ang,axis)
        else:
            ob.normal = __rotvec(ob.normal,ang,axis)
            if medmove == 'with':
                ob.out = __rottens(ob.out,ang,axis)
                ob.into = __rottens(ob.into,ang,axis)
            elif medmove == 'out':
                ob.out = __rottens(ob.out,ang,axis)
            elif medmove == 'into':
                ob.into = __rottens(ob.into,ang,axis)
        if verbose == True:
            print('New normal :',ob.normal)
            print('New out :',ob.into)
            print('New int :',ob.out)
       
    ### Handle medium instances        
    elif isinstance(ob,wavecalc.classes.medium):
        if not __goodtest(ob):
            raise Exception('Your wavecalc medium has improper attributes')
        if medmove not in set([None]):
            str1 ="For wavecalc media, 'medmove' option has no meaning. \n"
            str2 = "Variable 'medmove' will be set to None for following calculation."
            print(str1+str2)
        ob.epsilon = __rottens(ob.epsilon,ang,axis)
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
def __rotmatrix(ang,axis=None):
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
def __rottens(tens,ang,axis=None):
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
    Rinv = __rotmatrix(-ang,axis)
    Ro = __rotmatrix(ang,axis)
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
def __rotvec(vec,ang,axis=None):
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
    
    Ro = __rotmatrix(ang,axis)
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
def __waveinterf(k,s,ep,k0,act=None,verbose=None,roots=None):
    ''' Outputs reflection or transmission wave vectors as an array of (3,1) arrays '''
    
    
    ####################################################################################################
    #                                                                                                  #
    # The wave interface function                                                                      #
    #                                                                                                  #
    # INPUTS:                                                                                          #
    #       k - The input wave vector, given as a (3,1) array.                                         #
    #       s - The surface normal vector which defines interface surface, given as a (3,1) array.     #
    #      ep - The dielectric tensor of either the reflection or transmission medium, given as a      #
    #           (3,3) numpy array.                                                                     # 
    #      k0 - The magnitude of the wave vector in vacuum, useful for normalizing the results, given  #
    #           as an int or a float.                                                                  #
    #     act - The action of interest, either a reflection denoted by setting the variable to 'refl', #
    #           or a transmission denoted by setting the variable to 'trans'. Leaving the variable     #
    #           unspecified results in both outputs, but note that they are not simultaneously         #
    #           meaningful.                                                                            #
    # verbose - If set to True, prints more information about the calculation.                         #
    #                                                                                                  # 
    #                                                                                                  #
    # Outputs either a (2,) or (4,) array of (3,1) arrays, corresponding to the output wave vectors.   #
    #                                                                                                  # 
    #                                                                                                  #
    # Last Updated: May 26, 2019                                                                       #
    #                                                                                                  #
    ####################################################################################################
    
    

    
    ####################################################################################################
    ####################################################################################################
    ##################################### CHANGE OF COORDINATES ########################################
    ####################################################################################################
    ####################################################################################################
 
    snorm = numpy.sqrt((s.T @ s)[0,0])
    S = s/snorm
    zp = S
    kdotS = (k.T @ s)[0,0] 
    kstar = numpy.conj(k)
    knorm = numpy.sqrt((kstar.T @ k)[0,0].real)
    
    ### Check if the wave is incident on the surface #############
    if kdotS<=0:
        str1 = 'Error: Wave vector is not incident on interface'
        raise Exception(str1)
    ##############################################################
    
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
    kp = Uinv @ k
    
    
    ### The transformed dielectric tensor ###
    epp = Uinv @ ep @ U
    #########################################
    
    ####################################################################################################
    ####################################################################################################
    ####################################################################################################
    
    if verbose==True:
        print('k_hat . s_hat =',kdotS/knorm)
        print("k.T' =",kp.T)
        print("x.T' =",xp.T)
        print("y.T' =",yp.T)
        print("z.T' =",zp.T)
        print('U =',U)
    
    
    
    ####################################################################################################
    ####################################################################################################
    ###################################### SOLVE THE QUARTIC ###########################################
    ####################################################################################################
    ####################################################################################################
    
    ### Normalize the k components
    #norm = np.sqrt(kp[0]**2+kp[1]**2+kp[2]**2)
    kx = kp[0,0]/k0
    ###
    
    ### Compute the minors of epsilon in the solution coordinates
    M = numpy.zeros([3,3],dtype=complex)
    for i in range(3):
        for j in range(3):
            minora = numpy.delete(epp,i,axis=0)
            minor = numpy.delete(minora,j,axis=1)
            M[i,j] = numpy.linalg.det(minor)
    ###
    
    ### Some convenient quantities to calculate
    ID = numpy.array([[1,0,0],[0,1,0],[0,0,1]])
    detep = epp[0,0]*M[0,0]-epp[0,1]*M[0,1]+epp[0,2]*M[0,2]
    sig = epp+epp.T
    delt = (epp[0,0]+epp[1,1]+epp[2,2])*ID-epp
    ###
    
    ### The coefficients of the quartic in kz/k0
    DELTA = kx*sig[1,2]
    SIGMA = (kx**2)*delt[1,1]-(M[0,0]+M[1,1])
    PSI = (kx**3)*sig[0,2]+kx*(M[0,2]+M[2,0])
    GAMMA = (kx**4)*epp[0,0]-(kx**2)*(M[1,1]+M[2,2])+detep
    ###
    
    if verbose==True:
        print('DELTA =',DELTA)
        print("SIGMA =",SIGMA)
        print("PSI =",PSI)
        print("GAMMA =",GAMMA)
    
    
    
    ### Roots of the quartic
    coeffs_low_to_high = [GAMMA,PSI,SIGMA,DELTA,epp[2,2]]
    kzs = numpy.polynomial.polynomial.polyroots(coeffs_low_to_high)
    kzs = numpy.sort(kzs)
    ###
    
    if verbose==True:
        print("Booker coefficients low to high order: ",coeffs_low_to_high)
        print("Quartic roots are approximated as: ",kzs)
        
        
    if not __quarttest(coeffs_low_to_high,kzs):
        raise Exception("Quartic root solver failed to sufficiently approximate roots")
    
    
    ####################################################################################################
    ####################################################################################################
    ####################################################################################################
    
    
    
    
    ### The four (normalized) solution wave vectors in solver coordinates
    kbdp = numpy.array([[kx,0,kzs[0]]])
    kadp = numpy.array([[kx,0,kzs[1]]])
    kaup = numpy.array([[kx,0,kzs[2]]])
    kbup = numpy.array([[kx,0,kzs[3]]])
    ###
    
   # if verbose==True:
   #     print("kbd'",kbdp)
   #     print("kad'",kadp)
   #     print("kau'",kaup)
   #     print("kbu'",kbup)
    
    ### Transform the solutions back to lab coordinates
    kbd = U @ kbdp.T
    kad = U @ kadp.T 
    kau = U @ kaup.T
    kbu = U @ kbup.T
    

    #kaustar = numpy.conj(kau)
    #kadstar = numpy.conj(kad)
    #kbustar = numpy.conj(kbu)
    #kbdstar = numpy.conj(kbd)
    
    ### Build the array of solutions and un-normalize
    sol = numpy.array([kbd,kad,kau,kbu])*k0
    refls = numpy.array([kad,kbd])*k0
    transs = numpy.array([kau,kbu])*k0
    
    ### Calculate effective indices of refraction
    kaure = kau.real
    kadre = kad.real
    kbure = kbu.real
    kbdre = kbd.real
    
    
    aunorm = numpy.sqrt((kaure.T @ kaure)[0,0])    
    adnorm = numpy.sqrt((kadre.T @ kadre)[0,0])
    bunorm = numpy.sqrt((kbure.T @ kbure)[0,0])
    bdnorm = numpy.sqrt((kbdre.T @ kbdre)[0,0])
    
    nau = aunorm/k0
    nad = adnorm/k0
    nbu = bunorm/k0
    nbd = bdnorm/k0
    
    if roots=='only':
        return kzs
    
    if act=='refl':
        if verbose==True:
            print('a-wave effective index=',nad)
            print('b-wave effective index=',nbd)
        return refls
    elif act=='trans':
        if verbose==True:
            print('a-wave effective index=',nau)
            print('b-wave effective index=',nbu)
        return transs
    else:
        if verbose==True:
            print('a-wave refl effective index=',nad)
            print('b-wave refl effective index=',nbd)
            print('a-wave trans effective index=',nau)
            print('b-wave trans effective index=',nbu)
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
def __ferrari(coeffs):
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
    
    return print("Not if this is necessary yet")
#
#
#
#
#
#         
#
#
#
#
#
#
#         
#
#
#
#
#
#
#         
#
#
#
#
#
#
#         
#
def __root_order(lis):
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


