# -*- coding: utf-8 -*-
import numpy
import wavecalc
import warnings
#from wavecalc import classes



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
    
    
    
    

def __rotvec(vec,ang,axis=None):
    ''' A behind-the-scenes function to rotate vectors around a specified axis: 'x', 'y', or 'z' '''
    
    
    ####################################################################################################
    #                                                                                                  #
    # The vector rotation function                                                                     #
    #                                                                                                  #
    # INPUTS:                                                                                          #
    #  vec - The input vector to be rotated, given as a (3,) array.                                    #
    #  ang - The angle of rotation in degrees, given as an int or a float.                             #
    # axis - The axis about which the rotation is to be performed, either 'x', 'y', or 'z'.            #
    #                                                                                                  # 
    #                                                                                                  #
    # Outputs the rotated vector as a (3,) numpy array.                                                #
    #                                                                                                  # 
    #                                                                                                  #
    # Last Updated: May 21, 2019                                                                       #
    #                                                                                                  #
    ####################################################################################################
    
    Ro = __rotmatrix(ang,axis)
    sol = Ro @ vec
    return sol
    





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
    
    
    
    
    
    
def rotate(ob,ang,axis=None,medmove=None):
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
    #                                                                                                  #
    #                                                                                                  #
    # Rotates the object in accordance with its transformation properties.                             #
    #                                                                                                  # 
    #                                                                                                  #
    # Last Updated: May 21, 2019                                                                       #
    #                                                                                                  #
    ####################################################################################################
    
    ### Raise exception for invalid 'medmove' settings
    medmove_opts = set([None,'with','only','onlyinto','onlyout','into','out'])
    if medmove not in medmove_opts:
        str1 ="If specified, 'medmove' must be set to one of the following: \n"
        str2 = " 'with', 'only', 'into', 'out', 'onlyinto', 'onlyout' "
        str3 = " Variable 'medmove' will be set to None for following calculation. "
        warnings.warn(str1+str2+str3,stacklevel=1)
        medmove = None
        
    ### Handle (3,) numpy arrays
    if isinstance(ob,numpy.ndarray) and numpy.shape(ob) == (3,1):
        ob = __rotvec(ob,ang,axis)
        
    ### Handle (3,3) numpy arrays
    if isinstance(ob,numpy.ndarray) and numpy.shape(ob) == (3,3):
        ob = __rottens(ob,ang,axis)
        
        
    ### Handle wave instances
    if isinstance(ob,wavecalc.classes.wave):
        if not __goodtest(ob):
            warnings.warn('Your wavecalc wave has improper attributes',stacklevel=1)
        if medmove in medmove_opts-set([None,'with','only']):
            str1 ="For wavecalc waves, if specified, 'medmove' must be set to one of the following: 'with', 'only'.\n"
            str2 = "Variable 'medmove' will be set to None for following calculation."
            warnings.warn(str1+str2,stacklevel=1)
            medmove = None
        if medmove == 'only':
            ob.medium = __rottens(ob.medium,ang,axis)
        else:
            ob.kvec = __rotvec(ob.kvec,ang,axis)
            ob.efield = __rotvec(ob.efield,ang,axis)
            if medmove == 'with':
                ob.medium = __rottens(ob.medium,ang,axis)
    
    ### Handle surface instances
    elif isinstance(ob,wavecalc.classes.surface):
        if not __goodtest(ob):
            warnings.warn('Your wavecalc surface has improper attributes',stacklevel=1)
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
       
    ### Handle medium instances        
    elif isinstance(ob,wavecalc.classes.medium):
        if not __goodtest(ob):
            warnings.warn('Your wavecalc medium has improper attributes',stacklevel=1)
        if medmove not in set([None]):
            str1 ="For wavecalc media, 'medmove' option has no meaning. \n"
            str2 = "Variable 'medmove' will be set to None for following calculation."
            warnings.warn(str1+str2,stacklevel=1)
        ob.epsilon = __rottens(ob.epsilon,ang,axis)
    
    ### Handle unsupported object instances
    else:
        raise Exception("Argument 'ob' must be a (3,1) numpy array, (3,3) numpy array, or a wavecalc wave, surface, or medium")
    
    
    
    

def __waveinterf(k,s,ep,k0,act=None,verbose=None):
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
    # Last Updated: May 20, 2019                                                                       #
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
    kdotS = (k.T @ s)[0,0] #numpy.sum(k*S)
    kstar = numpy.conj(k)
    knorm = numpy.sqrt((kstar.T @ k)[0,0].real)
    
    ### Check if the wave is incident on the surface #############
    if kdotS<=0:
        str1 = 'Error: Wave vector is not incident on interface'
        return print(str1)
    ##############################################################
    
    xpa = k-kdotS*S
    xpastar = numpy.conj(xpa)
    xnorm = numpy.sqrt((xpastar.T @ xpa)[0,0].real)#numpy.sqrt(numpy.sum(xpa*xpastar).real)
    
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
                newnorm = numpy.sqrt((xpastar.T @ xpa)[0,0].real)   #numpy.sqrt(numpy.sum(xpa*xpastar).real)
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
        #print('k=',k)
        print('k_hat . s_hat =',kdotS/knorm)
        #print('S=',S)
        print("k.T' =",kp.T)
        #print('xnorm=',xnorm)
        #print('xpa=',xpa)
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
    kzs = numpy.roots([epp[2,2],DELTA,SIGMA,PSI,GAMMA])
    kzs = numpy.sort(kzs)
    ###
    
    ####################################################################################################
    ####################################################################################################
    ####################################################################################################
    
    
    
    
    ### The four (normalized) solution wave vectors in solver coordinates
    kbdp = numpy.array([[kx,0,kzs[0]]])
    kadp = numpy.array([[kx,0,kzs[1]]])
    kaup = numpy.array([[kx,0,kzs[2]]])
    kbup = numpy.array([[kx,0,kzs[3]]])
    ###
    
    ### Transform the solutions back to lab coordinates
    kbd = U @ kbdp.T
    kad = U @ kadp.T #* (numpy.asmatrix(kadp).T)
    kau = U @ kaup.T
    kbu = U @ kbup.T

    kaustar = numpy.conj(kau)
    kadstar = numpy.conj(kad)
    kbustar = numpy.conj(kbu)
    kbdstar = numpy.conj(kbd)
    
    ### Build the array of solutions and un-normalize
    sol = numpy.array([kbd,kad,kau,kbu])*k0
    refls = numpy.array([kad,kbd])*k0
    transs = numpy.array([kau,kbu])*k0
    
    ### Calculate effective indices of refraction
    aunorm = numpy.sqrt((kaustar.T @ kau)[0,0].real)    # numpy.sqrt(numpy.sum(kau*kaustar).real)   # numpy.sqrt((xpastar.T @ xpa)[0,0].real)
    adnorm = numpy.sqrt((kadstar.T @ kad)[0,0].real)
    bunorm = numpy.sqrt((kbustar.T @ kbu)[0,0].real)
    bdnorm = numpy.sqrt((kbdstar.T @ kbd)[0,0].real)
    
    nau = aunorm/k0
    nad = adnorm/k0
    nbu = bunorm/k0
    nbd = bdnorm/k0
    
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






def reflect(wav,surf,med=None,verbose=None,k0=None,HR=None):
    ''' Outputs reflection waves '''
    
     
    ####################################################################################################
    #                                                                                                  #
    # The wave interface function                                                                      #
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
                print('Using the wave medium as the reflection medium')
        else:
            reflmed = surf.out
            print("Using the surface 'out' medium as the reflection medium")
    else:
        reflmed = med.epsilon
    
    refla = wavecalc.classes.wave(efield=False,medium=reflmed)
    reflb = wavecalc.classes.wave(efield=False,medium=reflmed)
    
    if k0 is None:
        k0 = 1
        print("Assuming k0 = 1")
    
    kout = __waveinterf(wav.kvec,surf.normal,reflmed,k0,act='refl',verbose=verbose)
    
    refla.kvec = kout[0]
    reflb.kvec = kout[1]
    
    ''' Now handle the E-Field '''
    
    
    sol = [refla,reflb]
    
    return sol
            
                
        









