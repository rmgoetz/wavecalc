# -*- coding: utf-8 -*-
import numpy
from math import pi




def rotvec(vec,ang,axis=None):
    
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
    # Outputs the rotated vector as a (3,1) array.                                                     #
    #                                                                                                  # 
    #                                                                                                  #
    # Last Updated: May 20, 2019                                                                       #
    #                                                                                                  #
    ####################################################################################################
    
    angle = ang*pi/180
    if axis =='x':
        Ro = numpy.matrix([[1,0,0],
                           [0,numpy.cos(angle),-numpy.sin(angle)],
                           [0,numpy.sin(angle),numpy.cos(angle)]])
        mvec = numpy.transpose(numpy.asmatrix(vec))
        rot = Ro*mvec
        sol = numpy.asarray(numpy.transpose(rot))[0]
        return sol
    elif axis =='y':
        Ro = numpy.matrix([[numpy.cos(angle),0,numpy.sin(angle)],
                            [0,1,0],
                            [-numpy.sin(angle),0,numpy.cos(angle)]])
        mvec = numpy.transpose(numpy.asmatrix(vec))
        rot = Ro*mvec
        sol = numpy.asarray(numpy.transpose(rot))[0]
        return sol
    elif axis =='z':
        Ro = numpy.matrix([[numpy.cos(angle),-numpy.sin(angle),0],
                            [numpy.sin(angle),numpy.cos(angle),0],
                            [0,0,1]])
        mvec = numpy.transpose(numpy.asmatrix(vec))
        rot = Ro*mvec
        sol = numpy.asarray(numpy.transpose(rot))[0]
        return sol
    else:
        str1 = 'Error: axis variable must be set to one of the following strings: \n'
        str2 = " 'x', 'y', 'z' "
        return print(str1+str2)
    
    
    
    
    
    
    
    
    

def waveinterf(k,s,ep,k0,act=None,verbose=None):
    
    ####################################################################################################
    #                                                                                                  #
    # The wave interface function                                                                      #
    #                                                                                                  #
    # INPUTS:                                                                                          #
    #       k - The input wave vector, given as a (3,) array.                                          #
    #       s - The surface normal vector which defines interface surface, given as a (3,) array.      #
    #      ep - The dielectric tensor of either the reflection or transmission medium, given as a      #
    #           (3,3) numpy matrix.                                                                    # 
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
 
    snorm = numpy.sqrt(numpy.sum(s*s))
    kstar = numpy.conj(k)
    S = s/snorm
    zp = S
    kdotS = numpy.sum(k*S)
    
    ### Check if the wave is incident on the surface #############
    if kdotS<=0:
        str1 = 'Error: Wave vector is not incident on interface '
        return print(str1)
    ##############################################################
    
    xpa = k-kdotS*S
    xpastar = numpy.conj(xpa)
    knorm = numpy.sqrt(numpy.real(numpy.sum(k*kstar)))
    xnorm = numpy.sqrt(numpy.real(numpy.sum(xpa*xpastar)))
    
    ### Handle when k and s are parallel ############################
    if xnorm<1e-14:
        ### Consider when a component of zp is zero:
        if len(numpy.where(zp<1e-14)[0])>0:
            ### Case where 1 is zero:
            if len(numpy.where(zp<1e-14)[0])==0:
                goods = numpy.where(zp>1e-14)[0]
                xpa = numpy.zeros([3])
                xpa[goods[0]] = -zp[goods[1]]/zp[goods[2]]
                xpa[goods[1]] = 1
                xpastar = numpy.conj(xpa)
                newnorm = numpy.sqrt(numpy.real(numpy.sum(xpa*xpastar)))
                xp = xpa/newnorm
            ### Other case is 2 are zero (all three can't be zero):
            else:
                good = numpy.where(zp<1e-14)[0]
                xpa = numpy.zeros([3])
                xpa[good[0]] = 1
                xp = xpa
                
    #################################################################
    
    else:
        xp=xpa/xnorm
        
        
    yp = numpy.cross(zp,xp)
    U = numpy.matrix([[xp[0],yp[0],zp[0]],[xp[1],yp[1],zp[1]],[xp[2],yp[2],zp[2]]])
    Uin = numpy.linalg.inv(U)
    kpa = Uin*numpy.transpose(numpy.asmatrix(k))
    
    ### The transformed k array ##########
    kp = numpy.asarray(numpy.transpose(kpa))[0]
    ######################################
    
    ### The transformed dielectric tensor ###
    epp = Uin*ep*U
    #########################################
    
    ####################################################################################################
    ####################################################################################################
    ####################################################################################################
    
    if verbose==True:
        #print('k=',k)
        print('k . S =',kdotS)
        #print('S=',S)
        print("k' =",kp)
        #print('xnorm=',xnorm)
        #print('xpa=',xpa)
        print("x' =",xp)
        print("y' =",yp)
        print("z' =",zp)
        print('U =',U)
    
    
    
    ####################################################################################################
    ####################################################################################################
    ###################################### SOLVE THE QUARTIC ###########################################
    ####################################################################################################
    ####################################################################################################
    
    ### Normalize the k components
    #norm = np.sqrt(kp[0]**2+kp[1]**2+kp[2]**2)
    kx = kp[0]/k0
    ###
    
    ### Compute the minors of epsilon in the solution coordinates
    M = numpy.asmatrix(numpy.zeros([3,3]),dtype=complex)
    for i in range(3):
        for j in range(3):
            minora = numpy.delete(epp,i,axis=0)
            minor = numpy.delete(minora,j,axis=1)
            M[i,j] = numpy.linalg.det(minor)
    ###
    
    ### Some convenient quantities to calculate
    ID = numpy.matrix([[1,0,0],[0,1,0],[0,0,1]])
    detep = epp[0,0]*M[0,0]-epp[0,1]*M[0,1]+epp[0,2]*M[0,2]
    sig = epp+numpy.transpose(epp)
    delt = (epp[0,0]+epp[1,1]+epp[2,2])*ID-epp
    ###
    
    ### The coefficients of the quartic in kz/k0
    DELTA = kx*sig[1,2]
    SIGMA = (kx**2)*delt[1,1]-(M[0,0]+M[1,1])
    PSI = (kx**3)*sig[0,2]+kx*(M[0,2]+M[2,0])
    GAMMA = (kx**4)*epp[0,0]-(kx**2)*(M[1,1]+M[2,2])+detep
    ###
    
    ### Roots of the quartic
    kzs = numpy.roots([epp[2,2],DELTA,SIGMA,PSI,GAMMA])
    kzs = numpy.sort(kzs)
    ###
    
    ####################################################################################################
    ####################################################################################################
    ####################################################################################################
    
    
    
    
    ### The four (normalized) solution wave vectors in solver coordinates
    kbdp = numpy.array([kx,0,kzs[0]])
    kadp = numpy.array([kx,0,kzs[1]])
    kaup = numpy.array([kx,0,kzs[2]])
    kbup = numpy.array([kx,0,kzs[3]])
    ###
    
    ### Transform the solutions back to lab coordinates
    kbdm = U*numpy.transpose(numpy.asmatrix(kbdp))
    kadm = U*numpy.transpose(numpy.asmatrix(kadp))
    kaum = U*numpy.transpose(numpy.asmatrix(kaup))
    kbum = U*numpy.transpose(numpy.asmatrix(kbup))
    kbd = numpy.asarray(numpy.transpose(kbdm))[0]
    kad = numpy.asarray(numpy.transpose(kadm))[0]
    kau = numpy.asarray(numpy.transpose(kaum))[0]
    kbu = numpy.asarray(numpy.transpose(kbum))[0]
    kaustar = numpy.conj(kau)
    kadstar = numpy.conj(kad)
    kbustar = numpy.conj(kbu)
    kbdstar = numpy.conj(kbd)
    
    ### Build the array of solutions and un-normalize
    sol = numpy.array([kbd,kad,kau,kbu])*k0
    refls = numpy.array([kad,kbd])*k0
    transs = numpy.array([kau,kbu])*k0
    
    ### Calculate effective indices of refraction
    aunorm = numpy.sqrt(numpy.real(numpy.sum(kau*kaustar)))
    adnorm = numpy.sqrt(numpy.real(numpy.sum(kad*kadstar)))
    bunorm = numpy.sqrt(numpy.real(numpy.sum(kbu*kbustar)))
    bdnorm = numpy.sqrt(numpy.real(numpy.sum(kbd*kbdstar)))
    
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
