#!/usr/bin/env python
# -*- coding: utf-8 -*-   

############################################################
import os
import sys
import numpy as np
import multiprocessing
from vaspwfc import vaspwfc

############################################################

def orthogon(cic):
    
    S = np.dot(cic.conj(),cic.T)

    Dsqrt= np.zeros_like(S, dtype=np.complex)
    T = np.zeros_like(S, dtype=np.complex)
    cio = np.zeros_like(cic, dtype=np.complex)
    
    D,V=np.linalg.eig(S)
    for ii in range(np.size(S,0)):
        Dsqrt[ii,ii]=1/np.sqrt(D[ii])
    T=np.dot(np.dot(V.conj(),Dsqrt.conj()),V.T)
    
    cio = np.dot(T,cic)
    
    return cio

def init_wav_from_vaspwfc(wave0, gamma=True,icor=1,
                     bmin=None, bmax=None, omin=None, omax=None,
                     ikpt=1, ispin=1, dt=1.0):
    '''
    Initialize Wavefunction for phase correction
    
    inputs:
        wave0:  path of reference WAVECAR 
        gamma:  gamma version wavecar
        ikpt:   k-point index, starting from 1 to NKPTS
        ispin:  spin index, 1 or 2
        dt:     ionic time step, in [fs]          

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! Note, this method is much slower than fortran code. !!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! Now, It is much faster than fortran code :)         !!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    '''

    phi_0 = vaspwfc(wave0)

    bmin = 1 if bmin is None else bmin
    bmax = phi_0._nbands if bmax is None else bmax
    omin = bmin if icor==12 else omin
    omax = bmax if icor==12 else omax
    nbasis = bmax - bmin + 1
    obasis = omax - omin + 1    # Orthogonalized basis


    ci_t   = phi_0.readBandCoeff(ispin, ikpt, omax, norm=True)
    cic_0 = np.zeros([obasis] + list(ci_t.shape),dtype=np.complex)
    cio_0 = np.zeros([obasis] + list(ci_t.shape),dtype=np.complex)
    
    print cic_0.shape
    
    file_path = "cio.npy"    #Reference Wavefunction for phase correction
 
    if ( os.path.isfile(file_path)):
        cio_0 = np.load(file_path).view(complex)
        if( np.size(cio_0,0) != obasis):
            print("# of Bands from cio.npy not equal")
	    sys.exit(0)

    else:
        print("Reference Wavefunction not found")
        for ii in range(obasis):
            ib1 = ii + omin
            cic_0[ii,:]   = phi_0.readBandCoeff(ispin, ikpt, ib1, norm=True)
        cio_0 = orthogon(cic_0) if icor==12 else cic_0
        np.save(file_path,cio_0.view(float))
 	print(wave0,"used as reference wavefunction")



def nac_from_vaspwfc(waveA, waveB, gamma=True,icor=1,
                     bmin=None, bmax=None, omin=None, omax=None,
                     ikpt=1, ispin=1, dt=1.0):
    '''
    Calculate Nonadiabatic Couplings (NAC) from two WAVECARs
    <psi_i(t)| d/dt |(psi_j(t))> ~=~
                                    (<psi_i(t)|psi_j(t+dt)> -
                                     <psi_i(t+dt)|psi_j(t)>) / (2dt)
    In fact, the output are RAW NACs which are need to multiply by hbar/(2*dt)
    inputs:
        waveA:  path of WAVECAR A
        waveB:  path of WAVECAR B
        gamma:  gamma version wavecar
        dt:     ionic time step, in [fs]          
        ikpt:   k-point index, starting from 1 to NKPTS
        ispin:  spin index, 1 or 2

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! Note, this method is much slower than fortran code. !!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! Now, It is much faster than fortran code :)         !!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    '''

    phi_i = vaspwfc(waveA)      # wavecar at t
    phi_j = vaspwfc(waveB)      # wavecar at t + dt

    print 'Calculating NACs between <%s> and <%s>' % (waveA, waveB)

    assert phi_i._nbands == phi_j._nbands, '#bands not match!'
    assert phi_i._nplws[ikpt-1] == phi_j._nplws[ikpt-1], '#nplws not match!'

    

    bmin = 1 if bmin is None else bmin
    bmax = phi_i._nbands if bmax is None else bmax
    omin = bmin if icor!=12 else omin
    omax = bmax if icor!=12 else omax
    nbasis = bmax - bmin + 1
    obasis = omax - omin + 1    # Basis for orthogonalization

    nacType = np.float if gamma else np.complex
    nacs = np.zeros((nbasis, nbasis), dtype=nacType)
    #nacs2 = np.zeros((nbasis, nbasis), dtype=nacType)
    pij = np.zeros((nbasis, nbasis), dtype=np.complex)
    pji = np.zeros((nbasis, nbasis), dtype=np.complex)
    
    
    

    ci_t   = phi_i.readBandCoeff(ispin, ikpt, omax, norm=True)
    cic_0 = np.zeros([obasis] + list(ci_t.shape),dtype=np.complex)
    cio_0 = np.zeros([obasis] + list(ci_t.shape),dtype=np.complex)
    cic_t = np.zeros([obasis] + list(ci_t.shape),dtype=np.complex)
    cic_tdt = np.zeros([obasis] + list(ci_t.shape),dtype=np.complex)
           
    cio_t = np.zeros([obasis] + list(ci_t.shape),dtype=np.complex)
    cio_tdt = np.zeros([obasis] + list(ci_t.shape),dtype=np.complex)
    print cic_t.shape
    
    from time import time
    t1 = time()
    
    file_path = "cio.npy"    #Orthogonalized Wavefunction for time=0
 
    if (icor != 1):
        if (os.path.isfile(file_path)):
            cio_0 = np.load(file_path).view(complex)
            if (np.size(cio_0,0) != obasis):
                print("# of Bands from cio.npy not equal")
	        sys.exit(0)

        else:
            print("Reference Wavefunction not found")
	    sys.exit(0)


    for ii in range(obasis):
        ib1 = ii + omin
        cic_t[ii,:]   = phi_i.readBandCoeff(ispin, ikpt, ib1, norm=True)
        cic_tdt[ii,:]   = phi_j.readBandCoeff(ispin, ikpt, ib1, norm=True)
    
    cio_t = orthogon(cic_t) if icor==12 else cic_t
    cio_tdt = orthogon(cic_tdt) if icor==12 else cic_tdt
    #print "cio_t",cio_t.shape,np.dot(cio_t.conj()[1,:],cio_t[2,:]),np.dot(cio_t.conj()[2,:],cio_t[2,:]),np.dot(cio_t.conj()[1,:],cio_t[3,:])
    
    
    #np.savetxt('ci1.txt',np.dot(cio_tdt.conj(),cio_tdt.T))
    #np.savetxt('ci2.txt',np.dot(cio_t.conj(),cio_t.T))
    
    if icor == 1 :
        cio_0=cio_t


    cor1 =cio_0.conj()*cio_t
    cor1 =np.sum(cor1,axis=1)
    cc1= cor1/abs(cor1)
    #print "cor1",abs(cor1)

    cor2 =cio_0.conj()*cio_tdt
    cor2 =np.sum(cor2,axis=1)
    #print cor2.shape,cor2
    cc2= cor2/abs(cor2)
    
    if gamma:
    	for ii in range(obasis):
	    cc1[ii]= 1.0 if cor1[ii].real>0 else -1.0
	    cc2[ii]= 1.0 if cor2[ii].real>0 else -1.0
	    
    if icor==1:
        file_path = "cc1.npy"    #Phase correction factor for time=0
 
        if (not os.path.isfile(file_path)):
       	    np.save(file_path,cc1) 

        cc1=np.load('cc1.npy')
    	#print cc1
	cor2 =cio_t.conj()*cio_tdt
    	cor2 =np.sum(cor2,axis=1)
    	cc2= cor2/abs(cor2)
    		
    	if gamma:
    	    for ii in range(obasis):
	    	cc2[ii]= 1.0 if cor2[ii].real>0 else -1.0
	
        cc2=cc2*cc1
	np.save('cc1.npy',cc2)
	     
    for ii in range(nbasis):
        for jj in range(ii):
       # for jj in range(nbasis):

            ibi = ii + bmin - omin 
            ibj = jj + bmin - omin  

            # Non Phase Correction           

            #pij[ii,jj] = np.sum(cio_t[ibi,:].conj()*cio_tdt[ibj,:])
            #pji[ii,jj] = np.sum(cio_tdt[ibi,:].conj()*cio_t[ibj,:])
           
            #Phase Correction
            
            pij[ii,jj] = np.sum(cio_t[ibi,:].conj()*cio_tdt[ibj,:])*cc1[ibi]*cc2[ibj].conj()
            pji[ii,jj] = np.sum(cio_tdt[ibi,:].conj()*cio_t[ibj,:])*cc2[ibi]*cc1[ibj].conj()
            
	    # Four terms phase correction
	    #pij[ii,jj] = np.sum(cio_t[ibi,:].conj()*cio_tdt[ibj,:])*cc1[ibi]*cc2[ibj].conj()-np.sum(cio_t[ibi,:].conj()*cio_t[ibj,:])*cc1[ibi]*cc1[ibj].conj()
            #pji[ii,jj] = np.sum(cio_tdt[ibi,:].conj()*cio_t[ibj,:])*cc2[ibi]*cc1[ibj].conj()-np.sum(cio_tdt[ibi,:].conj()*cio_tdt[ibj,:])*cc2[ibi]*cc2[ibj].conj()
           
            
            tmp = pij[ii,jj]-pji[ii,jj]
            #tmp2 = abs(pij[ii,jj].real)-abs(pji[ii,jj].real)
            #tmp = np.sum(cic_t[ibi,:].conj() * cic_tdt[ibj,:])*cor1[ibi]*cor2.conj[ibj] - np.sum(cic_tdt[ibi,:].conj() * cic_tdt[ibj,:])
            
            nacs[ii,jj] = tmp.real if gamma else tmp
            #nacs2[ii,jj] = tmp2
            if gamma:
                nacs[jj,ii] = -nacs[ii,jj]
                #nacs2[jj,ii] = -nacs2[ii,jj]
            else:
                nacs[jj,ii] = -np.conj(nacs[ii,jj])

    t2 = time()
    print '1. Elapsed Time: %.4f [s]' % (t2 - t1)
  # EnT = (phi_i._bands[ispin-1,ikpt-1,:] + phi_j._bands[ispin-1,ikpt-1,:]) / 2.
    EnT = phi_i._bands[ispin-1,ikpt-1,bmin-1:bmax]
    
    # close the wavecar
    phi_i._wfc.close()
    phi_j._wfc.close()

    # return EnT, nacs / (2 * dt)
    return EnT, nacs,pij,pji
    #return EnT, nacs,pij,pji,nacs2


def parallel_nac_calc(runDirs, nproc=None, gamma=False,icor=1,
                      bmin=None, bmax=None,omin=None, omax=None,
                      ikpt=1, ispin=1, dt=1.0):
    '''
    Parallel calculation of NACs using python multiprocessing package.
    '''
    import multiprocessing
    import time 
    nproc = multiprocessing.cpu_count() if nproc is None else nproc
    pool = multiprocessing.Pool(processes=nproc)
    results = []
    
    for w1, w2 in zip(runDirs[:-1], runDirs[1:]):
        res = pool.apply_async(nac_from_vaspwfc, (w1, w2, gamma, icor, bmin, bmax, omin, omax, ikpt, ispin, dt))
        results.append(res)

    for ii in range(len(runDirs)-1):
        et, nc , pij, pji= results[ii].get()

        prefix = os.path.dirname(runDirs[ii])
        np.savetxt(prefix + '/eig%02d.txt' % (icor), et[np.newaxis, :])
        np.savetxt(prefix + '/nac%02d_re.txt' % (icor), nc.real.flatten()[np.newaxis, :])
        #np.savetxt(prefix + '/nacr%02d_re.txt' % (icor), nac2.real.flatten()[np.newaxis, :])
        np.savetxt(prefix + '/nac%02d_im.txt' % (icor), nc.imag.flatten()[np.newaxis, :])

############################################################
############################################################

if __name__ == '__main__':
    T_start = 1 
    T_end   = 4000 
    T_ref   = 1
    nproc   = 1
    gamma   = False
    icor    = 1 
    bmin    = 316
    bmax    = 350
    omin    = bmin
    omax    = bmax
    ikpt    = 1
    ispin   = 1

    WaveCars = ['./%04d/WAVECAR' % (ii + 1) for ii in range(T_start-1, T_end)]
    
    if T_start == T_ref:
        if os.path.exists("cio.npy"):
            os.remove("cio.npy")
        if os.path.exists("cc1.npy") and icor == 1:
            os.remove("cc1.npy")

    if icor == 1:
        if nproc != 1:
            print("This Phase correction scheme does not support parallel computing")
            sys.exit(0)
        if T_start != T_ref:
            print("Calculation is not started from the first step, be aware of what you are running!")


    if icor == 2 or icor == 12:

        wave0 = './%04d/WAVECAR' % (T_ref)
        print (wave0) 
        init_wav_from_vaspwfc(wave0, gamma, icor, bmin, bmax, omin, omax, ikpt, ispin)


    parallel_nac_calc(WaveCars, nproc, gamma, icor, bmin, bmax, omin, omax, ikpt, ispin)
   

    
