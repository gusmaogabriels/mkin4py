# -*- coding: utf-8 -*-

from __future__ import division, absolute_import, print_function
from . import np as __np
from .params import convergence_params
from ..bases import mkmodel as __mk

def analytical(param=1,cov=[]):
    """ Derivative Generator
     para__mk.ms(1) : param as int
         param = 0 -> returns reaction rates power-law on concentrations parcel (psi)
         param = 1 -> returns param = 0 and Jacobian Matrix
         param = 2 -> returns param = 1 and Hessian Matrix
    """
    if cov == []:
        pass
    else:
        __mk.coverage = cov
    pos = __mk.coverage == 0
    pos[__mk.maps['stoichs']] = 0
    __mk.coverage[pos != 0] = convergence_params['delta_min'] # Zero approximation for __mk.surface species
    __mk.coverage[__mk.maps['surface']] = __mk.coverage[__mk.maps['surface']]/sum(__mk.coverage[__mk.maps['surface']]) # Renormalization
    __mk.coverage_inv = __np.array(__mk.coverage[__mk.maps['surface']]**(-1),ndmin=2).T
    mfp = (__mk.ms<0).T # reactant species mapping in __mk.ms
    mfp_coeff = -(__mk.ms*(mfp.T)).T # reactant species power-dependency mapping in __mk.ms
    psi_q = __np.dot(mfp,__np.diag(__mk.coverage))
    for i in range(1,int(abs((__mk.ms)).max()+1)):
        psi_q *= __np.dot((__mk.ms<-i).T,__np.diag(__mk.coverage))+((__mk.ms<-i).T==0)*1
    psi = __np.empty([psi_q.shape[0],1])
    for i in range(0,psi_q.shape[0]):
        psi[i] = __np.prod(psi_q[i,mfp[i,:]==1])
    if param>=1: # Jacobian Creator
        jacobian = __np.dot(psi,__mk.coverage_inv.T)
        jacobian *= jacobian>convergence_params['delta_min']*__mk.ms.max()
        jacobian *= mfp_coeff[:,__mk.maps['surface']]
        hessian = []
        sitepos = __np.where(__mk.maps['surface']==__mk.maps['ndof'])
        otherpos = __np.array(range(0,len(__mk.maps['surface'])))[__np.array(range(0,len(__mk.maps['surface'])))!=sitepos[0]]
        if param == 2:
            hessianM = __np.empty([len(__mk.maps['surface']),len(__mk.maps['surface']),len(psi_q)])
            psi[psi==0] = convergence_params['delta_min']
            psi_inv = psi**-1
            f_hessian = lambda i : __np.dot(__np.dot(jacobian[i,:].T,psi_inv),jacobian[i,:])\
            -__np.diag(__np.dot(__np.diag(__mk.coverage_inv),jacobian[i,:].T))
            for i in range(0,len(psi)):
                hessianM[:,:,i] = f_hessian(i)
            hessian = __np.empty([hessianM.shape[0]-1,hessianM.shape[1]-1,hessianM.shape[2]])
            for i in range(0,hessian.shape[2]):
                hessianT = hessianM[:,:,0]
                sline = hessianT[sitepos,:]
                hessianT = hessianT[otherpos,:]
                hessianT -= __np.repeat(sline,hessianT.shape[0],0)
                sline = hessianT[:,sitepos]
                hessianT = hessianT[:,otherpos]
                hessianT -= __np.repeat(sline,hessianT.shape[1],1)
                hessian[:,:,i] = hessianT
                hessianM = hessianM[:,:,1::]
        sline = jacobian[:,sitepos]
        jacobian = jacobian[:,otherpos]
        jacobian -= __np.repeat(sline[:,:,0],jacobian.shape[1],1)
        return psi, jacobian, hessian
    else:
        pass
    return psi
