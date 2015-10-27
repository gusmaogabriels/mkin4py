# -*- coding: utf-8 -*-
from __future__ import division, absolute_import, print_function

from . import np as __np, linalg as __linalg, derivatives as __derivatives
from .params import convergence_params
from ..bases import mkmodel as __mk

def newton_type(param):
    """ lin_solver
     params(1) : param as int
         param = 1 -> Jacobian solution
         param = 2 -> returns param = 1 and Hessian Matrix
    """
    psi, jacobian, hessian = __derivatives.analytical(param)
    jacobian = __np.dot(__mk.maps['msas'],jacobian).T
    f = __np.dot(__mk.maps['msas'],psi)
    # add white noise to the Jacobian with std of the linear problem solver converfence criteria
    #jacobian += __np.random.normal(0,convergence_params['criteriaqmr'],jacobian.shape)
    # solution algorithm
    dcoverage = __np.array(__linalg.qmr(jacobian.T,-f*convergence_params['hfun'],\
    tol = convergence_params['inner_criteria'],maxiter=convergence_params['inner_convtol']))[0]
    if param == 2 and max(abs(dcoverage))<1:
        count = 0
        fhessp = lambda dconv,M : __np.dot(M.T,dconv)
        dhess = __np.empty([len(psi),len(__mk.xsurface)])
        vhess = __np.empty([len(psi),1])
        while count <= convergence_params['convtolH']:
            for i in range(0,len(psi)):
                dhess[i,:] = fhessp(dcoverage,hessian[:,:,i])
                vhess[i,:] = __np.dot(dhess[i,:],dcoverage)
            mhess = __np.dot(__mk.maps['msas'],dhess).T
            vhess = __np.dots(__mk.maps['msas'],vhess)
            dcoverage2 = __np.array(__linalg.qmr((jacobian+mhess).T,\
            -(f+__np.dot(jacobian.T,dcoverage)+0.5*vhess)))[0]
            dcoverage += dcoverage2
            count+=1
    else:
        pass
    for i in range(len(dcoverage)):
        if __np.isnan(dcoverage[i]):
            dcoverage[i] = convergence_params['delta_min']
        elif __np.isinf(dcoverage[i]):
            if dcoverage[i]>0:
                dcoverage[i] = convergence_params['delta_min']
            else:
                dcoverage[i] = -convergence_params['delta_min']
        else:
            pass
    return dcoverage
