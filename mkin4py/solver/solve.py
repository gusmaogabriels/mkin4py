# -*- coding: utf-8 -*-
from __future__ import division, absolute_import, print_function

from . import derivatives as __derivatives, coverage as __cfix, linsolver as __linsolver, np as __np, dc as __dc, time as __time
from .params import convergence_params
from ..bases import mkmodel as __mk

def rk4(param=1):
        """ Optimitzer Handler
         params(1) : param as int
             param = 1 -> Jacobian solution
             param = 2 -> returns param = 1 and Hessian Matrix
        """
        t0 = __time()
        dcoverage = __np.ones([len(__mk.maps['xsurface'])])
        dcoverage0 = dcoverage # store initial __mk.coverage for recursive comparispn
        count = 0 # total iteration counter
        restarts = 0 # total restarts counter
        psi = __derivatives.analytical(0)
        rate_sol_elem = __np.zeros([len(__mk.coverage),1])
        rate_sol = __np.dot(__mk.maps['msa'],psi)
        dcoverage[__np.isnan(dcoverage)] = convergence_params['delta_min']
        h = __dc(convergence_params['h']) # reinitialize integration step-size
        # iterative loop
        while any(abs(dcoverage-dcoverage0)>convergence_params['criteria']**2) or\
        any(abs(rate_sol[__mk.maps['surface']])>convergence_params['criteria']) or \
        min(abs(rate_sol[__mk.maps['stoichs']]))<max(abs(rate_sol[__mk.maps['surface']])):
            # 4th order RK method on the linear solver
            dcoverage0 = __dc(dcoverage)
            dcovrk = __np.zeros([len(__mk.maps['xsurface']),4])
            covrk = __dc(__mk.coverage)
            dcovrk[:,0] = __linsolver.newton_type(param)
            __mk.coverage = __cfix.coverage_update(dcovrk[:,0],0.5,covrk) # Change __mk.coverage for next step calculation
            dcovrk[:,1] = __linsolver.newton_type(param)
            __mk.coverage = __cfix.coverage_update(dcovrk[:,1],0.5,covrk)
            dcovrk[:,2] = __linsolver.newton_type(param)
            __mk.coverage = __cfix.coverage_update(dcovrk[:,2],1,covrk)
            dcovrk[:,3] = __linsolver.newton_type(param)
            dcoverage = (1/6.0)*__np.sum(dcovrk[:,[0,3]],1)+(1/3.0)*(__np.sum(dcovrk[:,[1,2]],1))
            __mk.coverage = __cfix.coverage_update(dcoverage,h,covrk)
            if any(__np.isnan(__mk.coverage)):
                __mk.__init___mk.coverage()
            else:
                pass
            psi = __derivatives.analytical(param=0)
            rate_sol = __np.dot(__mk.maps['msa'],psi)
            rate_sol_elem = __np.array(__mk.kinetic_parameters['k'],ndmin=2).T*__np.array(psi,ndmin=2)
            #print __np.sqrt(sum(rate_sol[surface]**2))
            #print __np.concatenate((__np.array(MKmodel.splabels,ndmin=2).T,MKmodel.ms.dot(psi)),axis=1)
            #print __mk.coverage[surface]
            count += 1
            if count>=convergence_params['convtol'] or __time()-t0>convergence_params['max_time']:
                msg = 'Convergence NOT achieved.'
                __mk.init_coverage()
                count = 0
                restarts += 1
                if restarts > convergence_params['max_restarts']:
                    return {'coverage':__np.array(__mk.coverage,ndmin=2).T,'rates':rate_sol,\
                    'elem_rate':rate_sol_elem,'msg':msg,'time':__time()-t0}
                else:
                    pass
        msg = 'Convergence achieved'
        return {'coverage':__np.array(__mk.coverage,ndmin=2).T,'rates':rate_sol,\
                    'elem_rate':rate_sol_elem,'msg':msg,'time':__time()-t0}
