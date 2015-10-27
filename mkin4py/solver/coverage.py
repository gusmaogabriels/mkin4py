# -*- coding: utf-8 -*-

from __future__ import division, absolute_import, print_function
from . import dc as __dc
from .params import convergence_params
from ..bases import mkmodel as __mk

def coverage_update(delta,h=1,cov=[]):
    """ coverage_update
    params(2) :
        delta coverage -> numeric of size of xsurface
        h = stepsize -> numeric value
        coverage -> reference to __mk.coverage or equivalent __mk.coverage
    return a adjusted version of the __mk.coverage
    """
    if cov == []:
        __mk.coverage = __dc(__mk.coverage)
    else:
        __mk.coverage = __dc(cov)
    __mk.coverage[__mk.maps['xsurface']] += h*delta # Iteration update
    __mk.coverage[__mk.maps['ndof']] -= h*sum(delta) # DOF Site iteration update
    __mk.coverage[__mk.maps['surface']] = __mk.coverage[__mk.maps['surface']]*(__mk.coverage[__mk.maps['surface']]>0)\
        + convergence_params['delta_min']*(__mk.coverage[__mk.maps['surface']]<=0)
    __mk.coverage[__mk.maps['surface']] /= sum(__mk.coverage[__mk.maps['surface']]) # Renormalization
    return __mk.coverage
