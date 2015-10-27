# -*- coding: utf-8 -*-
from __future__ import division, absolute_import, print_function

default = {}
i = default
i.__setitem__('h',1) # Integration Step Size; 1 = full step
i.__setitem__('hfun',0.995) # Linspace approach fraction to solution; must be <= 1
i.__setitem__('delta_min',1e-30) # Zero Approximation, for stiff problems, the lower the better
i.__setitem__('criteria',1e-8) # Convergence Criteria
i.__setitem__('inner_criteria',i['criteria']/10) # Linear Problem ALgorithm Convergence Criteria
i.__setitem__('convtol',100) # Maximum Number of Convergence Attempts
i.__setitem__('convtolH',20) # Hessian Loop Maximum Iterations
i.__setitem__('inner_convtol',300) # QMR Maximum Step Numbers
i.__setitem__('max_restarts',100) # maximum number of coverages restarts
i.__setitem__('max_time',60) # maximum attempt period (in seconds)

convergence_params = default
