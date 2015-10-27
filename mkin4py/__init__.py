# -*- coding: utf-8 -*-
from __future__ import division, absolute_import, print_function

__author__ = {'Gabriel S. Gusmao' : 'gusmaogabriels@gmail.com'}
__version__ = '1.0'

"""

    *mkin4py*

    A microkinetic package tranlated from Linearized Microkinetic Catalytic System Solver

    By the author:

    *Gabriel S. Gusmão* <gusmaogabriels@gmail.com>
    under *Dr. Phillip Christopher* <christopher@engr.ucredu> advisement.

    For detailed information, refer to code comments or associated publication.
    Gusmão, G. S. & Christopher, P., A general and robust approach for defining and solving
    microkinetic catalytic systems. AIChE J. 00, (2014).; http://dx.doi.org/10.1002/aic.14627

    The 17-Step Ethylene Epoxidation by Stegelmann et al. has been used as example.
    Stegelmann, C., Schiødt, N. C., Campbell, C. T. & Stoltze, P.
    Microkinetic modeling of ethylene oxidation over silver. J. Catal. 221, 630–649 (2004).

    ~~~~

    :copyright: (c) 2015 Gabriel S. Gusmão
    :license: MIT, see LICENSE for more details.

"""

import numpy as np
from scipy.sparse import linalg
from time import time
from copy import copy as dc
__status__ = False

from . import solver
from . import bases
environment = bases.environment
mkmodel = bases.mkmodel

__all__ = ['__author__','__version__','np','linalg','time','dc']
