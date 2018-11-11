# -*- coding: utf-8 -*-
from __future__ import division, absolute_import, print_function
from . import np

class Environment(object):
    """ Environment
    Class that encompasses environment variables that parametrizes the MK model
    """
    def __init__(self):
        self.gas_constant = []
        self.temperature = []
        self.pressure = []
        self.__status__ = False

    def set_gas_constant(self,x):
        if not isinstance(x,(int,long,float)):
            raise TypeError('Gas constant must be numeric (int,long,float)')
        else:
            self.gas_constant = x
            self.__validator()

    def set_temperature(self,x):
        if not isinstance(x,(int,long,float)):
            raise TypeError('Temperature must be numeric (int,long,float)')
        else:
            self.temperature = x
            self.__validator()

    def set_pressure(self,x):
        if not isinstance(x,(int,long,float)):
            raise TypeError('Pressure must be numeric (int,long,float)')
        else:
            self.pressure = x
            self.__validator()

    def __validator(self):
        if not any([x == [] for x in [self.gas_constant, self.temperature, self.pressure]]):
            self.__status__ = True
        else:
            pass

class MKmodel(object):

    def __init__(self,environment):
        self.environment = environment
        self.__mkinit_status__ = True
        self.dims = {}
        i = self.dims
        i.__setitem__('num_sp',[])
        i.__setitem__('num_reac',[])
        self.ms = []  # stoichsiometry Matrix
        self.kinetic_parameters = {}
        i = self.kinetic_parameters
        i.__setitem__('va', []) # Pre-exponential factors
        i.__setitem__('vea', []) # Activation barriers
        i.__setitem__('k', []) #  rate constants
        self.maps = {}
        i = self.maps
        i.__setitem__('stoichs',[])  # Non-adsorbed species (positions)
        i.__setitem__('msa',[])
        i.__setitem__('msas',[])
        i.__setitem__('surface',[]) # Adsorbed species position
        i.__setitem__('xsurface',[])
        i.__setitem__('ndof',[]) # Degree of freedom adsorbate (typically chosen randomly)
        self.splabels = []
        self.concs = []
        self.coverage = [] # initiate coverage
        self.__status__ = False
        self.__mkinit_status__ = False
        self.__param_status__ = {'ms':False,'kinetic_params':False,'concentrations':False,'splabels':False}

    def create(self,num_species, num_reactions, stoichs):
        if not isinstance(num_species, int):
            raise TypeError('Numnber of species (num_species) must be integer.')
        elif not isinstance(num_reactions, int):
            raise TypeError('Numnber of reactions (num_reactions) must be integer.')
        elif not (isinstance(stoichs,(tuple,list,np.ndarray)) and all([isinstance(x,(int)) for x in stoichs])):
            raise TypeError('stoichs must be a list, tuple or numpy.ndarray of integer values')
        else:
            pass
            """ Microkinetic model initialization
            params(2) :
                num_species as integer
                num_reactions as integer
            """
            self.dims['num_sp'] = num_species
            self.dims['num_reac'] = num_reactions
            self.maps['stoichs'] = np.array(stoichs).astype(int)
            self.maps['surface'] = np.setdiff1d(np.array(range(0,self.dims['num_sp'])),self.maps['stoichs']).astype(int)
            # Memory allocations
            self.ms = np.empty([self.dims['num_sp'],self.dims['num_reac']])  # stoichsiometry Matrix
            self.kinetic_parameters['va'] = np.empty([self.dims['num_reac'],1]) # Pre-exponential factors
            self.kinetic_parameters['vea'] = np.empty([self.dims['num_reac'],1]) # Activation barriers
            self.splabels = np.empty([self.dims['num_sp'],1])
            self.__mkinit_status__ = True

    def set_ms(self,ms):
        """ set_ms - Stochiometric Matrix parser
         params(1) : Stochiometric Matrix as numpy.ndarray, tuple or list
        """
        if not self.__mkinit_status__:
            raise Exception('MK model has not been initialized')
        elif not isinstance(ms,(np.ndarray,tuple,list)):
            raise TypeError('The stoichsiometric matrix must be parsed as numpy.ndarray, tuple or list.')
        elif np.array(ms).shape != (self.dims['num_sp'],self.dims['num_reac']):
            raise ValueError('Dimension mismatch. The stoichsiometry matrix must of dimensions {}-by{}'.format(self.num_sp,self.num_reac))
        else:
            self.ms = np.array(ms).astype(int)
            self.__param_status__['ms'] = True
            self.__validator()

    def set_kinetic_params(self,va,vea):
        """ set_kinetic_params - Kinetic parameters parser
         params(2) :
                 Pre-exponential factors as numpy.ndarray, tuple or list
                 Activation energies as numpy.ndarray, tupl or list
        """
        if not self.__mkinit_status__:
            raise Exception('MK model has not been initialized')
        elif not all([isinstance(x,(np.ndarray,tuple,list)) for x in [va,vea]]):
            raise TypeError('Kinetic parameters must be parsed as numpy.ndarray, tuple or list.')
        elif not all([np.array(x).shape == (self.dims['num_reac'],1) for x in [va,vea]]):
            raise ValueError('Dimension mismatch. Kinetic parameters must of dimensions {}-by-1'.format(self.dims['num_reac']))
        else:
            self.kinetic_parameters['va'] = np.array(va)
            self.kinetic_parameters['vea'] = np.array(vea)
            self.__param_status__['kinetic_params'] = True
            self.__validator()

    def set_concentrations(self,x):
        if not self.__mkinit_status__:
            raise Exception('MK model has not been initialized')
        elif not (isinstance(x,(tuple,list,np.ndarray)) and all([isinstance(i,(int,long,float)) for i in x])):
            raise TypeError('concs must be a list, tuple or numpy.ndarray of numeric values')
        else:
            self.concs = np.array(x)
            self.__param_status__['concentrations'] = True
            self.__validator()

    def set_splabels(self,splabels):
        """ set_splabels - Species labels parser
         params(1) : Species labels as numpy.ndarray, tuple or list
        """
        if not self.__mkinit_status__:
            raise Exception('MK model has not been initialized')
        elif not (isinstance(splabels,(tuple,list,np.ndarray)) and all([isinstance(x,(str)) for x in splabels])):
            raise TypeError('labels must be a list, tuple or numpy.ndarray of string values.')
        else:
            self.splabels = np.array(splabels)
            self.__param_status__['splabels'] = True
            self.__validator()

    def __validator(self):
        if all(self.__param_status__.values()):
            self.__status__ = True
            self.update_model()
            self.reset_model()
        else:
            pass

    def reset_model(self):
        if not self.environment.__status__:
            raise Exception('Environment conditions have not been initialized.')
        else:
            self.init_coverage()
            self.update_model()

    def update_model(self):
        if not self.environment.__status__:
            raise Exception('Environment conditions have not been initialized.')
        else:
            kpar = self.kinetic_parameters
            kpar['k'] = np.empty([self.dims['num_reac']] )
            for i in range(0,len(kpar['k'])):
                kpar['k'][i] = kpar['va'][i]*np.exp(-kpar['vea'][i]/(self.environment.gas_constant*self.environment.temperature))
            self.maps['msa'] = np.dot(self.ms,np.diag(kpar['k']))
            self.maps['msas'] = self.maps['msa'][self.maps['xsurface'],:]

    def init_coverage(self):
        if not self.environment.__status__:
            raise Exception('Environment conditions have not been initialized.')
        else:
            self.coverage = np.empty([self.dims['num_sp']])
            self.coverage[self.maps['stoichs']] = self.concs*self.environment.pressure
            self.coverage[self.maps['surface']] = np.random.rand(len(self.maps['surface']))
            self.coverage[self.maps['surface']] /= sum(self.coverage[self.maps['surface']])
            self.maps['ndof'] = self.__dof()
            self.maps['xsurface'] = self.maps['surface'][self.maps['surface']!=self.maps['ndof']]
            self.maps['msas'] = self.maps['msa'][self.maps['xsurface'],:]

    def __dof(self):
        return self.maps['surface'][np.random.randint(0,len(self.maps['surface']))]

environment = Environment()
mkmodel = MKmodel(environment)
