"""Environment and microkinetic model configuration, as in the original API."""
import math
from numbers import Integral, Real

import jax
import jax.numpy as np


def _positive(x, name, zero=False):
    if not isinstance(x, Real) or isinstance(x, bool):
        raise TypeError(f"{name} must be numeric")
    if not math.isfinite(x) or (x < 0 if zero else x <= 0):
        raise ValueError(f"{name} must be finite and {'nonnegative' if zero else 'positive'}")
    return float(x)


def _array(x, shape, name):
    value = np.asarray(x, dtype=float)
    if value.shape != shape:
        raise ValueError(f"{name} must have shape {shape}, got {value.shape}")
    if not bool(np.all(np.isfinite(value))):
        raise ValueError(f"{name} must contain finite numbers")
    return value


class Environment:
    """Temperature, pressure and gas constant (use consistent units)."""
    def __init__(self):
        self.gas_constant = self.temperature = self.pressure = None
        self.__status__ = False

    def set_gas_constant(self, x):
        self.gas_constant = _positive(x, "Gas constant")
        self.__validator()

    def set_temperature(self, x):
        self.temperature = _positive(x, "Temperature")
        self.__validator()

    def set_pressure(self, x):
        self.pressure = _positive(x, "Pressure", zero=True)
        self.__validator()

    def __validator(self):
        self.__status__ = all(x is not None for x in
                              (self.gas_constant, self.temperature, self.pressure))


class MKmodel:
    """Original stoichiometric model with one normalized surface-site balance.

    Configuration methods run on the host. Numerical arrays and solver kernels
    use JAX; each solve receives the current parameters explicitly.
    """
    def __init__(self, environment):
        self.environment = environment
        self.dims = {'num_sp': 0, 'num_reac': 0}
        self.ms = None
        self.kinetic_parameters = dict.fromkeys(('va', 'vea', 'k'))
        self.maps = dict.fromkeys(('stoichs', 'msa', 'msas', 'surface', 'xsurface', 'ndof'))
        self.splabels = []
        self.concs = self.coverage = None
        self.__status__ = self.__mkinit_status__ = False
        self.__param_status__ = dict.fromkeys(('ms', 'kinetic_params', 'concentrations', 'splabels'), False)
        self._key = jax.random.PRNGKey(0)

    def create(self, num_species, num_reactions, stoichs):
        for size in (num_species, num_reactions):
            if not isinstance(size, Integral) or isinstance(size, bool) or size < 1:
                raise ValueError("Species and reaction counts must be positive integers")
        indices = np.asarray(stoichs)
        if indices.ndim != 1 or (indices.size and not np.issubdtype(indices.dtype, np.integer)):
            raise TypeError("stoichs must be a one-dimensional sequence of integer indices")
        gas = [int(x) for x in indices.tolist()]
        if len(set(gas)) != len(gas) or any(x < 0 or x >= num_species for x in gas):
            raise ValueError("stoichs indices must be distinct and in range")
        surface = [i for i in range(num_species) if i not in gas]
        if not surface:
            raise ValueError("At least one surface species is required")
        self.dims = {'num_sp': int(num_species), 'num_reac': int(num_reactions)}
        self.maps['stoichs'] = np.asarray(gas, dtype=int)
        self.maps['surface'] = np.asarray(surface, dtype=int)
        self.maps['ndof'] = np.asarray(surface[0], dtype=int)
        self.maps['xsurface'] = np.asarray(surface[1:], dtype=int)
        self.__param_status__ = dict.fromkeys(self.__param_status__, False)
        self.__mkinit_status__ = True
        self.__status__ = False
        self.coverage = None

    def _require_created(self):
        if not self.__mkinit_status__:
            raise ValueError("MK model has not been initialized")

    def set_ms(self, ms):
        self._require_created()
        value = _array(ms, (self.dims['num_sp'], self.dims['num_reac']), 'Stoichiometric matrix')
        if not bool(np.all(value == np.floor(value))):
            raise ValueError("Stoichiometric coefficients must be integers")
        self.ms = value.astype(int)
        self.__param_status__['ms'] = True
        self.__validator()

    def set_kinetic_params(self, va, vea):
        self._require_created()
        shape = (self.dims['num_reac'], 1)
        va, vea = _array(va, shape, 'Pre-exponential factors'), _array(vea, shape, 'Activation energies')
        if bool(np.any(va < 0)):
            raise ValueError("Pre-exponential factors must be nonnegative")
        self.kinetic_parameters.update(va=va, vea=vea)
        self.__param_status__['kinetic_params'] = True
        self.__validator()

    def set_concentrations(self, x):
        self._require_created()
        value = _array(x, (len(self.maps['stoichs']),), 'Concentrations')
        if bool(np.any(value < 0)):
            raise ValueError("Concentrations must be nonnegative")
        self.concs = value
        self.__param_status__['concentrations'] = True
        self.__validator()

    def set_splabels(self, splabels):
        self._require_created()
        labels = list(splabels)
        if len(labels) != self.dims['num_sp'] or not all(isinstance(x, str) for x in labels):
            raise ValueError("Provide one string label per species")
        self.splabels = labels
        self.__param_status__['splabels'] = True
        self.__validator()

    def __validator(self):
        self.__status__ = all(self.__param_status__.values())
        if self.__status__:
            self.reset_model()

    def reset_model(self, seed=None):
        self.update_model()
        self.init_coverage(seed)

    def update_model(self):
        if not self.__status__ or not self.environment.__status__:
            raise ValueError("Initialize all model and environment parameters first")
        kpar = self.kinetic_parameters
        kpar['k'] = kpar['va'].ravel() * np.exp(-kpar['vea'].ravel() /
                    (self.environment.gas_constant * self.environment.temperature))
        if not bool(np.all(np.isfinite(kpar['k']))):
            raise ValueError("Rate constants overflow at the configured temperature")
        self.maps['msa'] = self.ms * kpar['k'][None, :]
        self.maps['msas'] = self.maps['msa'][self.maps['xsurface']]
        if self.coverage is not None:
            self.coverage = self.coverage.at[self.maps['stoichs']].set(self.concs * self.environment.pressure)

    def init_coverage(self, seed=None):
        if not self.__status__ or not self.environment.__status__:
            raise ValueError("Initialize all model and environment parameters first")
        if seed is not None:
            self._key = jax.random.PRNGKey(seed)
        self._key, key, dof_key = jax.random.split(self._key, 3)
        surface = self.maps['surface']
        values = jax.random.uniform(key, (len(surface),), dtype=self.concs.dtype, minval=0.01, maxval=1.)
        self.coverage = np.zeros(self.dims['num_sp'], dtype=self.concs.dtype)
        self.coverage = self.coverage.at[self.maps['stoichs']].set(self.concs * self.environment.pressure)
        self.coverage = self.coverage.at[surface].set(values / np.sum(values))
        dof = int(jax.random.randint(dof_key, (), 0, len(surface)))
        self.maps['ndof'] = surface[dof]
        self.maps['xsurface'] = np.concatenate((surface[:dof], surface[dof+1:]))
        self.maps['msas'] = self.maps['msa'][self.maps['xsurface']]
        return self.coverage


environment = Environment()
mkmodel = MKmodel(environment)
