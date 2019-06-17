# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.3'
#       jupytext_version: 1.0.2
#   kernelspec:
#     display_name: ML
#     language: python
#     name: ml
# ---

# +
import mkin4py
import numpy as np

# Environment Conditions
T = 500; #K
P = 1; #bar
gas_constant = 8.31456e-3 # Gas Constant - kJ/(mol×K)

# Set the environment conditions
mkin4py.environment.set_temperature(T)
mkin4py.environment.set_gas_constant(gas_constant)
mkin4py.environment. set_pressure(P)

# Stoichsiometric Matrix
ms = [
[-1, 1, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],\
[ 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],\
[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],\
[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],\
[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0],\
[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1],\
[-1, 1,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1,-1, 1, 1,-1, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1],\
[ 1,-1,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],\
[ 0, 0, 2,-2,-2, 2,-1, 1, 1,-1,-1, 1, 0, 0, 0, 0, 1,-1,-6, 6, 0, 0,-1, 1,-1, 1,-5, 5, 1,-1, 0, 0, 0, 0],\
[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4,-4, 0, 0, 0, 0, 1,-1, 3,-3,-2, 2, 0, 0, 0, 0],\
[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0,-1, 1],\
[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0,-1, 1, 0, 0],\
[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],\
[ 0, 0, 0, 0, 2,-2, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],\
[ 0, 0, 0, 0, 0, 0, 1,-1,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],\
[ 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0,-1, 1,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],\
[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],\
[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1,-1, 1,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],\
[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0],\
[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1,-1, 1, 0, 0, 0, 0, 0, 0],\
];

nreac = [0, 1] # Reactant Rows in mS
nprod = [2, 3, 4, 5] # Product Rows in mS
stoichs = np.concatenate((nreac,nprod)) # Reactants and Products are not under PSSA
mkin4py.mkmodel.create(np.shape(ms)[0],np.shape(ms)[1],stoichs) # Initialize the model
mkin4py.mkmodel.set_ms(ms) # Set the stoichiometry matrix

# Species labels 
splabels = ['O2','C2H4','C2H4O','CH3CHO','CO2','H2O','*','O2*','O*','OH*',\
'H2O*','CO2*','C2H4*','O∙O*','C2H4∙O*','CH2CH2O∙O*','C2H4O∙O*','CH3CHO∙O*',\
'CH2CHOH∙O*','CH2CHO∙O*']

mkin4py.mkmodel.set_splabels(splabels) # Set species labels

# Pre-exponential Factors of Eelementary Reactions (1/s)
va =[2.71e5, 1.1e12, 4.0e12, 8.0e14, 2.0e7, 1.3e15, 7.2e7,\
2.2e11, 9.0e14, 5.3e14, 1.95e8, 4.8e12, 1.13e13, 2.11e12,\
9.0e12, 4.5e10, 2.9e13, 2.6e9, 2.0e20, 5.3e13, 7.2e7, 2.2e11,\
4.0e11, 3.1e14, 2.6e13, 1.3e9, 1.0e20, 5.5e13, 1.4e10, 1.0e11,\
3.6e14, 1.0e8, 5.9e14, 1.4e9]

# Activation Barriers for Elementary ReactionS (kJ/mol)
vea = [5.7000, 47.3000, 75.0000, 157.5000, 20.0000, 96.9000, 0, 37.1000, 112.0000,\
183.3000, 0, 39.1000, 95.0000, 93.5000, 95.0000, 204.3000, 41.9000, 4.4000,\
11.0000, 791.6000, 0, 30.1000, 32.0000, 42.8000, 86.0000, 106.1000, 0, 906.6000,\
65.6000, 50.0000, 38.9000, 0, 46.6000, 0]

# Set the kinetic parameters
mkin4py.mkmodel.set_kinetic_params(np.array(va,ndmin=2).T,np.array(vea,ndmin=2).T)

y = [0.5, 0.5, 0, 0, 0, 0] # Reactants and Products Initial Fraction
mkin4py.mkmodel.set_concentrations(y) # Set the *free*-species concentrations
# -

sol = mkin4py.solver.solve.rk4() # 4th-order Runge-Kutta method coupled within the LP solved via QMR
# Outupts
print ('...')
print (sol['msg'], 'time: ', sol['time'])
print ('Coverage')
print (sol['coverage'])
print ('Rates')
print (sol['rates'])
