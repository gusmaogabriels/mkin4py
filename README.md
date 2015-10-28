============
**mkin4py** 
============
*mkin(microkinetics) 4 py(thon)*

A general package for linearly defining and solvig microkinetic catalytic systems.

================
**Description**
================

  A microkinetic package translated from Linearized Microkinetic Catalytic System Solver

  By the author:

  *Gabriel S. Gusmão* <gusmaogabriels@gmail.com>
  under *Dr. Phillip Christopher* <christopher@engr.ucredu> advisement.

  For detailed information, refer to code comments or associated publication.
  Gusmão, G. S. & Christopher, P., *A general and robust approach for defining and solving 
  microkinetic catalytic systems.* AIChE J. 00, (2014).; http://dx.doi.org/10.1002/aic.14627

  The 17-Step Ethylene Epoxidation by Stegelmann et al. has been used as example.
  Stegelmann, C., Schiødt, N. C., Campbell, C. T. & Stoltze, P.
  *Microkinetic modeling of ethylene oxidation over silver*. J. Catal. 221, 630–649 (2004).

  1. Set-up the environment conditions (temperature, pressure, gas constant)
  2. Create a MK (microkinetic model object)
     - Define its dimensions: number of reactants (rows) and elementary reactions (columns) involved in the stoichiometry matrix, and parse the rows that refer to *free*-species (non-adsorbed)
     - Parse the stoichiometry matrix (must be of size number of reactants × number of elementary reactions)
     - Set the kinetic parameters: Activation Energies and Pre-exponential factors (must be of the size of the involved elementary reactions)
     - Set the fixed concentration of *free*-species (molar fraction in non-adsorbed phase)
     - Parse the string-labels of involved species (array of size of number of species) 
  3. Solve the decurrent ensuing LP (linear problem)
     - For now, there is only available a *Newton*-type method.
     - Standard iterative-procedure adopted for solving the inner-loop LP (Quasi-minimum residue)

  The convergence parameters are set as default in the module `solver` in `.params`

================
**Features**
================

   **Linearization**

  The project makes use of explicit routines for the calculation of the MK model derivatives
  
    - *Jacobian*: Available as standard.
    - *Hessian*: Used in the convex two-step method (details in the aforementioned reference)

================
**On the way**
================

  1. Additional LP solvers in "switchable" fashion.
  2. Evolutionary methods for the definition of best convergence parameters for *stiff* problems (when TOF`s are close to the machine precision)

================
**Instructions**
================

  - **Installation**

        pip install mk4py==version_no

  - **Example**: Stoltze's 17-Step Ethylene Epoxidation MK system

        import mkin4py
        import numpy as np
        
        # Environment Conditions
        T = 500; #K
        P = 2; #bar
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
        
  - **Evaluation**:
       
        sol = mkin4py.solver.solve.rk4() # 4th-order Runge-Kutta method coupled within the LP solved via QMR
        # Outupts
        print '...'
        print sol['msg'], 'time: ', sol['time']
        print 'Coverage'
        print sol['coverage']
        print 'Rates'
        print sol['rates']

  - **Output**:

        ...
        Convergence achieved time:  2.25999999046
        Coverage
        [[  5.00000000e-01]
        [  5.00000000e-01]
        [  0.00000000e+00]
        [  0.00000000e+00]
        [  0.00000000e+00]
        [  0.00000000e+00]
        [  4.39342950e-01]
        [  1.19743307e-03]
        [  1.07992516e-01]
        [  1.10447591e-01]
        [  2.99730332e-09]
        [  7.70711567e-10]
        [  1.00256049e-01]
        [  9.78269843e-02]
        [  1.32727193e-01]
        [  9.64368428e-03]
        [  3.28419897e-08]
        [  4.59118359e-13]
        [  5.65562897e-04]
        [  1.12150752e-15]]
        Rates
        [[ -4.24252190e+01]
        [ -2.49532633e+01]
        [  1.29732696e+01]
        [  5.58723011e-04]
        [  2.39588699e+01]
        [  2.39588699e+01]
        [  0.00000000e+00]
        [  3.65929509e-13]
        [  0.00000000e+00]
        [  4.32857086e-11]
        [  2.76796815e-16]
        [  2.87485591e-11]
        [  0.00000000e+00]
        [ -2.27373675e-13]
        [ -4.65661287e-10]
        [  9.86479981e-16]
        [  0.00000000e+00]
        [ -1.60491195e-13]
        [  4.65661287e-10]
        [ -1.42115222e-11]]
