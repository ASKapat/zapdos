dom0Scale=1

[Mesh]
  [./file]
    type = FileMeshGenerator
    file = 'wall_refined.msh'
  [../]

  [./left]
    type = SideSetsFromNormalsGenerator
    normals = '-1 0 0'
    new_boundary = 'left'
    input = file
  [../]
  [./right]
    type = SideSetsFromNormalsGenerator
    normals = '1 0 0'
    new_boundary = 'right'
    input = left
  [../]
[]

[GlobalParams]
  potential_units = kV
  use_moles = true
[]

[Problem]
  type = FEProblem
  # kernel_coverage_check = false
[]

[Debug]
  show_var_residual_norms = true
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  automatic_scaling = true
  compute_scaling_once = false
  end_time = 10
  petsc_options = '-snes_converged_reason -snes_linesearch_monitor'
  solve_type = pjfnk
  line_search = 'basic'
  petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'lu NONZERO 1.e-10'
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-10
  l_tol = 1e-3
  l_max_its = 100
  nl_max_its = 25
  dtmin = 1e-14
  dtmax = 1
  steady_state_detection = true

  [./TimeStepper]
    type = IterationAdaptiveDT
    cutback_factor = 0.4
    dt = 1e-11
    growth_factor = 1.2
    optimal_iterations = 10
  [../]

[]

[Outputs]
  perf_graph = true
  [out_01]
    type = Exodus
  []
[]

[UserObjects]
#   [./data_provider]
#     type = ProvideMobility
#     e = 1.6e-19
#   [../]
[]

[Variables]
  [potential]
  []
  
  [E_D+]
    block = 0
  []
  
  [E_D]
    block = 0
  []
  
  [E_Li+]
    block = 0
  []

#   [E_Li]
#     block = 0
#   []
  
[]

[DriftDiffusionAction]
  [plasma]
    electrons = 'em'
    charged_particle = 'Li+ D+'
    Neutrals = 'Li D'
    potential = potential
    mean_energy = mean_en
    using_offset = true
    offset = 40
    use_ad = true
    position_units = ${dom0Scale}
    Additional_Outputs = 'ElectronTemperature'
    block = 0
  []
[]

[AuxVariables] 
  [T_D+]
    block = 0
    order = CONSTANT
    family = MONOMIAL
  []
  
  [T_D]
    block = 0
    order = CONSTANT
    family = MONOMIAL
  []
  
  [T_Li+]
    block = 0
    order = CONSTANT
    family = MONOMIAL
  []
#   
#   [T_Li]
#     block = 0
#     order = CONSTANT
#     family = MONOMIAL
#   []  
[]

[Kernels]
# ########################D+ energy terms######################
  [E_D+_time_deriv]
    type = ElectronTimeDerivative
    variable = E_D+
    block = 0
  []
  
  [E_D+_Efield_adv]
    type = EFieldAdvection
    variable = E_D+
    potential = potential
    position_units = ${dom0Scale}
    block = 0
  []
  
  [E_D+_diff]
    type = CoeffDiffusion
    variable = E_D+
    position_units = ${dom0Scale}
    block = 0
  []
  
  [E_D+_ionize]
    type = ADEEDFEnergyLog_HS
    variable = E_D+
    electrons = em
    target = D
    threshold_energy = E_D
    reaction = 'D         + em  ->  D+     + em + em'
    gain_loss = 1.0
    number = 0
  []
  
  [E_D+_recomb]
    type = ADEEDFEnergyLog_HS
    variable = E_D+
    electrons = em
    target = D+
    threshold_energy = E_D+
    reaction = 'D+        + em  ->  D'
    gain_loss = -1.0
    number = 1
  []
# 
# ########################D energy terms######################
  [E_D_time_deriv]
    type = ElectronTimeDerivative
    variable = E_D
    block = 0
  []
  
  
  [E_D_diff]
    type = CoeffDiffusion
    variable = E_D
    position_units = ${dom0Scale}
    block = 0
  []
  
  [E_D_ionize]
    type = ADEEDFEnergyLog_HS
    variable = E_D
    electrons = em
    target = D
    threshold_energy = E_D
    reaction = 'D         + em  ->  D+     + em + em'
    gain_loss = -1.0
    number = 0
  []
  
  [E_D_recomb]
    type = ADEEDFEnergyLog_HS
    variable = E_D
    electrons = em
    target = D+
    threshold_energy = E_D+
    reaction = 'D+        + em  ->  D'
    gain_loss = 1.0
    number = 1
  []
  
########################Li+ energy terms######################
  [E_Li+_time_deriv]
    type = ElectronTimeDerivative
    variable = E_Li+
    block = 0
  []
#   
#   [E_Li+_Efield_adv]
#     type = EFieldAdvection
#     variable = E_Li+
#     potential = potential
#     position_units = ${dom0Scale}
#     block = 0
#   []
#   
#   [E_Li+_diff]
#     type = CoeffDiffusion
#     variable = E_Li+
#     position_units = ${dom0Scale}
#     block = 0
#   []
  
#   [E_Li+_ionize]
#     type = ADEEDFEnergyLog_HS
#     variable = E_Li+
#     electrons = em
#     target = Li
#     threshold_energy = E_Li
#     reaction = 'Li        + em  ->  Li+    + em + em'
#     gain_loss = 1.0
#     number = 2
#   []
#   
#   [E_Li+_recomb]
#     type = ADEEDFEnergyLog_HS
#     variable = E_Li+
#     electrons = em
#     target = Li+
#     threshold_energy = E_Li+
#     reaction = 'Li+       + em  ->  Li'
#     gain_loss = -1.0
#     number = 3
#   []


########################Li energy terms######################
#   [E_Li_time_deriv]
#     type = ElectronTimeDerivative
#     variable = E_Li
#     block = 0
#   []
#   
#   
#   [E_Li_diff]
#     type = CoeffDiffusion
#     variable = E_Li
#     position_units = ${dom0Scale}
#     block = 0
#   []
  
#   [E_Li_ionize]
#     type = ADEEDFEnergyLog_HS
#     variable = E_Li
#     electrons = em
#     target = Li
#     threshold_energy = E_Li
#     reaction = 'Li        + em  ->  Li+    + em + em'
#     gain_loss = -1.0
#     number = 2
#   []
#   
#   [E_Li_recomb]
#     type = ADEEDFEnergyLog_HS
#     variable = E_Li
#     electrons = em
#     target = Li+
#     threshold_energy = E_Li+
#     reaction = 'Li+       + em  ->  Li'
#     gain_loss = 1.0
#     number = 3
#   []
[]
########################Aux Kernels ######################
[AuxKernels]
	[D+_temp]
		type = ElectronTemperature
		variable = T_D+
		electron_density = D+
		mean_en = E_D+
		block = 0
	[]

	[D_temp]
		type = ElectronTemperature
		variable = T_D
		electron_density = D
		mean_en = E_D
		block = 0
	[]

	[T_Li+]
		type = ElectronTemperature
		variable = T_Li+
		electron_density = Li+
		mean_en = E_Li+
	[]
# 
# 	[T_Li]
# 		type = ElectronTemperature
# 		variable = T_Li
# 		electron_density = Li
# 		mean_en = E_Li
# 	[]
[]


[BCs]
########################################### Boundary conditions##########################
############################################potential BCs################################

#Current methodology: Far enough upstream so that E-field is zero, zero voltage for now

  [./potential_left]
    type = DirichletBC
    variable = potential 
    boundary = 'left'
    value = 0.0
  [../]
  
  [./potential_right]
    type = NeumannBC
    variable = potential
    boundary = 'right'
    value = 0.0
  [../]
  
  
############################################Electron BCs##################################
 
 #Current methodology: Constant flux (Neumann) upstream, abosorbing (Hagelaar) at wall
  

  # Setting no BC on the left means the code automatically
  # applies a NeumannBC with zero flux. 

  [./em_physical_left]
    type = HagelaarElectronBC
    variable = em
    boundary = 'left'
    potential = potential
    mean_en = mean_en
    r = 0.0
    position_units = ${dom0Scale}
    
  [../]

  [./em_physical_right]
    type = NeumannBC
    variable = em
    boundary = 'right'
    potential = potential
    value = 0.065094652896
    position_units = ${dom0Scale}
  []

	
	

##############################################Electron energy############################

#Current methodology: 20eV Dirichlet upstream, absorbing (Hagelaar) at wall 


  
  [mean_en_physical_right]
    type = ElectronTemperatureDirichletBC
    variable = mean_en
    boundary = 'right'
    potential = potential
    em = em
    value = 10
    position_units = ${dom0Scale}
  []

  # Electron energy should be absorbed by the wall. 
  [mean_en_physical_left]
    type = HagelaarEnergyBC
    variable = mean_en
    r = 0.0
    boundary = 'left'
    potential = potential
    em = em
    position_units = ${dom0Scale}
  []
  

  
##############################D+ BCs######################################################



  
  [D+_physical_left_adv]
    type = HagelaarIonAdvectionBC
    variable = D+
    boundary = 'left'
    potential = potential
    r = 0.0
    position_units = ${dom0Scale}
  []
  [D+_physical_left_diff]
    type = HagelaarIonDiffusionBC_Ti_Te
    variable = D+
    boundary = 'left'
    r = 0.0
    mean_en = E_D+
    position_units = ${dom0Scale}
  []
  
  [./D+_physical_right]
    type = NeumannBC
    variable = D+
    boundary = 'right'
    potential = potential
    value = 0.065094652896
    position_units = ${dom0Scale}
  []
  
  
##############################################D+ energy############################

#Current methodology: 10eV Dirichlet upstream, absorbing (Hagelaar) at wall 


  
  [E_D+_physical_right]
    type = ElectronTemperatureDirichletBC
    variable = E_D+
    boundary = 'right'
    potential = potential
    em = D+
    value = 10
    position_units = ${dom0Scale}
  []

  # Electron energy should be absorbed by the wall. 
  [E_D+_physical_left]
    type = HagelaarEnergyBC
    variable = E_D+
    boundary = 'left'
    potential = potential
    em = D+
    r = 0.0
    position_units = ${dom0Scale}
  []  
  
  
  
  
##################################D BCs###################################################

  

    
  [D_physical_left_diff]
    type = HagelaarIonDiffusionBC_Ti_Te
    variable = D
    boundary = 'left'
    r = 0.0
    mean_en = E_D
    position_units = ${dom0Scale}
  []
  
  [D_physical_right]
    type = NeumannBC
    variable = D
    boundary = 'right'
    value = -0.006
  []
##############################################D energy############################

#Current methodology: 20eV Dirichlet upstream, absorbing (Hagelaar) at wall 


  
  [E_D_physical_right]
    type = ElectronTemperatureDirichletBC
    variable = E_D
    boundary = 'right'
    potential = potential
    em = D
    value = 10
    position_units = ${dom0Scale}
  []

  # Electron energy should be absorbed by the wall. 
  [E_D_physical_left]
    type = HagelaarEnergyBC
    variable = mean_en
    boundary = 'left'
    potential = potential
    em = D
    r = 0.0
    position_units = ${dom0Scale}
  []  
  
  
##############################Li+ BCs#####################################################

  # Setting absorbing BCs for Li+ on both left and right boundaries

 #  [Li+_physical_left_adv]
#     type = HagelaarIonAdvectionBC
#     variable = Li+
#     boundary = 'left'
#     potential = potential
#     r = 0.0
#     position_units = ${dom0Scale}
#   []
#   [Li+_physical_left_diff]
#     type = HagelaarIonDiffusionBC_Ti_Te
#     variable = Li+
#     boundary = 'left'
#     r = 0.0
#     mean_en = mean_en
#     position_units = ${dom0Scale}
#   []

  [Li+_physical_left]
    type = DirichletBC
    variable = Li+
    boundary = 'left'
    value =  -20 #2.734E18
  []
  
#   [Li+_physical_left]
#     type = NeumannBC
#     variable = Li+
#     boundary = 'left'
#     potential = potential
#     value = 0.065094652896
#     position_units = ${dom0Scale}
#   []
  



  [Li+_physical_right_adv]
    type = HagelaarIonAdvectionBC
    variable = Li+
    boundary = 'right'
    potential = potential
    r = 0.0
    position_units = ${dom0Scale}
  []
  [Li+_physical_right_diff]
    type = HagelaarIonDiffusionBC_Ti_Te
    variable = Li+
    boundary = 'right'
    r = 0.0
    mean_en = mean_en
    position_units = ${dom0Scale}
  []


# ##############################################Li+ energy############################

  
  
  # [E_Li+_physical_left]
#     type = HagelaarEnergyBC
#     variable = E_Li+
#     boundary = 'left'
#     potential = potential
#     em = Li+
#     r = 0.0
#     position_units = ${dom0Scale}
#   []  
  
  [E_Li+_physical_left]
    type = ElectronTemperatureDirichletBC
    variable = E_Li+
    boundary = 'left'
    potential = potential
    em = Li+
    value = 1
    position_units = ${dom0Scale}
  []
# 
# 
# 
  [E_Li+_physical_right]
    type = HagelaarEnergyBC
    variable = E_Li+
    boundary = 'right'
    potential = potential
    em = Li+
    r = 0.0
    position_units = ${dom0Scale}
  []  
  

######################################Li BCs##############################################



  [Li_physical_left]
    type = DirichletBC
    variable = Li
    boundary = 'left'
    value =  -10 #2.734E19
  []
  
#   [Li_physical_left]
#     type = NeumannBC
#     variable = Li
#     boundary = 'left'
#     value = -0.006
#   []



  [Li_physical_right_diff]
    type = HagelaarIonDiffusionBC_Ti_Te
    variable = Li
    boundary = 'right'
    r = 0.0
    mean_en = mean_en
    position_units = ${dom0Scale}
  []
 


# ##############################################Li energy############################

  # Electron energy should be absorbed by the wall. 
  

  # Li energy is evap temp. 
#   [E_Li_physical_left]
#     type = ElectronTemperatureDirichletBC
#     variable = E_Li
#     boundary = 'left'
#     potential = potential
#     em = Li
#     value = 10
#     position_units = ${dom0Scale}
#   []
#   
#   [E_Li_physical_right]
#     type = HagelaarEnergyBC
#     variable = E_Li
#     boundary = 'right'
#     potential = potential
#     em = Li
#     r = 0.0
#     position_units = ${dom0Scale}
#   []  
[]
########################################Initial conditions############################

[ICs]
  [em_ic]
    type = ConstantIC
    variable = em
    value = -7.07393416931 #n_e = 5.1E20
    block = 0
  []
  [mean_en_ic]
    type = ConstantIC
    variable = mean_en
    value = -5.301977
    block = 0
  []
  
  [D+_ic]
    type = ConstantIC
    variable = D+
    value = -7.07393416931 #n_D+ = 5.1E20
    block = 0
  []
  [E_D+_ic]
    type = ConstantIC
    variable = E_D+
    value = -5.301977
    block = 0
  []
  
  [D_ic]
    type = ConstantIC
    variable = D
    value = -6
    block = 0
  []
  [E_D_ic]
    type = ConstantIC
    variable = E_D
    value = -5.301977
    block = 0
  []
  
  [Li+_ic]
    type = ConstantIC
    variable = Li+
    value = -20
    block = 0
  []
  [E_Li+_ic]
    type = ConstantIC
    variable = E_Li+
    value = -19.59
    block = 0
  []
  
#   [Li_ic]
#     type = ConstantIC
#     variable = Li
#     value = -6 #n_Li+ = 2.734E18
#     block = 0
#   []
#   [E_Li_ic]
#     type = ConstantIC
#     variable = E_Li
#     value = -19.49
#     block = 0
#   []
  [Li_ic]
    type = FunctionIC
    variable = Li
    function = 'lithium_ic'
    block = 0
  []
#   [E_Li_ic]
#     type = FunctionIC
#     variable = E_Li
#     function = 'lithiumE_ic'
#     block = 0
#   []
[]

###################################Functions#######################################
[Functions]

  [lithium_ic]
    type = ParsedFunction
    value = '-100*x - 10'
  []
#   
#   [lithiumE_ic]
#     type = ParsedFunction
#     value = '-100*x - 10.9' #lithium_ic times T_Li = 0.4
#   []
[]

###################################Materials#######################################
[Materials]
  [se_left]
    type = GenericConstantMaterial
    boundary = 'left'
    prop_names = 'se_coeff'
    prop_values = 0.0
  []
  [se_right]
    type = GenericConstantMaterial
    boundary = 'right'
    prop_names = 'se_coeff'
    prop_values = 0.0
  []
  
  [electron_moments]
    type = ADGasElectronMoments
    block = 0
    em = em
    mean_en = mean_en
    property_tables_file = 'e_mob_test.txt'
    interp_trans_coeffs = true
    ramp_trans_coeffs = false
    interp_elastic_coeff = false
    user_p_gas = 0.133
  []

  [gas_constants]
    type = GenericConstantMaterial
    block = 0
    prop_names  =  'e         N_A       k_boltz   eps         se_coeff    se_energy   T_gas     massem   p_gas'
    prop_values =  '1.6e-19   6.022e23  1.38e-23   8.854e-12   0.0         3.          766.094   9.11e-31 0.133'
  []

  [gas_phase]
    type = ADGenericConstantMaterial
    prop_names = 'diffpotential'
    prop_values = '8.85e-12'
  []
  
  [Deuterium_ion]
    type = ADHeavySpecies_E
    heavy_species_name = D+
    heavy_species_mass = 3.34358e-27
    heavy_species_charge = 1.0
    heavy_species_energy = E_D+
    heavy_species_density = D+
    diffusivity = 0.1
    block = 0
  []
  [Deuterium]
    type = ADHeavySpecies_E
    heavy_species_name = D
    heavy_species_mass = 3.343458e-27
    heavy_species_charge = 0.0
    heavy_species_energy = E_D
    heavy_species_density = D
    diffusivity = 0.1 
    block = 0
  []
  
  [lithium_ion]
    type = ADHeavySpecies_E
    heavy_species_name = Li+
    heavy_species_mass = 1.15258e-26
    heavy_species_charge = 1.0
    heavy_species_energy = E_Li+
    heavy_species_density = Li+
    diffusivity = 0.1
    block = 0
  []
#   [lithium]
#     type = ADHeavySpecies_E
#     heavy_species_name = Li
#     heavy_species_mass = 1.15258e-26
#     heavy_species_charge = 0.0
#     heavy_species_energy = E_Li
#     heavy_species_density = Li
#     diffusivity = 0.1
#     block = 0
#   []
#   [Deuterium_ion]
#     type = ADHeavySpecies
#     heavy_species_name = D+
#     heavy_species_mass = 3.34358e-27
#     heavy_species_charge = 1.0
#     diffusivity = 0.1
#     block = 0
#   []
#   [Deuterium]
#     type = ADHeavySpecies
#     heavy_species_name = D
#     heavy_species_mass = 3.343458e-27
#     heavy_species_charge = 0.0
#     diffusivity = 0.1 
#     block = 0
#   []
#   
#   [lithium_ion]
#     type = ADHeavySpecies
#     heavy_species_name = Li+
#     heavy_species_mass = 1.15258e-26
#     heavy_species_charge = 1.0
#     diffusivity = 0.1
#     block = 0
#   []
  [lithium]
    type = ADHeavySpecies
    heavy_species_name = Li
    heavy_species_mass = 1.15258e-26
    heavy_species_charge = 0.0
    diffusivity = 0.1
    block = 0
  []    
[]

[Reactions]
  active = 'all'
  [all]
    species = 'em Li+ D+ Li D'
    reaction_coefficient_format = 'rate'
    file_location = 'RATES_Zapdos'
    electron_energy = 'mean_en'
    electron_density = 'em'
    include_electrons = true
    potential = 'potential'
    use_log = true
    use_ad = true
    position_units = ${dom0Scale}
    # These are parameters required equation-based rate coefficients
    block = 0
    reactions = 'D         + em  ->  D+     + em + em : EEDF [-15.4666]   (scd12_h_IPRT=_1_IGRD=_1_Z1=_1.dat)      #0
                 D+        + em  ->  D                : EEDF              (acd12_h_IPRT=_1_IGRD=_1_Z1=_1.dat)      #1
                 Li        + em  ->  Li+    + em + em : EEDF [-5.39]      (scd96r_li_IPRT=_1_IGRD=_1_Z1=_1.dat)    #2
                 Li+       + em  ->  Li               : EEDF              (acd96r_li_IPRT=_1_IGRD=_1_Z1=_1.dat)'   #3
                  

  []


[]