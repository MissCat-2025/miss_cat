# mpirun -n 9 /home/yp/projects/miss_cat/miss_cat-opt -i lagrangemult.i 

L = 5
M = 5
C_alpha = 0.3
C_beta = 0.7
kappa_c = 3
kappa_eta = 3
rho2 = 2
w = 1
alpha = 5


[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 100
  ny = 100
  xmax = 200
  ymax = 200
[]
[AuxVariables]
  [./bnds]
    order = FIRST
    family = LAGRANGE
  [../]
[]
[AuxKernels]
  # AuxKernel block, defining the equations used to calculate the auxvars
  [bnds_aux]
    # AuxKernel that calculates the GB term
    type = BndsCalcAux
    variable = bnds #计算晶界(Grain Boundary)位置，通过计算相邻晶粒之间的界面能量梯度。
    execute_on = 'initial timestep_end'
    v = 'eta0 eta1 eta2 eta3'
  []
[]
[Variables]
  [./w]
    order = FIRST
    family = LAGRANGE
  [../]
  [./c]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = FunctionIC
      function = cIC
    [../]
  [../]

  [./eta0]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = FunctionIC
      function = eta_0IC
    [../]
  [../]

  [./eta1]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = FunctionIC
      function = eta_1IC
    [../]
  [../]

  [./eta2]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = FunctionIC
      function = eta_2IC
    [../]
  [../]

  [./eta3]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = FunctionIC
      function = eta_3IC
    [../]
  [../]
[]

[Kernels]
  [./deta0dt]
    type = TimeDerivative
    variable = eta0
  [../] 
  [./ACBulk0]
    type = AllenCahn
    variable = eta0
    coupled_variables = 'c eta1 eta2 eta3'
    f_name = F
  [../]
  [./ACInterface0]
    type = ACInterface
    variable = eta0
    kappa_name = kappa_eta
  [../]

  [./deta1dt]
    type = TimeDerivative
    variable = eta1
  [../]
  [./ACBulk1]
    type = AllenCahn
    variable = eta1
    coupled_variables = 'c eta2 eta0 eta3'
    f_name = F
  [../]
  [./ACInterface1]
    type = ACInterface
    variable = eta1
    kappa_name = kappa_eta
  [../]

  [./deta2dt]
    type = TimeDerivative
    variable = eta2
  [../]
  [./ACBulk2]
    type = AllenCahn
    variable = eta2
    coupled_variables = 'c eta1 eta0 eta3'
    f_name = F
  [../]
  [./ACInterface2]
    type = ACInterface
    variable = eta2
    kappa_name = kappa_eta
  [../]

    [./deta3dt]
      type = TimeDerivative
      variable = eta3
    [../]
    [./ACBulk3]
      type = AllenCahn
      variable = eta3
      coupled_variables = 'c eta1 eta0 eta2'
      f_name = F
    [../]
    [./ACInterface3]
      type = ACInterface
      variable = eta3
      kappa_name = kappa_eta
    [../]

  [./c_res]
    type = SplitCHParsed
    variable = c
    f_name = F
    kappa_name = kappa_c
    w = w
    coupled_variables = 'eta0 eta1 eta2 eta3'
  [../]
  [./w_res]
    type = SplitCHWRes
    variable = w
    mob_name = M
  [../]
  [./time1]
    type = CoupledTimeDerivative
    variable = w
    v = c
  [../]
[]

[BCs]
  [./Periodic]
    [./c_bcs]
      auto_direction = 'x y'
    [../]
  [../]
[]

[Materials]
  [Constants]
    type = GenericConstantMaterial
    prop_names = 'L M kappa_c kappa_eta'
    prop_values = '${L} ${M} ${kappa_c} ${kappa_eta}'
  []

  [./H]
    type = DerivativeParsedMaterial
    property_name = h
    coupled_variables = 'eta0 eta1 eta2 eta3'
    expression = 'eta0^3*(6*eta0^2-15*eta0+10)+eta1^3*(6*eta1^2-15*eta1+10)+eta2^3*(6*eta2^2-15*eta2+10)+eta3^3*(6*eta3^2-15*eta3+10)'

  [../]

  [./f_alpha_beta]
    type = DerivativeParsedMaterial
    property_name = f_alpha_beta
    coupled_variables = 'c eta0 eta1 eta2 eta3'
    constant_names = 'rho2 c_alpha c_beta'
    constant_expressions = '${rho2} ${C_alpha} ${C_beta}'
    material_property_names = 'h(eta0,eta1,eta2,eta3)'
    expression = 'rho2*(c-c_alpha)^2*(1-h)+rho2*(c_beta-c)^2*h'
    derivative_order = 2
  [../]
    [./f1]
      type = DerivativeParsedMaterial
      property_name = f1
      coupled_variables = 'eta0 eta1 eta2 eta3'
      constant_names = 'w'
      constant_expressions = '${w}'
      derivative_order = 2
      expression = 'w*(eta0^2*(1-eta0)^2 + eta1^2*(1-eta1)^2 + eta2^2*(1-eta2)^2 + eta3^2*(1-eta3)^2)'
    [../]
    [./f2]
      type = DerivativeParsedMaterial
      property_name = f2
      derivative_order = 2
      coupled_variables = 'eta0 eta1 eta2 eta3'
      constant_names = 'alpha w'
      constant_expressions = '${alpha} ${w}'
      expression = 'w*(alpha * (eta0^2*eta1^2 + eta0^2*eta2^2 + eta0^2*eta3^2 + eta1^2*eta0^2 + eta1^2*eta2^2 + eta1^2*eta3^2 + eta2^2*eta0^2 + eta2^2*eta1^2 + eta2^2*eta3^2 + eta3^2*eta0^2 + eta3^2*eta1^2 + eta3^2*eta2^2))'
    [../]

  [./F]
    type = DerivativeParsedMaterial
    property_name = F
    derivative_order = 2
    outputs = exodus
    coupled_variables = 'c eta0 eta1 eta2 eta3'
    material_property_names = 'f_alpha_beta(c,eta0,eta1,eta2,eta3) f1(eta0,eta1,eta2,eta3) f2(eta0,eta1,eta2,eta3)'
    expression = 'f_alpha_beta + f1 + f2'
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  #nl_abs_tol = 1e-11
  nl_rel_tol = 1e-08
  end_time = 1e+6
  #petsc_options_iname = '-pc_type -ksp_type'
  #petsc_options_value = 'bjacobi  gmres'
  # petsc_options_iname = '-pc_type -pc_hypre_type'
  # petsc_options_value = 'hypre  boomeramg'
   petsc_options_iname = '-pc_type -sub_pc_type'
   petsc_options_value = 'asm lu'    
  dt = 0.01
  [./Adaptivity]
    max_h_level = 1
    initial_adaptivity = 1
    refine_fraction = 0.9
    coarsen_fraction = 0.1
  [../]
[]

[Outputs]
  execute_on = 'timestep_end'
  exodus = true
[]
[Postprocessors]
  [./total_energy]          # Total free energy at each timestep
    type = ElementIntegralVariablePostprocessor
    variable = F
    execute_on = 'initial timestep_end'
  [../]
  []
[Functions]
  [./eta_0IC]
    type = ParsedFunction
    vars = 'epsl_eta fi i'
    vals = '0.1 1.5 1'
    value = 'epsl_eta*(cos((0.01*i)*x-4)*cos((0.007+0.01*i)*y)+cos((0.11+0.01*i)*x)*cos((0.11+0.01*i)*y)+fi*(cos((0.046+0.001*i)*x+(0.0405+0.001*i)*y)*cos((0.031+0.001*i)*x-(0.004+0.001*i)*y))^2)^2'
  [../]
  [./eta_1IC]
    type = ParsedFunction
    vars = 'epsl_eta fi i'
    vals = '0.1 1.5 2'
    value = 'epsl_eta*(cos((0.01*i)*x-4)*cos((0.007+0.01*i)*y)+cos((0.11+0.01*i)*x)*cos((0.11+0.01*i)*y)+fi*(cos((0.046+0.001*i)*x+(0.0405+0.001*i)*y)*cos((0.031+0.001*i)*x-(0.004+0.001*i)*y))^2)^2'
  [../]
  [./eta_2IC]
    type = ParsedFunction
    vars = 'epsl_eta fi i'
    vals = '0.1 1.5 3'
    value = 'epsl_eta*(cos((0.01*i)*x-4)*cos((0.007+0.01*i)*y)+cos((0.11+0.01*i)*x)*cos((0.11+0.01*i)*y)+fi*(cos((0.046+0.001*i)*x+(0.0405+0.001*i)*y)*cos((0.031+0.001*i)*x-(0.004+0.001*i)*y))^2)^2'
  [../]
  [./eta_3IC]
    type = ParsedFunction
    vars = 'epsl_eta fi i'
    vals = '0.1 1.5 4'
    value = 'epsl_eta*(cos((0.01*i)*x-4)*cos((0.007+0.01*i)*y)+cos((0.11+0.01*i)*x)*cos((0.11+0.01*i)*y)+fi*(cos((0.046+0.001*i)*x+(0.0405+0.001*i)*y)*cos((0.031+0.001*i)*x-(0.004+0.001*i)*y))^2)^2'
  [../]
  [./cIC]
    type = ParsedFunction
    vars = 'c0 epsl'
    vals = '0.5 0.01'
    value = 'c0+epsl*(cos(0.105*x)*cos(0.11*y)+(cos(0.13*x)*cos(0.087*y))^2+cos(0.025*x-0.15*y)*cos(0.07*x-0.02*y))'
  [../]
[]
