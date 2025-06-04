# mpirun -n 9 /home/yp/projects/miss_cat/miss_cat-opt -i grain_growth_2D_graintracker.i
[Mesh]
  type = GeneratedMesh
  dim = 2
  xmax = 200
  ymax = 200
  nx = 50
  ny = 50
[]

[Modules]
  [./PhaseField]
    [./Conserved]
      [./c]
        args = 'eta_1 eta_2 eta_3 eta_0'
        free_energy = f_chem
        kappa = kappa_c
        mobility = M
        solve_type = REVERSE_SPLIT
      [../]
    [../]
    [./Nonconserved]
      [./eta_1]
        args = 'c eta_2 eta_3 eta_0'
        free_energy = f_chem
        kappa = kappa_c
        mobility = L
      [../]
      [./eta_2]
        args = 'c eta_1 eta_3 eta_0'
        free_energy = f_chem
        kappa = kappa_c
        mobility = L
      [../]
      [./eta_3]
        args = 'c eta_1 eta_2 eta_0'
        free_energy = f_chem
        kappa = kappa_c
        mobility = L
      [../]
      [./eta_0]
        args = 'c eta_1 eta_2 eta_3'
        free_energy = f_chem
        kappa = kappa_c
        mobility = L
      [../]
    [../]
  [../]
[]

[ICs]
  [PolycrystalICs]
    [PolycrystalColoringIC]
      polycrystal_ic_uo = voronoi
      var_name_base = eta_
      op_num = 4
    []
  []
[]

[BCs]
  [./Periodic]
    [./c_bcs]
      auto_direction = 'x y'
    [../]
  [../]
[]

[AuxVariables]
  [./bnds]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[AuxKernels]
  [bnds_aux]
    # AuxKernel that calculates the GB term
    type = BndsCalcAux
    variable = bnds #计算晶界(Grain Boundary)位置，通过计算相邻晶粒之间的界面能量梯度。
    v = 'eta_1 eta_2 eta_3 eta_0'
    execute_on = 'initial timestep_end'
  []
[]

[Materials]
  [./constants]
    type = GenericConstantMaterial
    prop_names = 'M L kappa_c kappa_eta'
    prop_values = '1.0 1.0 1.0 1.0'
  [../]
  [./f_chem]           # Local free energy function (eV/mol)
    type = DerivativeParsedMaterial
    f_name = f_chem
    args = 'c eta_1 eta_2 eta_3 eta_0'
    constant_names = 'roh ca cb w alpha'
    constant_expressions = '1.41421356 0.3 0.7 1 5'
    function = 'roh*roh*(c-ca)^2*(1-eta_1^3*(6*eta_1^2-15*eta_1+10)-eta_2^3*(6*eta_2^2-15*eta_2+10)-eta_3^3*(6*eta_3^2-15*eta_3+10)-eta_0^3*(6*eta_0^2-15*eta_0+10))
                +roh*roh*(cb-c)^2*(eta_1^3*(6*eta_1^2-15*eta_1+10)+eta_2^3*(6*eta_2^2-15*eta_2+10)+eta_3^3*(6*eta_3^2-15*eta_3+10)+eta_0^3*(6*eta_0^2-15*eta_0+10))
                +w*((eta_1^2*(1-eta_1)^2+eta_2^2*(1-eta_2)^2+eta_3^2*(1-eta_3)^2+eta_0^2*(1-eta_0)^2)
                +alpha*2*(eta_1^2*eta_2^2+eta_1^2*eta_3^2+eta_1^2*eta_0^2+eta_2^2*eta_3^2+eta_2^2*eta_0^2+eta_3^2*eta_0^2))'
    derivative_order = 2
  [../]
  [./free_energy]
    # add the chemical and nucleation free energy contributions together
    type = DerivativeSumMaterial
    derivative_order = 2
    coupled_variables = 'c eta_1 eta_2 eta_3 eta_0'
    sum_materials = 'f_chem Fn'
  [../]
  [./probability]
    type = ParsedMaterial
    property_name = P
    coupled_variables = c
    expression = max(-10*c,0)
    outputs = exodus
  [../]
  [./nucleation]
    type = DiscreteNucleation
    property_name = Fn
    op_names  = c
    op_values = 0.90
    penalty = 5
    penalty_mode = MIN
    map = map
    outputs = exodus
  [../]
[]

[UserObjects]
  [./inserter]
    type = DiscreteNucleationInserter
    hold_time = 0.1
    probability = P
    radius = 3
  [../]
  [./map]
    type = DiscreteNucleationMap
    periodic = c
    inserter = inserter
  [../]
  [voronoi]
    type = PolycrystalVoronoi
    grain_num = 4 # Number of grains
    rand_seed = 2
    int_width = 5
    variable = 'eta_1 eta_2 eta_3 eta_0'
  []
  [grain_tracker]
    type = GrainTracker
    variable = 'eta_1 eta_2 eta_3 eta_0'
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  nl_rel_tol = 1e-08
  end_time = 1e+6
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
  exodus = true
  file_base = 'results4/1'
[]
