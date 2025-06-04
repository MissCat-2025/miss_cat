#
# Test the DiscreteNucleation material in a toy system. The global
# concentration is above the solubility limit, but below the spinodal.
# Without further intervention no nucleation will occur in a phase
# field model. The DiscreteNucleation material will locally modify the
# free energy to coerce nuclei to grow.
#
[GlobalParams]
  op_num = 3
  var_name_base = gr
  grain_num = 3
[]
[Variables]
  [./PolycrystalVariables]
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
    v = 'gr0 gr1 gr2 eta1'
    execute_on = 'initial timestep_end'
  []
[]
[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 120
  ny = 120
  xmax = 500
  ymax = 500
  elem_type = QUAD
[]

[Modules]
  [./PhaseField]
    [./Conserved]
      [./c]
        args = 'eta1 gr0 gr1 gr2'
        free_energy = F
        mobility = M
        kappa = kappa_c
        solve_type = REVERSE_SPLIT
      [../]
    [../]
    [./Nonconserved]
      [./eta1]
        args = 'c gr0 gr1 gr2'
        free_energy = F
        kappa = kappa_eta
        mobility = L
      [../]
    []
  [../]
[]

[ICs]
  [./c_IC]
    type = RandomIC
    variable = c
    min = 0.001
    max = 0.01
  [../]
  [PolycrystalICs]
    [PolycrystalColoringIC]
      polycrystal_ic_uo = voronoi
      var_name_base = gr
      # inactive = 
    []
  []
  [./eta1_IC]
    type = RandomIC
    variable = eta1
    min = 0.01
    max = 0.1
  [../]
[]
[Kernels]
  # 添加以下NullKernel配置阻止晶粒变量演化
  [gr0_null]
    type = NullKernel
    variable = gr0
  []
  [gr1_null]
    type = NullKernel
    variable = gr1
  []
  [gr2_null]
    type = NullKernel
    variable = gr2
  []
[]
[Materials]
  [./pfmobility]
    type = GenericConstantMaterial
    prop_names = 'M L kappa_c kappa_eta'
    prop_values = '1.0 1.0 1.0 1.0'
  [../]
    [./f_chem]           # Local free energy function (eV/mol)
    type = DerivativeParsedMaterial
    f_name = f_chem
    args = 'c eta1 gr0 gr1 gr2'
    constant_names = 'roh ca cb w alpha'
    constant_expressions = '1.41421356 0.3 0.7 1 5'
    function = 'roh*roh*(c-ca)^2*(1-eta1^3*(6*eta1^2-15*eta1+10)-gr0^3*(6*gr0^2-15*gr0+10)-gr1^3*(6*gr1^2-15*gr1+10)-gr2^3*(6*gr2^2-15*gr2+10))
                +roh*roh*(cb-c)^2*(eta1^3*(6*eta1^2-15*eta1+10)+gr0^3*(6*gr0^2-15*gr0+10)+gr1^3*(6*gr1^2-15*gr1+10)+gr2^3*(6*gr2^2-15*gr2+10))
                +w*((eta1^2*(1-eta1)^2+gr0^2*(1-gr0)^2+gr1^2*(1-gr1)^2+gr2^2*(1-gr2)^2)
                +alpha*2*(eta1^2*gr0^2+eta1^2*gr1^2+eta1^2*gr2^2+gr0^2*gr1^2+gr0^2*gr2^2+gr1^2*gr2^2))'
    derivative_order = 2
  [../]
  [./probability]
    # This is a made up toy nucleation rate it should be replaced by
    # classical nucleation theory in a real simulation.
    type = ParsedMaterial
    property_name = P
    coupled_variables = bnds
    expression = max((0.65-bnds)*1e-4,0)
    outputs = exodus
  [../]
  [./nucleation]
    # The nucleation material is configured to insert nuclei into the free energy
    # tht force the concentration to go to 0.95, and holds this enforcement for 500
    # time units.
    type = DiscreteNucleation
    property_name = Fn
    op_names  = c
    op_values = 0.90
    penalty = 5
    penalty_mode = MIN
    map = map
    outputs = exodus
  [../]
    [./nucleation2]
      # The nucleation material is configured to insert nuclei into the free energy
      # tht force the concentration to go to 0.95, and holds this enforcement for 500
      # time units.
      type = DiscreteNucleation
      property_name = Fn2
      op_names  = eta1
      op_values = 0.95
      penalty = 5
      penalty_mode = MIN
      map = map2
      outputs = exodus
    [../]
  [./free_energy]
    # add the chemical and nucleation free energy contributions together
    type = DerivativeSumMaterial
    derivative_order = 2
    coupled_variables = 'c bnds eta1 gr0 gr1 gr2'
    sum_materials = 'f_chem Fn Fn2'
  [../]
[]

[UserObjects]
  [./inserter]
    # The inserter runs at the end of each time step to add nucleation events
    # that happend during the timestep (if it converged) to the list of nuclei
    type = DiscreteNucleationInserter
    hold_time = 50
    probability = P
    radius = 3
  [../]
  [./map]
    # The map UO runs at the beginning of a timestep and generates a per-element/qp
    # map of nucleus locations. The map is only regenerated if the mesh changed or
    # the list of nuclei was modified.
    # The map converts the nucleation points into finite area objects with a given radius.
    type = DiscreteNucleationMap
    periodic = gr2
    inserter = inserter
  [../]
    [./inserter2]
      # The inserter runs at the end of each time step to add nucleation events
      # that happend during the timestep (if it converged) to the list of nuclei
      type = DiscreteNucleationInserter
      hold_time = 100
      probability = P
      radius = 2
    [../]
    [./map2]
      # The map UO runs at the beginning of a timestep and generates a per-element/qp
      # map of nucleus locations. The map is only regenerated if the mesh changed or
      # the list of nuclei was modified.
      # The map converts the nucleation points into finite area objects with a given radius.
      type = DiscreteNucleationMap
      periodic = eta1
      inserter = inserter2
    [../]
  [voronoi]
    type = PolycrystalVoronoi
    rand_seed = 1
    int_width = 15
  []
  [grain_tracker]
    type = GrainTracker
  []
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[BCs]
  [./Periodic]
    [./all]
      auto_direction = 'x y'
    [../]
  [../]
[]

[Executioner]
  type = Transient
  scheme = bdf2
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -sub_pc_type'
  petsc_options_value = 'asm      lu          '

  nl_max_its = 20
  l_tol = 1.0e-4
  nl_rel_tol = 1.0e-10
  nl_abs_tol = 1.0e-10
  start_time = 0.0
  num_steps = 5000
  dt = 1
[]

[Outputs]
  exodus = true
[]
