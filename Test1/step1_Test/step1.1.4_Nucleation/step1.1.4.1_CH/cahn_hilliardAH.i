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
    # v = ''
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
        free_energy = F
        mobility = M
        kappa = kappa_c
        solve_type = REVERSE_SPLIT
      [../]
    [../]
  [../]
[]

[ICs]
  [./c_IC]
    type = RandomIC
    variable = c
    min = 0.2
    max = 0.21
  [../]
  [PolycrystalICs]
    [PolycrystalColoringIC]
      polycrystal_ic_uo = voronoi
    []
  []
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
    prop_names  = 'M kappa_c'
    prop_values = '1 0.1'
  [../]
  [./chemical_free_energy]
    # simple double well free energy
    type = DerivativeParsedMaterial
    property_name = Fc
    coupled_variables = 'c'
    constant_names       = 'barr_height  cv_eq'
    constant_expressions = '0.1          0'
    expression = 16*barr_height*c^2*(1-c)^2 # +0.01*(c*plog(c,0.005)+(1-c)*plog(1-c,0.005))
    derivative_order = 2
    outputs = exodus
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
  [./free_energy]
    # add the chemical and nucleation free energy contributions together
    type = DerivativeSumMaterial
    derivative_order = 2
    coupled_variables = 'c bnds'
    sum_materials = 'Fc Fn'
  [../]
[]

[UserObjects]
  [./inserter]
    # The inserter runs at the end of each time step to add nucleation events
    # that happend during the timestep (if it converged) to the list of nuclei
    type = DiscreteNucleationInserter
    hold_time = 50
    probability = P
    radius = 1
  [../]
  [./map]
    # The map UO runs at the beginning of a timestep and generates a per-element/qp
    # map of nucleus locations. The map is only regenerated if the mesh changed or
    # the list of nuclei was modified.
    # The map converts the nucleation points into finite area objects with a given radius.
    type = DiscreteNucleationMap
    periodic = c
    inserter = inserter
  [../]
  [voronoi]
    type = PolycrystalVoronoi
    grain_num = 3 # Number of grains
    rand_seed = 2
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
