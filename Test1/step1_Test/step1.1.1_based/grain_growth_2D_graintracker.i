# This simulation predicts GB migration of a 2D copper polycrystal with 100 grains represented with 8 order parameters
# Mesh adaptivity and time step adaptivity are used
# An AuxVariable is used to calculate the grain boundary locations
# Postprocessors are used to record time step and the number of grains

T = 450 # Constant temperature of the simulation (for mobility calculation)
wGB = 14 # Width of the diffuse GB
GBmob0 = 2.5e-6 #m^4(Js) for copper from schonfelder1997molecular bibtex entry
Q = 0.23 #eV for copper from schonfelder1997molecular bibtex entry
GBenergy = 0.708 #J/m^2 from schonfelder1997molecular bibtex entry
length_scale4 = 1e-36  #模拟的尺度为纳米
length_scale2 = 1e-18  #模拟的尺度为纳米
time_scale = 1e-9  #模拟的尺度为纳秒
kb = 8.617343e-5 # Boltzmann constant in eV/K
JtoeV = 6.24150974e18 # Joule to eV conversion


M0 = ${fparse GBmob0*time_scale/JtoeV/length_scale4}
GBmob = ${fparse M0*exp(-Q/T/kb)}
L = ${fparse 4/3*GBmob/wGB}
gamma = 1.5
kappa = ${fparse 3/4*GBenergy*wGB*JtoeV*length_scale2}  #为什么出现了*6.24150974
mu = ${fparse 6*GBenergy/wGB*JtoeV*length_scale2}  #为什么出现了*6.24150974


[Mesh]
  # Mesh block.  Meshes can be read in or automatically generated
  type = GeneratedMesh
  dim = 2 # Problem dimension
  nx = 40 # Number of elements in the x-direction
  ny = 40 # Number of elements in the y-direction
  xmax = 1000 # maximum x-coordinate of the mesh
  ymax = 1000 # maximum y-coordinate of the mesh
  elem_type = QUAD4 # Type of elements used in the mesh
  uniform_refine = 2 # Initial uniform refinement of the mesh
[]

[GlobalParams]
  # Parameters used by several kernels that are defined globally to simplify input file
  op_num = 3 # Number of order parameters used
  var_name_base = gr # Base name of grains
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
  []
[]

[Variables]
  [./PolycrystalVariables]
  [../]
[]
[ICs]
  [PolycrystalICs]
    [PolycrystalColoringIC]
      polycrystal_ic_uo = voronoi
    []
  []
[]
[BCs]
  # Boundary Condition block
  [Periodic]
    [All]
      auto_direction = 'x y' # Makes problem periodic in the x and y directions
    []
  []
[]


[Kernels]
  # Kernel block, where the kernels defining the residual equations are set up.
  [TimeDerivative]
    type = ADTimeDerivative
    variable = 'gr0'
  []
  [ACGrGrPoly]
    type = ADACGrGrMulti
    variable = 'gr0'
    v = 'gr1 gr2'
    gamma_names = 'gamma_asymm gamma_asymm'
  []
  [ACInterface]
    type = ADACInterface
    variable = 'gr0'
    variable_L = false
  []

  [TimeDerivative2]
    type = ADTimeDerivative
    variable = 'gr1'
  []
  [ACGrGrPoly2]
    type = ADACGrGrMulti
    variable = 'gr1'
    v = 'gr0 gr2'
    gamma_names = 'gamma_asymm gamma_asymm'
  []
  [ACInterface2]
    type = ADACInterface
    variable = 'gr1'
    variable_L = false
  []

  [TimeDerivative3]
    type = ADTimeDerivative
    variable = 'gr2'
  []
  [ACGrGrPoly3]
    type = ADACGrGrMulti
    variable = 'gr2'
    v = 'gr0 gr1'
    gamma_names = 'gamma_asymm gamma_asymm'
  []
  [ACInterface3]
    type = ADACInterface
    variable = 'gr2'
    variable_L = false
  []
[]

[UserObjects]
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

[Materials]
[Constants]
  type = ADGenericConstantMaterial
  prop_names = 'L gamma_asymm kappa_op mu'
  prop_values = '${L} ${gamma} ${kappa} ${mu}'
  outputs     = exodus
[]
[]

[Postprocessors]
  # Scalar postprocessors
  [dt]
    # Outputs the current time step
    type = TimestepSize
  []
[]

[Executioner]
  type = Transient # Type of executioner, here it is transient with an adaptive time step
  scheme = bdf2 # 指定时间积分格式为二阶后向差分格式(backward differentiation formula)。

  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'

  # Uses newton iteration to solve the problem.
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'

  l_max_its = 50 # Max number of linear iterations
  l_tol = 1e-4 # Relative tolerance for linear solves
  nl_max_its = 10 # Max number of nonlinear iterations

  end_time = 2000

  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 20 # Initial time step.  In this simulation it changes.
    optimal_iterations = 6 # 使每个时间步的非线性求解需要大约6次迭代才能收敛。
  []

  [Adaptivity]
    # Block that turns on mesh adaptivity. Note that mesh will never coarsen beyond initial mesh (before uniform refinement)
    initial_adaptivity = 2 # Number of times mesh is adapted to initial condition
    refine_fraction = 0.8 # Fraction of high error that will be refined
    coarsen_fraction = 0.05 # Fraction of low error that will coarsened
    max_h_level = 2 # Max number of refinements used, starting from initial mesh (before uniform refinement)
  []
[]

[Outputs]
  exodus = true
  file_base = 'results/grain_growth'  # 自定义文件名基础
[]
