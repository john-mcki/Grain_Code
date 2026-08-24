!! ---------------------------------------------------------------------------------------------------------------------------------
!! file:      corr_sim_param_mod.f90
!! authors:   Alexander Chadwick (University of Michigan)
!! date:      October 20, 2017
!! email:     afchadwi@umich.edu
!! std:       Fortran 2008
!!
!! This file contains all of the parameter modules for the 3D corrosion simulation for LIFT. It also contains all of the subroutines
!! needed for reading the input file and setting these parameters. Refer to the README for more information.
!! ---------------------------------------------------------------------------------------------------------------------------------

module variables
  use, intrinsic :: iso_fortran_env, only : real32, real64, int32, int64

  implicit none

  !! ----- Set public/private attributes for the subroutines ----- !!
  public :: read_input_file
  private :: read_vec_arg
  private :: read_vec_arg_i32, read_vec_arg_i64
  private :: read_vec_arg_f32, read_vec_arg_f64
  private :: read_vec_arg_char
  private :: read_vec_arg_logical
  interface read_vec_arg
    module procedure read_vec_arg_i32, read_vec_arg_i64
    module procedure read_vec_arg_f32, read_vec_arg_f64
    module procedure read_vec_arg_char
    module procedure read_vec_arg_logical
  end interface

  !! ----- Basic variables related to the simulation ----- !!
  character(len=:), allocatable :: run_name              !! Root name of the folder where data will be stored
  logical :: restart                                          !! If true, we are restarting a previous simulation!
  logical :: pit_simulation                       !! Is this a pit simulation? Sets some parts of BC
  logical :: well_stirred_ion
  logical :: well_stirred_elect
  logical :: supporting_elect
  logical :: solve_conc
  logical :: solve_pot
  logical :: solve_cur
  logical :: disable_deposition
  logical :: show_progress
  integer(int32) :: master_seed                            !! Master seed for RNG

  !! ----- Parameters for MPI ----- !!
  integer(int32), dimension(1:3) :: ndivs                !! Processors in each dimension
  logical, dimension(1:3) :: periodic   !! Periodic BC toggle

  !! ----- Variables for microstructure generation ----- !!
  character(len=3) :: mxc_mode                              !! Mode for the microstructure generation (vor = voronoi)
  character(len=3) :: int_mode                              !! Mode for the interfacial shape
  logical :: pre_smooth                                     !! Whether or not to use voxel (1/0) or smooth (tanh) microstructures
  integer(int32) :: pit_rad                                 !! Radius of the initial pit, in grid points
  integer(int32) :: n_per                                   !! Number of complete periods for the sinusoidal perturbations
  real(real64) :: s_amp                                     !! Amplitude of the sinusoidal perturbations

  !! ----- Spatial and temporal discretization parameters ----- !!
  !! Size of Numerical Grid / Mesh
  integer(int32) :: nx
  integer(int32) :: ny
  integer(int32) :: nz
  integer(int32) :: lithSize
  integer(int32) :: prisSize
  real(real64) :: dh                            !! Grid spacing (nm)
  real(real64) :: dh2                                 !! Double the grid spacing
  real(real64) :: dhsq                                   !! Grid spacing squared
  real(real64) :: dt                                !! Time step size (s)
  integer(int32) :: tsteps                               !! Total number of time steps
  integer(int32) :: out_freq                          !! How often to write data to disk

  !! ----- Parameters for iterative solvers ----- !!
  logical :: SimultaneousSolverMode                     !! If true, all variables are solved simultaneously
  logical :: ConUseSOR
  integer(int32) :: MaxIters                           !! Maximum number of iterations for any solver
  real(real64) :: MaxError_C                           !! Max residual for concentration solver(s)
  real(real64) :: Omega_C                           !! Over-relaxation paramter ( 0 < omega < 2 )
  real(real64) :: MaxError_V                           !! Max residual for potential solver(s)
  real(real64) :: Omega_V                            !! Over-relaxation parameter for potential (between 0 and 2)
  real(real64) :: MaxError_P                           !! Max difference in iteration for order parameter  solver(s)
  real(real64) :: Omega_P
  real(real64) :: MaxError_D
  real(real64) :: Omega_D
  real(real64) :: SBMCutoff                             !! Order parameter cutoff value for SBM solvers
  real(real64) :: PhiCutoff                             !! Order parameter cutoff for AC/CH solvers
  integer(int32) :: res_skip

  !! ----- General physical constants ----- !!
  real(real64), parameter :: pi = ACOS(-1.0d0)                              !! Pi
  real(real64), parameter :: FrdCst = 96485.33289d0                         !! Faraday's constant (C/mol)
  real(real64), parameter :: GasCst = 8.3144598d0                           !! Ideal gas constatnt (J/mol K)
  real(real64) :: Temp                                 !! Absolute temperature (K)
  real(real64) :: FoRT                   !! F/RT

  !! ----- Phase field parameters ----- !!
  integer(int32) :: np                                       !! Number of solid order parameters
  real(real64) :: Ms                                      !! Cahn-Hilliard mobility (nm^2/s)
  character(len=1) :: mob_mode                            !! Toggle for constant (c) or Peclet-based (p) mobilities
  real(real64) :: eps_sq                                 !! Gradient energy coefficient (nm^2)
  real(real64) :: W                                        !! Double-well potential height
  real(real64) :: gamma                                  !! Interfacial overlap penalty coefficient (nm)

  !! Boundary Conditions: P = Periodic, F = Fixed at Initial Value, N = No-Flux
  character(len=1) :: Sol_W
  character(len=1) :: Sol_E
  character(len=1) :: Sol_N
  character(len=1) :: Sol_S
  character(len=1) :: Sol_U
  character(len=1) :: Sol_D

  !! ----- Concentration solver/transport properties ----- !!
  real(real64), parameter :: rRxn_units = 1.0d24                            !! Conversion from mol/nm^3 to mol/L
  real(real64) :: Vm                                  !! Molar volume of metal (nm^3/mol)
  integer(int32) :: ns                                  !! Number of NON-REFERENCE species in the electrolyte
  real(real64), allocatable :: z_elec(:)              !! Array for the charges of NON-REFERENCE species in the electrolyte
  real(real64), allocatable :: D_elec(:)              !! Array for the diffusivities of NON-REFERENCE species in electrolyte
  real(real64), allocatable :: c_bulk(:)              !! Array for the bulk concentration of NON-REFERENCE species in electrolyte
  real(real64), allocatable :: c_sat(:)               !! Array for the saturation concentrations of NON-REFERENCE species in electrolyte
  real(real64) :: z_ref                               !! Charge of REFERENCE species in electrolyte
  real(real64) :: D_ref                               !! Diffusivity of REFERENCE species in electrolyte
  real(real64) :: tp                                  !! Lithium transference number

  !! Ion Concentration Boundary Conditions: F = Fixed at Initial Value, N = No-Flux
  character(len=1) :: Ion_W
  character(len=1) :: Ion_E
  character(len=1) :: Ion_N
  character(len=1) :: Ion_S
  character(len=1) :: Ion_U
  character(len=1) :: Ion_D
  integer(int32) :: SBM_Robin_BC

  !! Electrostatic Potential Boundary Conditions: F = Fixed at Initial Value, N = No-Flux
  character(len=1) :: Pot_W
  character(len=1) :: Pot_E
  character(len=1) :: Pot_N
  character(len=1) :: Pot_S
  character(len=1) :: Pot_U
  character(len=1) :: Pot_D

  !! ----- Butler-Volmer kinetic properties (this version should work for multi-order parameter kinetics)  ----- !!
  real(real64) :: area_relative                       !! Relative increase in surface area
  character(len=3) :: sim_mode                        !! Toggle for galvanostatic vs potentiostatic
  real(real64) :: i_target                     !! Target total current density for the bisection search (A/???)
  real(real64) :: i_tolerance                   !! Tolerance on the bisection search (A/???)
  real(real64) :: vWE_init                                   !! Working electrode potential (V)
  real(real64) :: ne
  real(real64) :: beta
  real(real64) :: vOC
  real(real64) :: kDep
  real(real64) :: kDis
  real(real64) :: kAmp

  !! ----- Initial Allen-Cahn smoothing parameters ----- !!
  integer(int32) :: ssteps                                  !! Number of smoothing steps
  real(real64) :: sLdt                              !! Allen-Cahn mobility for smoothing (1/s)
  logical :: smooth_Li                                    !! Logical for whether or not we smooth OP 0

  !! ----- Traction-free boundary configuration ----- !!
  logical :: BC_TF_top, BC_TF_bot                         !! Toggle for whether we are using traction free boundaries at the top and bottom
  integer(int32) :: TF_size                               !! Size of the traction free layer(s)
  integer(int32) :: PadTop, PadBot                        !! Padding indices for the microstructure generation
  real(real64) :: zl, zh, zeta

  !! ----- Mechanical properties ----- !!
  real(real64) :: Emod_sub, nu_sub, Az_sub                !! Young's modulus, Poisson's ratio, and Zener anisotropy for substrate
  real(real64) :: Emod_lyr, nu_lyr, Az_lyr                !! ... for the protection layer
  real(real64), dimension(1:21) :: Cel_sub, Cel_lyr       !! Elastic constant vectors
  real(real64), dimension(1:6) :: Bel_sub, Bel_lyr        !! Volumetric expansion vectors
  logical :: Calc_Cel_vec_sub, Calc_Cel_vec_lyr           !! Logicals that toggle depending on whether we've passed the elastic constants or the vector
  real(real64) :: Euler_sub 
  real(real64), allocatable :: Euler_lyr(:)                    !! Manual override for the Euler angles
  logical :: Euler_sub_override, Euler_lyr_override       !! Logicals that automatically toggle the override

  !! Mechanical Boundary Conditions: F = Fixed at Initial Value, N = No-Flux
  character(len=1) :: Mec_W
  character(len=1) :: Mec_E
  character(len=1) :: Mec_N
  character(len=1) :: Mec_S
  character(len=1) :: Mec_U
  character(len=1) :: Mec_D

!! ----- Subroutines for reading input file and setting parameters -------------------------------------------------------------- !!

contains

  subroutine read_input_file( input_fname )
    !! ~~~~~~~~~~~~~~~~~~~~
    !! This subroutine is for reading in the input file
    !! ~~~~~~~~~~~~~~~~~~~~
    implicit none

    character(len=200), intent(in) :: input_fname
    character(len=200), allocatable :: input_file_args(:)
    character(len=200) :: rname_temp

    integer :: i , j , l , nlines , rr, run_name_len

    !! Set some default behavior
    int_mode = 'flt'    !! Default to flat interface
    Calc_Cel_vec_sub = .true.
    Calc_Cel_vec_lyr = .true.
    Euler_sub_override = .false.
    Euler_lyr_override = .false.
    smooth_Li = .true.

    !! Open the input file
    open( unit=10 , file=input_fname , form='formatted' , action='read' , status='old' )

    !! Count the number of lines in the input file
    nlines = 0
    do
      read(10,*,iostat=rr)
      nlines = nlines + 1
      if ( rr .lt. 0 ) exit
    enddo
    nlines = nlines-1
    !! Rewind the input file to actually read it
    rewind(10)
    !! Read the input file to an array of characters
    allocate( input_file_args(1:nlines) )
    do i = 1, nlines
      read(10,'(A200)') input_file_args(i)
    enddo
    close(10)

    !! Parse the input file array
    do j = 1, nlines
      l = scan( input_file_args(j) , '=' )
      !! ----- Read the general run variables ----- !!
      !! Read the run name
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'run_name' ) then
        read(input_file_args(j)(l+1:200),'(A200)') rname_temp
        run_name_len = len( trim(adjustl(rname_temp)) )
        allocate( character(len=run_name_len) :: run_name )
        run_name = trim(adjustl(rname_temp))
      endif
      !! Read the logical control variables
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'restart' ) read(input_file_args(j)(l+1:200),'(L1)') restart
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'pit_simulation' ) read(input_file_args(j)(l+1:200),'(L1)') pit_simulation
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'well_stirred_ion' ) read(input_file_args(j)(l+1:200),'(L1)') well_stirred_ion
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'well_stirred_elect' ) read(input_file_args(j)(l+1:200),'(L1)') well_stirred_elect
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'supporting_elect' ) read(input_file_args(j)(l+1:200),'(L1)') supporting_elect
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'solve_conc' ) read(input_file_args(j)(l+1:200),'(L1)') solve_conc
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'solve_pot' ) read(input_file_args(j)(l+1:200),'(L1)') solve_pot
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'solve_cur' ) read(input_file_args(j)(l+1:200),'(L1)') solve_cur
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'disable_deposition' ) read(input_file_args(j)(l+1:200),'(L1)') disable_deposition
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'show_progress' ) read(input_file_args(j)(l+1:200),'(L1)') show_progress
      !! Read the master seed for the RNG
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'master_seed' ) read(input_file_args(j)(l+1:200),'(I25)') master_seed

      !! ----- Read the variables related to MPI ----- !!
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'ndivs' ) call read_vec_arg( ndivs, input_file_args(j)(l+1:200) )
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'periodic' ) call read_vec_arg( periodic, input_file_args(j)(l+1:200) )

      !! ----- Read the variables for the microstructure generation ----- !!
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'mxc_mode' ) then
        read(input_file_args(j)(l+1:200),'(A200)') rname_temp
        mxc_mode = trim(adjustl(rname_temp))
      endif
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'int_mode' ) then
        read(input_file_args(j)(l+1:200),'(A200)') rname_temp
        int_mode = trim(adjustl(rname_temp))
      endif
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'pre_smooth' ) read(input_file_args(j)(l+1:200), *) pre_smooth
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'pit_rad' ) read(input_file_args(j)(l+1:200),'(I25)') pit_rad
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'n_per' ) read(input_file_args(j)(l+1:200), * ) n_per
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 's_amp' ) read(input_file_args(j)(l+1:200), * ) s_amp

      !! ----- Read the variables related to the domain size/discretization ----- !!
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'dh' ) read(input_file_args(j)(l+1:200),'(E25.12)') dh
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'nx' ) read(input_file_args(j)(l+1:200),'(I25)') nx
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'ny' ) read(input_file_args(j)(l+1:200),'(I25)') ny
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'nz' ) read(input_file_args(j)(l+1:200),'(I25)') nz
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'lithSize' ) read(input_file_args(j)(l+1:200),'(I25)') lithSize
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'prisSize' ) read(input_file_args(j)(l+1:200),'(I25)') prisSize

      !! ----- Read the temporal discretization parameters ----- !!
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'dt' ) read(input_file_args(j)(l+1:200),'(E25.12)') dt
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'tsteps' ) read(input_file_args(j)(l+1:200),'(I25)') tsteps
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'out_freq' ) read(input_file_args(j)(l+1:200),'(I25)') out_freq

      !! ----- Read the iterative solver parameters ----- !!
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'SimultaneousSolverMode' ) &
        read(input_file_args(j)(l+1:200),'(L1)') SimultaneousSolverMode
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'ConUseSOR' ) &
        read(input_file_args(j)(l+1:200),'(L1)') ConUseSOR
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'MaxIters' ) read(input_file_args(j)(l+1:200),'(I25)') MaxIters
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'MaxError_C' ) read(input_file_args(j)(l+1:200),'(E25.12)') MaxError_C
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'Omega_C' ) read(input_file_args(j)(l+1:200),'(E25.12)') Omega_C
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'MaxError_V' ) read(input_file_args(j)(l+1:200),'(E25.12)') MaxError_V
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'Omega_V' ) read(input_file_args(j)(l+1:200),'(E25.12)') Omega_V
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'MaxError_P' ) read(input_file_args(j)(l+1:200),'(E25.12)') MaxError_P
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'Omega_P' ) read(input_file_args(j)(l+1:200),'(E25.12)') Omega_P
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'SBMCutoff' ) read(input_file_args(j)(l+1:200),'(E25.12)') SBMCutoff
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'PhiCutoff' ) read(input_file_args(j)(l+1:200),'(E25.12)') PhiCutoff
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'res_skip' ) read(input_file_args(j)(l+1:200),'(I25)') res_skip

      !! ----- Read the physical constants ----- !!
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'Temp' ) read(input_file_args(j)(l+1:200),'(E25.12)') Temp

      !! ----- Read the phase-field parameters ----- !!
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'np' ) then
        read(input_file_args(j)(l+1:200),'(I25)') np
        allocate( Euler_lyr(1:np) )
      endif
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'Ms' ) read(input_file_args(j)(l+1:200),'(E25.12)') Ms
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'mob_mode' ) read(input_file_args(j)(l+1:200),'(A1)') mob_mode
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'eps_sq' ) read(input_file_args(j)(l+1:200),'(E25.12)') eps_sq
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'W' ) read(input_file_args(j)(l+1:200),'(E25.12)') W
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'gamma' ) read(input_file_args(j)(l+1:200),'(E25.12)') gamma

      !! ----- Read the phase-field boundary conditions ----- !!
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'Sol_W' ) read(input_file_args(j)(l+1:200),'(A1)') Sol_W
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'Sol_E' ) read(input_file_args(j)(l+1:200),'(A1)') Sol_E
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'Sol_N' ) read(input_file_args(j)(l+1:200),'(A1)') Sol_N
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'Sol_S' ) read(input_file_args(j)(l+1:200),'(A1)') Sol_S
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'Sol_U' ) read(input_file_args(j)(l+1:200),'(A1)') Sol_U
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'Sol_D' ) read(input_file_args(j)(l+1:200),'(A1)') Sol_D

      !! ----- Read parameters for the transport solver ----- !!
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'Vm' ) read(input_file_args(j)(l+1:200),*) Vm
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'ns' ) then
        read(input_file_args(j)(l+1:200),*) ns
        allocate( z_elec(1:ns), D_elec(1:ns), c_bulk(1:ns), c_sat(1:ns) )
      endif
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'z_elec' ) call read_vec_arg( z_elec, input_file_args(j)(l+1:200) )
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'D_elec' ) call read_vec_arg( D_elec, input_file_args(j)(l+1:200) )
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'c_bulk' ) call read_vec_arg( c_bulk, input_file_args(j)(l+1:200) )
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'c_sat' ) call read_vec_arg( c_sat, input_file_args(j)(l+1:200) )
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'z_ref' ) read(input_file_args(j)(l+1:200), * ) z_ref
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'D_ref' ) read(input_file_args(j)(l+1:200), * ) D_ref
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'tp' ) read(input_file_args(j)(l+1:200), * ) tp

      !! ----- Read the phase-field boundary conditions ----- !!
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'Ion_W' ) read(input_file_args(j)(l+1:200),'(A1)') Ion_W
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'Ion_E' ) read(input_file_args(j)(l+1:200),'(A1)') Ion_E
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'Ion_N' ) read(input_file_args(j)(l+1:200),'(A1)') Ion_N
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'Ion_S' ) read(input_file_args(j)(l+1:200),'(A1)') Ion_S
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'Ion_U' ) read(input_file_args(j)(l+1:200),'(A1)') Ion_U
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'Ion_D' ) read(input_file_args(j)(l+1:200),'(A1)') Ion_D
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'SBM_Robin_BC' ) read(input_file_args(j)(l+1:200),'(I25)') SBM_Robin_BC

      !! ----- Read the phase-field boundary conditions ----- !!
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'Pot_W' ) read(input_file_args(j)(l+1:200),'(A1)') Pot_W
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'Pot_E' ) read(input_file_args(j)(l+1:200),'(A1)') Pot_E
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'Pot_N' ) read(input_file_args(j)(l+1:200),'(A1)') Pot_N
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'Pot_S' ) read(input_file_args(j)(l+1:200),'(A1)') Pot_S
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'Pot_U' ) read(input_file_args(j)(l+1:200),'(A1)') Pot_U
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'Pot_D' ) read(input_file_args(j)(l+1:200),'(A1)') Pot_D

      !! ----- Read the Butler-Volmer kinetic parameters ----- !!
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'sim_mode' ) then
        read(input_file_args(j)(l+1:200),'(A200)') rname_temp
        sim_mode = trim(adjustl(rname_temp))
      endif
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'i_target' ) read(input_file_args(j)(l+1:200),'(E25.12)') i_target
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'i_tolerance' ) read(input_file_args(j)(l+1:200),'(E25.12)') i_tolerance
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'vWE_init' ) read(input_file_args(j)(l+1:200),'(E25.12)') vWE_init
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'ne' ) read(input_file_args(j)(l+1:200),'(E25.12)') ne
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'beta' ) read(input_file_args(j)(l+1:200),'(E25.12)') beta
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'vOC' ) read(input_file_args(j)(l+1:200),'(E25.12)') vOC
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'kDep' ) read(input_file_args(j)(l+1:200),'(E25.12)') kDep
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'kDis' ) read(input_file_args(j)(l+1:200),'(E25.12)') kDis
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'kAmp' ) read(input_file_args(j)(l+1:200),'(E25.12)') kAmp

      !! ----- Read Allen-Cahn smoothing parameters ----- !!
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'ssteps' ) read(input_file_args(j)(l+1:200),'(I25)') ssteps
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'sLdt' ) read(input_file_args(j)(l+1:200),'(E25.12)') sLdt
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'smooth_Li' ) read(input_file_args(j)(l+1:200), *) smooth_Li

      !! ----- Read the traction-free configuration ----- !!
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'BC_TF_top' ) read(input_file_args(j)(l+1:200),'(L1)') BC_TF_top
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'BC_TF_bot' ) read(input_file_args(j)(l+1:200),'(L1)') BC_TF_bot
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'TF_size' ) read(input_file_args(j)(l+1:200),'(I25)') TF_size

      !! ----- Read the mechanical properties ----- !!
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'Emod_sub' ) read(input_file_args(j)(l+1:200), * ) Emod_sub
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'nu_sub' ) read(input_file_args(j)(l+1:200), * ) nu_sub
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'Az_sub' ) read(input_file_args(j)(l+1:200), * ) Az_sub
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'Emod_lyr' ) read(input_file_args(j)(l+1:200), * ) Emod_lyr
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'nu_lyr' ) read(input_file_args(j)(l+1:200), * ) nu_lyr
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'Az_lyr' ) read(input_file_args(j)(l+1:200), * ) Az_lyr
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'Cel_sub' ) then
        call read_vec_arg( Cel_sub, input_file_args(j)(l+1:200) )
        Calc_Cel_vec_sub = .false.
      endif
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'Cel_lyr' ) then
        call read_vec_arg( Cel_lyr, input_file_args(j)(l+1:200) )
        Calc_Cel_vec_lyr = .false.
      endif
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'Bel_sub' ) call read_vec_arg( Bel_sub, input_file_args(j)(l+1:200) )
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'Bel_lyr' ) call read_vec_arg( Bel_lyr, input_file_args(j)(l+1:200) )
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'Euler_sub' ) then
        read(input_file_args(j)(l+1:200), * ) Euler_sub
        Euler_sub_override = .true.
      endif
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'Euler_lyr' ) then
        call read_vec_arg( Euler_lyr, input_file_args(j)(l+1:200))
        Euler_lyr_override = .true.
      endif


      !! ----- Read the mechanical boundary conditions ----- !!
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'Mec_W' ) read(input_file_args(j)(l+1:200),'(A1)') Mec_W
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'Mec_E' ) read(input_file_args(j)(l+1:200),'(A1)') Mec_E
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'Mec_N' ) read(input_file_args(j)(l+1:200),'(A1)') Mec_N
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'Mec_S' ) read(input_file_args(j)(l+1:200),'(A1)') Mec_S
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'Mec_U' ) read(input_file_args(j)(l+1:200),'(A1)') Mec_U
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'Mec_D' ) read(input_file_args(j)(l+1:200),'(A1)') Mec_D

      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'MaxError_D' ) read(input_file_args(j)(l+1:200),'(E25.18)') MaxError_D
      if ( adjustl(input_file_args(j)(1:l-1)) .eq. 'Omega_D' ) read(input_file_args(j)(l+1:200),'(E25.18)') Omega_D


    enddo

    deallocate( input_file_args )

    !! Set the remaining simulation parameters that aren't controlled by input file
    call set_remaining_parameters()

  end subroutine read_input_file

  subroutine set_remaining_parameters()
    !! This subroutine sets variables that need other parameters but aren't controlled in the input file
    implicit none

    integer(int32) :: i
    real(real64) :: C11_t, C12_t, C44_t
    real(real64) :: dc

    dh2 = 2.0d0*dh                                 !! Double the grid spacing
    dhsq = dh*dh                                   !! Grid spacing squared
    FoRT = FrdCst/(GasCst * Temp)                  !! F/RT

    if ( sim_mode .eq. 'cur' ) i_tolerance = i_tolerance*i_target

    !! Set some transport properties based on the transference number
    dc = c_bulk(1)*(1.d0-tp)*D_elec(1) / ( -(1.d0-tp)*(D_elec(2)-D_elec(1)) + tp*D_ref )
    c_bulk(1) = c_bulk(1) - dc
    c_bulk(2) = dc

    !! Determine if we are using any traction free boundary and configure variables accordingly
    PadTop = 0
    PadBot = 0
    if ( BC_TF_top .eqv. .true. ) then
      PadTop = TF_size
      nz = nz + TF_size
    endif

    if ( BC_TF_bot .eqv. .true. ) then
      PadBot = TF_size
      nz = nz + TF_size
    endif

    zeta = sqrt(2.d0*eps_sq/W)
    zl = -100.d0*dh
    zh = 100.d0*dh + real(nz,real64)*dh
    if ( BC_TF_bot .eqv. .true. ) zl = real(PadBot,real64)*dh
    if ( BC_TF_top .eqv. .true. ) zh = real(nz-PadTop,real64)*dh

    if ( Calc_Cel_vec_sub ) then
      Cel_sub = 0.d0
      C11_t = - Emod_sub*(11.d0 + 13.d0*Az_sub + Az_sub*Az_sub - 19.d0*nu_sub - 7.d0*Az_sub*nu_sub + Az_sub*Az_sub*nu_sub) &
        / ((3.d0 + 19.d0*Az_sub + 3.d0*Az_sub*Az_sub)*(2.d0*nu_sub*nu_sub + nu_sub - 1.d0))
      C12_t = - Emod_sub*(-4.d0 + 3.d0*Az_sub + Az_sub*Az_sub + 11.d0*nu_sub + 13.d0*Az_sub*nu_sub + Az_sub*Az_sub*nu_sub) &
        / ((3.d0 + 19.d0*Az_sub + 3.d0*Az_sub*Az_sub)*(2.d0*nu_sub*nu_sub + nu_sub - 1.d0))
      C44_t = Az_sub*(C11_t - C12_t)/2.d0
      Cel_sub(1:3) = C11_t
      Cel_sub(4:6) = C44_t
      Cel_sub(7:8) = C12_t
      Cel_sub(12) = C12_t
    endif
    if ( Calc_Cel_vec_lyr ) then
      Cel_lyr = 0.d0
      C11_t = - Emod_lyr*(11.d0 + 13.d0*Az_lyr + Az_lyr*Az_lyr - 19.d0*nu_lyr - 7.d0*Az_lyr*nu_lyr + Az_lyr*Az_lyr*nu_lyr) &
        / ((3.d0 + 19.d0*Az_lyr + 3.d0*Az_lyr*Az_lyr)*(2.d0*nu_lyr*nu_lyr + nu_lyr - 1.d0))
      C12_t = - Emod_lyr*(-4.d0 + 3.d0*Az_lyr + Az_lyr*Az_lyr + 11.d0*nu_lyr + 13.d0*Az_lyr*nu_lyr + Az_lyr*Az_lyr*nu_lyr) &
        / ((3.d0 + 19.d0*Az_lyr + 3.d0*Az_lyr*Az_lyr)*(2.d0*nu_lyr*nu_lyr + nu_lyr - 1.d0))
      C44_t = Az_lyr*(C11_t - C12_t)/2.d0
      Cel_lyr(1:3) = C11_t
      Cel_lyr(4:6) = C44_t
      Cel_lyr(7:8) = C12_t
      Cel_lyr(12) = C12_t
    endif

    print *, sLdt, W, eps_sq, np

  end subroutine set_remaining_parameters

  subroutine read_vec_arg_i32 ( var , arg )
    !! This subroutine reads in a vector-type parameter (i.e. the diffusivities of multiple species).
    !! If more arguments are given than there are needed, the last arguments are ignored.
    !! If fewer arguments are given than there are needed, the last argument that is given is copied to all of the following
    !! parameters. You have been warned...
    !! This is the INT32 version of the subroutine that will be in an overloaded interface.

    implicit none

    integer(INT32), intent(inout), dimension(1:) :: var
    character(len=*), intent(in) :: arg

    integer :: i, l, m

    i = 0
    l = 1
    m = 1

    do
      i = i + 1

      if ( i .gt. size(var) ) exit

      m = scan( arg(l:), ',' )

      if ( m .ne. 0 ) then
        read ( arg(l:l+m-2), * ) var(i)
      else
        read ( arg(l:), * ) var(i)
        exit
      endif

      l = l + m
    enddo

    if ( i .lt. size(var) ) var(i+1:) = var(i)

  end subroutine read_vec_arg_i32

  subroutine read_vec_arg_i64 ( var , arg )
    !! This subroutine reads in a vector-type parameter (i.e. the diffusivities of multiple species).
    !! If more arguments are given than there are needed, the last arguments are ignored.
    !! If fewer arguments are given than there are needed, the last argument that is given is copied to all of the following
    !! parameters. You have been warned...
    !! This is the INT64 version of the subroutine that will be in an overloaded interface.

    implicit none

    integer(INT64), intent(inout), dimension(1:) :: var
    character(len=*), intent(in) :: arg

    integer :: i, l, m

    i = 0
    l = 1
    m = 1

    do
      i = i + 1

      if ( i .gt. size(var) ) exit

      m = scan( arg(l:), ',' )

      if ( m .ne. 0 ) then
        read ( arg(l:l+m-2), * ) var(i)
      else
        read ( arg(l:), * ) var(i)
        exit
      endif

      l = l + m
    enddo

    if ( i .lt. size(var) ) var(i+1:) = var(i)

  end subroutine read_vec_arg_i64

  subroutine read_vec_arg_f32 ( var , arg )
    !! This subroutine reads in a vector-type parameter (i.e. the diffusivities of multiple species).
    !! If more arguments are given than there are needed, the last arguments are ignored.
    !! If fewer arguments are given than there are needed, the last argument that is given is copied to all of the following
    !! parameters. You have been warned...
    !! This is the REAL32 version of the subroutine that will be in an overloaded interface.

    implicit none

    real(REAL32), intent(inout), dimension(1:) :: var
    character(len=*), intent(in) :: arg

    integer :: i, l, m

    i = 0
    l = 1
    m = 1

    do
      i = i + 1

      if ( i .gt. size(var) ) exit

      m = scan( arg(l:), ',' )

      if ( m .ne. 0 ) then
        read ( arg(l:l+m-2), * ) var(i)
      else
        read ( arg(l:), * ) var(i)
        exit
      endif

      l = l + m
    enddo

    if ( i .lt. size(var) ) var(i+1:) = var(i)

  end subroutine read_vec_arg_f32

  subroutine read_vec_arg_f64 ( var , arg )
    !! This subroutine reads in a vector-type parameter (i.e. the diffusivities of multiple species).
    !! If more arguments are given than there are needed, the last arguments are ignored.
    !! If fewer arguments are given than there are needed, the last argument that is given is copied to all of the following
    !! parameters. You have been warned...
    !! This is the REAL64 version of the subroutine that will be in an overloaded interface.

    implicit none

    real(REAL64), intent(inout), dimension(1:) :: var
    character(len=*), intent(in) :: arg

    integer :: i, l, m

    i = 0
    l = 1
    m = 1

    do
      i = i + 1

      if ( i .gt. size(var) ) exit

      m = scan( arg(l:), ',' )

      if ( m .ne. 0 ) then
        read ( arg(l:l+m-2), * ) var(i)
      else
        read ( arg(l:), * ) var(i)
        exit
      endif

      l = l + m
    enddo

    if ( i .lt. size(var) ) var(i+1:) = var(i)

  end subroutine read_vec_arg_f64

  subroutine read_vec_arg_char ( var , arg )
    !! This subroutine reads in a vector-type parameter (i.e. the diffusivities of multiple species).
    !! If more arguments are given than there are needed, the last arguments are ignored.
    !! If fewer arguments are given than there are needed, the last argument that is given is copied to all of the following
    !! parameters. You have been warned...
    !! This is the REAL128 version of the subroutine that will be in an overloaded interface.

    implicit none

    character(len=*), intent(inout), dimension(1:) :: var
    character(len=*), intent(in) :: arg

    integer :: i, l, m

    i = 0
    l = 1
    m = 1

    do
      i = i + 1

      if ( i .gt. size(var) ) exit

      m = scan( arg(l:), ',' )

      if ( m .ne. 0 ) then
        read ( arg(l:l+m-2), * ) var(i)
      else
        read ( arg(l:), * ) var(i)
        exit
      endif

      l = l + m
    enddo

    if ( i .lt. size(var) ) var(i+1:) = var(i)

  end subroutine read_vec_arg_char

  subroutine read_vec_arg_logical ( var , arg )
    !! This subroutine reads in a vector-type parameter (i.e. the diffusivities of multiple species).
    !! If more arguments are given than there are needed, the last arguments are ignored.
    !! If fewer arguments are given than there are needed, the last argument that is given is copied to all of the following
    !! parameters. You have been warned...
    !! This is the REAL128 version of the subroutine that will be in an overloaded interface.

    implicit none

    logical, intent(inout), dimension(1:) :: var
    character(len=*), intent(in) :: arg

    integer :: i, l, m

    i = 0
    l = 1
    m = 1

    do
      i = i + 1

      if ( i .gt. size(var) ) exit

      m = scan( arg(l:), ',' )

      if ( m .ne. 0 ) then
        read ( arg(l:l+m-2), * ) var(i)
      else
        read ( arg(l:), * ) var(i)
        exit
      endif

      l = l + m
    enddo

    if ( i .lt. size(var) ) var(i+1:) = var(i)

  end subroutine read_vec_arg_logical

end module variables
