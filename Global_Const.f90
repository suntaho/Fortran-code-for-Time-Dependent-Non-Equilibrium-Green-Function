!======================================================================
!    Global parameters and variables
!======================================================================
    
module globaldef
  implicit none
  integer, parameter :: R_KIND = selected_real_kind(15)
  
  ! constants
  integer :: ncpus = 2                                                ! number of threads
  real(R_KIND), parameter :: pi = 3.1415926535898_R_KIND
  real(R_KIND), parameter :: delta = 0.00001_R_KIND
  real(R_KIND), parameter :: Hat2eV = 27.211396132_R_KIND             ! transform energy unit from Hatree to eV, default input H in eV for present version
  real(R_KIND), parameter :: hpa = 6.5821188926_R_KIND*10.0_R_KIND**(-16)
  real(R_KIND), parameter :: e_over_hpa_pi = 0.00007748091697940347_R_KIND
  real(R_KIND):: omega = 0.0_R_KIND*10.0_R_KIND**(15)                 ! angular frequency of AC voltage
  
  ! user-defined variables
  integer, parameter :: nitgl = 1999                                  ! maximum number of pints for integral function
  integer :: ntstep = 9999
  integer, parameter :: mdtstep = 99                                  ! time step in MD(Gromacs) setup, for inputHS
  integer :: ntstep_pre = 10
  real(R_KIND) :: dt_std = 0.000005_R_KIND                            ! time interval (ps) for standard time-dependent solution
  real(R_KIND) :: dt_pre = 0.000005_R_KIND                            ! time interval (ps) for equilibrium of initial charge density 
  real(R_KIND) :: Vb = 0.05_R_KIND                                     ! voltage for symmetric bias (V) 
  real(R_KIND) :: VLac = 0.0_R_KIND                                   ! amplitude of AC voltage for L lead
  real(R_KIND) :: VRac = 0.0_R_KIND                                   ! amplitude of AC voltage for R lead
  real(R_KIND) :: tau = 0.0002_R_KIND                                 ! rise time for lead bias (ps)  
  real(R_KIND) :: erng = 1.2_R_KIND, deng = 0.01_R_KIND               ! for spectrum of transmission in steady state analysis (eV)
  logical :: genTf = .true.                                          ! calculate transmission function for steady state or not
  logical :: TdHS = .false.                                           ! use time-dependent hamiltonian/overlap matrix or time-independent (initial) one, not work in windows version
  logical :: NADi = .false.                                           ! consider electron-nucleus (non-adiabatic) coupling or not, not work in present version
  logical :: readrho = .false.                                        ! read inital rho from existed file
  real(R_KIND) :: fEN = 1.0_R_KIND                                    ! scaling the strength of electron-nucleus (non-adiabatic) coupling
  real(R_KIND) :: elow = -10.0_R_KIND, ehigh, dee = 0.01_R_KIND       ! grid for energy integration 
  real(R_KIND) :: kbT = (297.0_R_KIND)*8.61734215_R_KIND*10.0_R_KIND**(-5)     ! thermal energy (eV)
  integer :: DForder = 1                                              ! the order for solving differential equation, 1-5 (1-4 for BDF method)
  logical :: BDF = .false.                                            ! Backward differentiation method(BDF), if .false., use Adams–Bashforth method   
  logical :: orthbol = .true.                                         ! generator hamiltoain having orthogonal basis by given h/s matrix, set true in this version
  logical :: t0_intgl_search = .true.                                 ! find integral subdivision for t0 functions by math. library
  logical :: t_intgl_search = .false.                                 ! find integral subdivision for t functions by math. library
  logical :: cv_model = .true.                                        ! switch to use capacitance network model or not, note that H/S will use time-independent (t=0) values in this version   
  ! switch integral approach
  integer :: diag_int = 1                                             ! (0): slow adaptive Gauss–Legendre quadrature for integration; SLOW speed (only consider states near fermi level)
                                                                      ! (1): diagonal integral with fermi-dirac fun; MEDIUM speed and stable
																	  ! (2): diagonal integral with zero-temperature fermi-dirac fun; FAST speed
  ! Gauss–Legendre quadrature 
  integer :: GLpnt = 3                                                ! Gauss rule for numerical integral of functions: 1-5
  
  ! system variable
  integer :: time(0:1), rate                                          ! timing
  character(1), parameter :: creturn=achar(13), tab=achar(9)          ! ascii character  
  real(R_KIND), parameter :: nil = 0.0_R_KIND, one = 1.0_R_KIND       ! frequently used numbers
  complex(R_KIND), parameter :: c_i = (0.0_R_KIND,1.0_R_KIND)
  complex(R_KIND), parameter :: c_nil = (0.0_R_KIND,0.0_R_KIND)
  complex(R_KIND), parameter :: c_one=(1.0_R_KIND,0.0_R_KIND)
  
  ! other for global variables
  integer :: Ddim
  real(R_KIND) :: Ef, t_glb, mag
  complex(R_KIND), allocatable, dimension(:,:) :: vect, vect0, valt, valt0
  complex(R_KIND), allocatable, dimension(:,:) :: Hglb, Sglb, iglb
  ! for multi-thread integrand 
  complex(R_KIND), allocatable, dimension(:,:) :: zbufs
  real(R_KIND), allocatable, dimension(:) :: sng_pnts
  real(R_KIND), allocatable, dimension(:,:) :: pts_m
  real(R_KIND), allocatable, dimension(:,:) :: alist_m,blist_m,rlist_m,elist_m
  integer, allocatable, dimension(:,:) :: iord_m, level_m, ndin_m
  
  ! filepath
  character(255) :: fi_t0 = "input//"
  character(255) :: fi_hs = "input//"
  character(255) :: fout = "output//"
  
  ! green function
  type :: GFtype
	integer :: tdim, Ldim, Rdim, Ddim                                 ! matrix dimensions for Total/L-lead/R-lead/Device partition
    integer, allocatable, dimension(:) :: numorb                      ! number of orbitals each atom
	real(R_KIND) :: Ef                                                ! Fermi energy
	real(R_KIND), allocatable, dimension(:) :: VL, VR                 ! time series for L/R bias
    real(R_KIND), allocatable, dimension(:) :: Lcup, Rcup             ! capacitance vector for the coupling between device and lead
    real(R_KIND), allocatable, dimension(:,:) :: ccinv                ! inverse of the capacitance-netwwork matrix
    real(R_KIND), allocatable, dimension(:) :: rhot0
    real(R_KIND), allocatable, dimension(:) :: rhoti, rhotj, rhotk
	complex(R_KIND), allocatable, dimension(:,:) :: HL, HR, SL, SR
	complex(R_KIND), allocatable, dimension(:,:) :: GD, QL, QR, QN, Q, rho
	complex(R_KIND), allocatable, dimension(:,:) :: GLs, GRs          ! surface green functions of electrodes 
	complex(R_KIND), allocatable, dimension(:,:) :: SEL, SER, SEN, SE ! self energy for L-lead/R-lead/Nucleus/Total
	complex(R_KIND), allocatable, dimension(:,:) :: HD, SD            ! initial hamiltonian/overalp matrix for device, exlcuding L/R electrodes
	complex(R_KIND), allocatable, dimension(:,:) :: Ht, St            ! initial H/S for total system
	complex(R_KIND), allocatable, dimension(:,:) :: Hi, Si            ! H(t)/S(t)
	complex(R_KIND), allocatable, dimension(:,:) :: LaL, GaL, LaR, GaR! the real part and imaginary part of self-energy  
	complex(R_KIND), allocatable, dimension(:,:) :: GaN               ! the imaginary part of self-energy for nonadiabatic coupling
	complex(R_KIND), allocatable, dimension(:,:) :: KL, KR, KN        ! components of dissipation term
	complex(R_KIND), allocatable, dimension(:,:) :: ULi, URi, UNi     ! the index/exponent of the propagator factor
	complex(R_KIND), allocatable, dimension(:,:) :: SD_rcd, Si_rcd, St_rcd
  end type GFtype
  
  ! Table for Gauss–Legendre quadrature
  real(R_KIND) :: GL1(1)=(/ 0.0_R_KIND /)                                    ! scaled coordination within -1 and 1
  real(R_KIND) :: GL1w(1)=(/ 2.0_R_KIND /)                                   ! weighting function
  real(R_KIND) :: GL2(2)=(/ -sqrt(1.0_R_KIND/3.0_R_KIND),sqrt(1.0_R_KIND/3.0_R_KIND) /) 
  real(R_KIND) :: GL2w(2)=(/ 1.0_R_KIND,1.0_R_KIND /)
  real(R_KIND) :: GL3(3)=(/ -sqrt(3.0_R_KIND/5.0_R_KIND),0.0_R_KIND,sqrt(3.0_R_KIND/5.0_R_KIND) /)
  real(R_KIND) :: GL3w(3)=(/ 5.0_R_KIND/9.0_R_KIND,8.0_R_KIND/9.0_R_KIND,5.0_R_KIND/9.0_R_KIND /)
  real(R_KIND) :: GL4(4)=(/ -sqrt(3.0_R_KIND/7.0_R_KIND+sqrt(6.0_R_KIND/5.0_R_KIND)*2.0_R_KIND/7.0_R_KIND), &
					& -sqrt(3.0_R_KIND/7.0_R_KIND-sqrt(6.0_R_KIND/5.0_R_KIND)*2.0_R_KIND/7.0_R_KIND), &
					& sqrt(3.0_R_KIND/7.0_R_KIND-sqrt(6.0_R_KIND/5.0_R_KIND)*2.0_R_KIND/7.0_R_KIND), &
					& sqrt(3.0_R_KIND/7.0_R_KIND+sqrt(6.0_R_KIND/5.0_R_KIND)*2.0_R_KIND/7.0_R_KIND) /) 
  real(R_KIND) :: GL4w(4)=(/ (18.0_R_KIND-sqrt(30.0_R_KIND))/36.0_R_KIND,(18.0_R_KIND+sqrt(30.0_R_KIND))/36.0_R_KIND, &
					& (18.0_R_KIND+sqrt(30.0_R_KIND))/36.0_R_KIND,(18.0_R_KIND-sqrt(30.0_R_KIND))/36.0_R_KIND /)
  real(R_KIND) :: GL5(5)=(/ -sqrt(5.0_R_KIND/9.0_R_KIND+sqrt(10.0_R_KIND/7.0_R_KIND)*2.0_R_KIND/9.0_R_KIND), &
					& -sqrt(5.0_R_KIND/9.0_R_KIND-sqrt(10.0_R_KIND/7.0_R_KIND)*2.0_R_KIND/9.0_R_KIND), &
					& 0.0_R_KIND,sqrt(5.0_R_KIND/9.0_R_KIND-sqrt(10.0_R_KIND/7.0_R_KIND)*2.0_R_KIND/9.0_R_KIND), &
					& sqrt(5.0_R_KIND/9.0_R_KIND+sqrt(10.0_R_KIND/7.0_R_KIND)*2.0_R_KIND/9.0_R_KIND) /) 
  real(R_KIND) :: GL5w(5)=(/ (322.0_R_KIND-13.0_R_KIND*sqrt(70.0_R_KIND))/900.0_R_KIND, &
					& (322.0_R_KIND+13.0_R_KIND*sqrt(70.0_R_KIND))/900.0_R_KIND, &
					& 128.0_R_KIND/225.0_R_KIND,(322.0_R_KIND+13.0_R_KIND*sqrt(70.0_R_KIND))/900.0_R_KIND, &
					& (322.0_R_KIND-13.0_R_KIND*sqrt(70.0_R_KIND))/900.0_R_KIND /)   
  
  end module globaldef

    
    
    