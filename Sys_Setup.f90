module sys_setup
  use globaldef
  use matrixoperation
  implicit none
  
  contains
  
  
  !....................................................................
  ! read parameters of initial setups
  !....................................................................
  subroutine read_t0(GF)
	type(GFtype) :: GF
	integer :: ierr, ndevatm
	integer :: i, j, k, m, n, o, p, ndim1, ndim2, ndim3
	real(R_KIND) :: dt
	real(R_KIND) :: db1, db2, db3, db4, db5
	real(R_KIND), allocatable, dimension(:,:) :: datain
	character(255) :: fpath, str1, str2, str3, str4, str5

	! set parameters from parameters.txt
	genTf = .false.
	readrho = .false.
	cv_model = .false.
	TdHS =  .false.
	fpath = trim(fi_t0) // "parameters.txt"             ! read parameters
	open(UNIT=431, FILE=trim(fpath), STATUS="OLD", IOSTAT=ierr)
	if(ierr /= 0) stop "sys_setup: parameters.txt is not found in default path!"
	read(431,*) str1
	read(431,*) str1, ncpus
	read(431,*) str1, omega
	read(431,*) str1, ntstep
	read(431,*) str1, ntstep_pre
	read(431,*) str1, dt_std
	read(431,*) str1, dt_pre
	read(431,*) str1, Vb
	read(431,*) str1, VLac
	read(431,*) str1, VRac
	read(431,*) str1, tau
	read(431,*) str1, erng
	read(431,*) str1, deng
	read(431,*) str1, elow
	read(431,*) str1, dee
	read(431,*) str1, m
	read(431,*) str1, n
	read(431,*) str1, p
	read(431,*) str1, o
	read(431,*) str1, diag_int
	! read(431,*) str1, fi_t0
	! read(431,*) str1, fi_hs
	! read(431,*) str1, fout
	if(m==1) genTf = .true.
	if(n==1) readrho = .true.
	if(p==1) TdHS = .true.
	if(o==1) cv_model = .true.
	close(431)
	write(*,*) "set num. of thread:                         ", ncpus
	write(*,*) "angular frequency of AC voltage(1/s):       ", omega
	write(*,*) "time steps for simulation:                  ", ntstep
	write(*,*) "time steps for pre-equilibrium:             ", ntstep_pre
	write(*,*) "time interval (ps) for solving TDNEGF:      ", dt_std
	write(*,*) "time interval (ps) for initial equilibrium: ", dt_pre
	write(*,*) "voltage for symmetric bias (V):             ", Vb
	write(*,*) "amplitude of AC voltage for L lead:         ", VLac
	write(*,*) "amplitude of AC voltage for R lead:         ", VRac
	write(*,*) "rise time for lead bias (ps):               ", tau
	write(*,*) "spectrum range of transmission (eV):        ", erng
	write(*,*) "spectrum interval of transmission (eV):     ", deng
	write(*,*) "low-bound (eV) for energy integration:      ", elow
	write(*,*) "energy interval for energy integration:     ", dee
	write(*,*) "bool for calculating transmission T:         ", genTf
	write(*,*) "bool for reading initial rho:                ", readrho
	write(*,*) "solve time-dependent H/S system:             ", TdHS
	write(*,*) "bool for switch on/off capacitance model:    ", cv_model
	write(*,*) "integral scheme:                             ", diag_int
	if(readrho .and. tau > dt_std) then
		write(*,*) ""
		write(*,*) ">>> Warning!!! if set initial rho, it is suggested for tau = dt_std."
		write(*,*) ""
	end if
	
	
	! set time series for electrode's bias
	dt = dt_std
	allocate(GF%VL(ntstep+1), GF%VR(ntstep+1))
	do i = 0, ntstep
		GF%VL(i+1) = Vb*(1.0_R_KIND-exp(-i*dt/tau))+VLac*sin((omega*10.0_R_KIND**(-12))*i*dt)
		GF%VR(i+1) = -Vb*(1.0_R_KIND-exp(-i*dt/tau))+VRac*sin((omega*10.0_R_KIND**(-12))*i*dt)
	end do

	! read parameters for the capacitance-network model
	if(cv_model) then
		fpath = trim(fi_t0) // "Lcouple.txt"         ! for L-coupled vector
	    open(UNIT=511, FILE=trim(fpath), STATUS="OLD", IOSTAT=ierr)
	    if(ierr /= 0) stop "sys_setup: Lcouple.txt is not found in default path!"
	    read(511,*) ndim1
		allocate(GF%Lcup(ndim1),GF%Rcup(ndim1),GF%ccinv(ndim1,ndim1),GF%numorb(ndim1))
		do i=1,ndim1
			read(511,*) GF%Lcup(i)
		end do
	    close(511)
		fpath = trim(fi_t0) // "Rcouple.txt"         ! for R-coupled vector
	    open(UNIT=411, FILE=trim(fpath), STATUS="OLD", IOSTAT=ierr)
	    if(ierr /= 0) stop "sys_setup: Rcouple.txt is not found in default path!"
	    read(411,*) ndim2
		if(ndim1 .ne. ndim2) stop "dimensions of vector for capacitance-network model is error!"
		do i=1,ndim2
			read(411,*) GF%Rcup(i)
		end do
	    close(411)
		fpath = trim(fi_t0) // "num_orb.txt"         ! for number of orbitals each atom
	    open(UNIT=211, FILE=trim(fpath), STATUS="OLD", IOSTAT=ierr)
	    if(ierr /= 0) stop "sys_setup: num_orb.txt is not found in default path!"
	    read(211,*) ndim2
		if(ndim1 .ne. ndim2) stop "dimensions of vector for orbital number is error!"
		do i=1,ndim2
			read(211,*) GF%numorb(i)
		end do
		close(211)
		fpath = trim(fi_t0) // "inv_c_matrix.txt"         ! for inverse matrix of capacitance network
	    open(UNIT=411, FILE=trim(fpath), STATUS="OLD", IOSTAT=ierr)
	    if(ierr /= 0) stop "sys_setup: inv_c_matrix is not found in default path!"
	    read(411,*) ndim2
		if(ndim1 .ne. ndim2) stop "dimensions of inverse matrix for capacitance network is error!"
		do i=1,ndim2
			read(411,*) GF%ccinv(i,:)
		end do
		close(411)
		! for diagonal rho
		allocate(GF%rhot0(ndim1),GF%rhoti(ndim1),GF%rhotj(ndim1))
		do i=1,ndim1
			GF%rhot0(i)=0.0
			GF%rhoti(i)=0.0
			GF%rhotj(i)=0.0
		end do
	end if
	
	! set Fermi energy
	GF%Ef = 0.0_R_KIND
	fpath = trim(fi_t0) // "Ef_source.txt"         ! read Ef for source(L) electrode
	open(UNIT=111, FILE=trim(fpath), STATUS="OLD", IOSTAT=ierr)
	if(ierr /= 0) stop "sys_setup: Ef_source.txt is not found in default path!"
	read(111,*) str1, str2, db1, str3, db2, str4
	close(111)
	GF%Ef = GF%Ef + db2/3.0_R_KIND
	fpath = trim(fi_t0) // "Ef_drain.txt"          ! read Ef for drain(R) electrode
	open(UNIT=111, FILE=trim(fpath), STATUS="OLD", IOSTAT=ierr)
	if(ierr /= 0) stop "sys_setup: Ef_drain.txt is not found in default path!"
	read(111,*) str1, str2, db1, str3, db2, str4
	close(111)
	GF%Ef = GF%Ef + db2/3.0_R_KIND
	fpath = trim(fi_t0) // "Ef_device.txt"         ! read Ef for total system, including electrodes
	open(UNIT=111, FILE=trim(fpath), STATUS="OLD", IOSTAT=ierr)
	if(ierr /= 0) stop "sys_setup: Ef_device.txt is not found in default path!"
	read(111,*) str1, str2, db1, str3, db2, str4
	close(111)
	!GF%Ef = GF%Ef + db2/3.0_R_KIND                ! use Ef by averaging the device, and the semi-infinte electrodes
	GF%Ef = db2                                    ! use Ef of the total system 
	write(*,*) "Fermi-energy = ", GF%Ef
	
	! set hamiltonian/overlap matrix for source(L) electrode, 2PL-only
	fpath = trim(fi_t0) // "hamsqr1_source.dat"	   ! read hamiltonian matrix
	open(UNIT=112, FILE=trim(fpath), STATUS="OLD", IOSTAT=ierr)
	if(ierr /= 0) stop "sys_setup: hamsqr1_source.dat is not found in default path!"
	read(112,*) GF%Ldim
	allocate(GF%HL(GF%Ldim,GF%Ldim), GF%SL(GF%Ldim,GF%Ldim))
	allocate(datain(1,GF%Ldim))
	do i = 1, GF%Ldim
		read(112,*) datain
		GF%HL(i,1:GF%Ldim)=datain(1,:)
	end do
	close(112)
	deallocate(datain)
	
	!GF%HL = GF%HL*Hat2eV                           ! Hatree -> eV
	fpath = trim(fi_t0) // "oversqr_source.dat"	   ! read overlap matrix
	open(UNIT=113, FILE=trim(fpath), STATUS="OLD", IOSTAT=ierr)
	if(ierr /= 0) stop "sys_setup: oversqr_source.dat is not found in default path!"
	read(113,*) m
	if(m .ne. GF%Ldim) stop "sys_setup: HL and SL have different dimensions!"
	allocate(datain(1,GF%Ldim))
	do i = 1, GF%Ldim
		read(113,*) datain
		GF%SL(i,1:GF%Ldim)=datain(1,:)
	end do
	close(113)
	deallocate(datain)
	
	! set hamiltonian/overlap matrix for drain(R) electrode, 2PL-only
	fpath = trim(fi_t0) // "hamsqr1_drain.dat"	   ! read hamiltonian matrix
	open(UNIT=112, FILE=trim(fpath), STATUS="OLD", IOSTAT=ierr)
	if(ierr /= 0) stop "sys_setup: hamsqr1_drain.dat is not found in default path!"
	read(112,*) GF%Rdim
	allocate(GF%HR(GF%Rdim,GF%Rdim), GF%SR(GF%Rdim,GF%Rdim))
	allocate(datain(1,GF%Rdim))
	do i = 1, GF%Rdim
		read(112,*) datain
		GF%HR(i,1:GF%Rdim)=datain(1,:)
	end do
	close(112)
	deallocate(datain)
	
	!GF%HR = GF%HR*Hat2eV                           ! Hatree -> eV
	fpath = trim(fi_t0) // "oversqr_drain.dat"	   ! read overlap matrix
	open(UNIT=113, FILE=trim(fpath), STATUS="OLD", IOSTAT=ierr)
	if(ierr /= 0) stop "sys_setup: oversqr_drain.dat is not found in default path!"
	read(113,*) m
	if(m .ne. GF%Rdim) stop "sys_setup: HR and SR have different dimensions!"
	allocate(datain(1,GF%Rdim))
	do i = 1, GF%Rdim
		read(113,*) datain
		GF%SR(i,1:GF%Rdim)=datain(1,:)
	end do
	close(113)
	deallocate(datain)
	
	! set hamiltonian/overlap matrix for total system, 1PL-source+device+1PL-drain
	fpath = trim(fi_t0) // "hamsqr1_device_ele.dat"	   ! read hamiltonian matrix
	open(UNIT=112, FILE=trim(fpath), STATUS="OLD", IOSTAT=ierr)
	if(ierr /= 0) stop "sys_setup: hamsqr1_device_ele.dat is not found in default path!"
	read(112,*) GF%tdim
	allocate(GF%Ht(GF%tdim,GF%tdim),GF%St(GF%tdim,GF%tdim),GF%St_rcd(GF%tdim,GF%tdim))
	allocate(datain(1,GF%tdim))
	do i = 1, GF%tdim
		read(112,*) datain
		GF%Ht(i,1:GF%tdim)=datain(1,:)
	end do
	close(112)
	deallocate(datain)
	
	!GF%Ht = GF%Ht*Hat2eV                           ! Hatree -> eV
	fpath = trim(fi_t0) // "oversqr_device_ele.dat"	   ! read overlap matrix
	open(UNIT=113, FILE=trim(fpath), STATUS="OLD", IOSTAT=ierr)
	if(ierr /= 0) stop "sys_setup: oversqr_device_ele.dat is not found in default path!"
	read(113,*) m
	if(m .ne. GF%tdim) stop "sys_setup: Ht and St have different dimensions!"
	allocate(datain(1,GF%tdim))
	do i = 1, GF%tdim
		read(113,*) datain
		GF%St(i,1:GF%tdim)=datain(1,:)
	end do
	close(113)
	deallocate(datain)
	GF%St_rcd = GF%St
	
	! set hamiltonian/overlap matrix for device, excluding L/R electrodes
	GF%Ddim = GF%tdim - GF%Ldim/2 - GF%Rdim/2
	allocate(GF%HD(GF%Ddim,GF%Ddim),GF%SD(GF%Ddim,GF%Ddim),GF%SD_rcd(GF%Ddim,GF%Ddim))
	GF%HD = GF%Ht(GF%Ldim/2+1:(GF%tdim-GF%Rdim/2),GF%Ldim/2+1:(GF%tdim-GF%Rdim/2))
	GF%SD = GF%St(GF%Ldim/2+1:(GF%tdim-GF%Rdim/2),GF%Ldim/2+1:(GF%tdim-GF%Rdim/2))
	GF%SD_rcd = GF%St_rcd(GF%Ldim/2+1:(GF%tdim-GF%Rdim/2),GF%Ldim/2+1:(GF%tdim-GF%Rdim/2))
	
	write(*,*) GF%Ddim,GF%Ldim,GF%Rdim,GF%tdim 
	
	! check h/s is hermitian or not
	if(.not. is_hermitian(GF%Ldim,GF%HL)) stop "sys_setup: HL is not hermitian"
	if(.not. is_hermitian(GF%Ldim,GF%SL)) stop "sys_setup: SL is not hermitian"
	if(.not. is_hermitian(GF%Rdim,GF%HR)) stop "sys_setup: HR is not hermitian"
	if(.not. is_hermitian(GF%Rdim,GF%SR)) stop "sys_setup: SR is not hermitian"
	if(.not. is_hermitian(GF%tdim,GF%Ht)) stop "sys_setup: Ht is not hermitian"
	if(.not. is_hermitian(GF%tdim,GF%St)) stop "sys_setup: St is not hermitian"
	if(.not. is_hermitian(GF%Ddim,GF%HD)) stop "sys_setup: HD is not hermitian"
	if(.not. is_hermitian(GF%Ddim,GF%SD)) stop "sys_setup: SD is not hermitian"
	if(.not. is_hermitian(GF%tdim,GF%St_rcd)) stop "sys_setup: St_rcd is not hermitian"
	if(.not. is_hermitian(GF%Ddim,GF%SD_rcd)) stop "sys_setup: SD_rcd is not hermitian"
	write(*,*) "hermitian-check for inputed H/S matrixes is done! "
	
	! for common cases, input of H/S are symmetry (real and hermitian)
	GF%HL = 0.5_R_KIND*(GF%HL+transpose(GF%HL))
	GF%SL = 0.5_R_KIND*(GF%SL+transpose(GF%SL))
	GF%HR = 0.5_R_KIND*(GF%HR+transpose(GF%HR))
	GF%SR = 0.5_R_KIND*(GF%SR+transpose(GF%SR))
	GF%Ht = 0.5_R_KIND*(GF%Ht+transpose(GF%Ht))
	GF%St = 0.5_R_KIND*(GF%St+transpose(GF%St))
	GF%St_rcd = 0.5_R_KIND*(GF%St_rcd+transpose(GF%St_rcd))
	GF%HD = 0.5_R_KIND*(GF%HD+transpose(GF%HD))
	GF%SD = 0.5_R_KIND*(GF%SD+transpose(GF%SD))
	GF%SD_rcd = 0.5_R_KIND*(GF%SD_rcd+transpose(GF%SD_rcd))

	! initialize array
	allocate(GF%GLs(GF%Ldim/2,GF%Ldim/2),GF%GRs(GF%Rdim/2,GF%Rdim/2))
	allocate(GF%GD(GF%Ddim,GF%Ddim),GF%rho(GF%Ddim,GF%Ddim),GF%Q(GF%Ddim,GF%Ddim))
	allocate(GF%QL(GF%Ddim,GF%Ddim),GF%QR(GF%Ddim,GF%Ddim),GF%QN(GF%Ddim,GF%Ddim))
	allocate(GF%LaL(GF%Ddim,GF%Ddim),GF%GaL(GF%Ddim,GF%Ddim))
	allocate(GF%LaR(GF%Ddim,GF%Ddim),GF%GaR(GF%Ddim,GF%Ddim))
	allocate(GF%Hi(GF%Ddim,GF%Ddim),GF%Si(GF%Ddim,GF%Ddim),GF%Si_rcd(GF%Ddim,GF%Ddim))
	allocate(GF%SEL(GF%Ddim,GF%Ddim),GF%SER(GF%Ddim,GF%Ddim))
	allocate(GF%SEN(GF%Ddim,GF%Ddim),GF%SE(GF%Ddim,GF%Ddim))
	allocate(GF%KL(GF%Ddim,GF%Ddim),GF%KR(GF%Ddim,GF%Ddim))
	allocate(GF%GaN(GF%Ddim,GF%Ddim),GF%KN(GF%Ddim,GF%Ddim))
	allocate(GF%ULi(GF%Ddim,GF%Ddim),GF%URi(GF%Ddim,GF%Ddim),GF%UNi(GF%Ddim,GF%Ddim))
	
	! initialize global variables
	allocate(Hglb(GF%Ddim,GF%Ddim),Sglb(GF%Ddim,GF%Ddim),iglb(GF%Ddim,GF%Ddim))
	forall(i=1:GF%Ddim,j=1:GF%Ddim) iglb(i,j) = c_nil
	forall(i=1:GF%Ddim) iglb(i,i) = c_one
	allocate(vect(GF%Ddim,GF%Ddim),vect0(GF%Ddim,GF%Ddim))
	allocate(valt(GF%Ddim,GF%Ddim),valt0(GF%Ddim,GF%Ddim))
	forall(i=1:GF%Ddim,j=1:GF%Ddim) valt(i,j) = c_nil
	forall(i=1:GF%Ddim,j=1:GF%Ddim) valt0(i,j) = c_nil
	allocate(zbufs(1,ncpus))
	allocate(alist_m(nitgl,ncpus),blist_m(nitgl,ncpus),rlist_m(nitgl,ncpus),elist_m(nitgl,ncpus))
	allocate(iord_m(nitgl,ncpus),level_m(nitgl,ncpus))
	write(*,*) "allocating-memory for green-function formulism is done! "
	
  end subroutine read_t0
  
  
  !....................................................................
  ! free array memory
  !....................................................................
  subroutine free_array(GF)
	type(GFtype) :: GF
	
	deallocate(GF%VL, GF%VR)
	deallocate(GF%HL, GF%SL)
	deallocate(GF%HR, GF%SR)
	deallocate(GF%Ht, GF%St)
	deallocate(GF%HD, GF%SD)
	deallocate(GF%GLs,GF%GRs)
	deallocate(GF%GD,GF%rho,GF%Q)
	deallocate(GF%QL,GF%QR,GF%QN)
	deallocate(GF%Hi,GF%Si)
	deallocate(GF%LaL,GF%GaL,GF%LaR,GF%GaR)
	deallocate(GF%SEL,GF%SER,GF%SEN,GF%SE)
	deallocate(GF%GaN,GF%KN,GF%KL,GF%KR)
	deallocate(GF%ULi,GF%URi,GF%UNi)
	deallocate(Hglb,Sglb,iglb)
	deallocate(GF%SD_rcd,GF%Si_rcd,GF%St_rcd)
	deallocate(vect,vect0,valt,valt0)
	deallocate(zbufs)
	deallocate(alist_m,blist_m,rlist_m,elist_m,iord_m,level_m)
	if(cv_model) then 
		deallocate(GF%Lcup,GF%Rcup,GF%ccinv,GF%numorb)
		deallocate(GF%rhot0,GF%rhoti,GF%rhotj)
	end if
  end subroutine free_array
	
end module sys_setup