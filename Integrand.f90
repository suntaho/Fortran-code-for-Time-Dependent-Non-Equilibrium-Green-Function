!======================================================================
!    define integrand functions
!======================================================================
module integrand_fun
  use globaldef
  use omp_lib
  use matrixoperation
  use gfun
  use Complex_Incomplete_Gamma
  implicit none

  contains


  !....................................................................
  ! calculate time dependent current by nonequilibrium green function
  !....................................................................
  subroutine TDNEGF(GF)
	type(GFtype) :: GF
	integer :: i, j, k, m, n, ne, net, tid, ierr, kpnt, ii, u, uu, ndevatm, pp
	integer :: n_export_orbit_rho
	real(R_KIND) :: ee, tt, IL, IR, IN, timeold, timenow
	real(R_KIND) :: VL, VR, thermaleng, fermi_apx
	real(R_KIND) :: dt, rhonew, rhoold, elowold, elowt, pn, db1, db2, db3, ttbuf
	real(R_KIND), allocatable, dimension(:,:) :: Itdata					! time-dependent current
	real(R_KIND), allocatable, dimension(:) :: eiibuf, eii, deii, eiit, deiit
	complex(R_KIND) :: zbuf1, zbuf2, zbuf3, zbuf4, zbuf5, ttbuf2
	complex(R_KIND) :: zbuf1a, zbuf2a, zbuf3a, zbuf4a
	complex(R_KIND), allocatable, dimension(:) :: eval
	complex(R_KIND), allocatable, dimension(:,:) :: heff, heff0, sm, sm0, smold
	complex(R_KIND), allocatable, dimension(:,:) :: heffbuf, smbuf
	complex(R_KIND), allocatable, dimension(:,:) :: buf1, buf2, buf3, idty
	complex(R_KIND), allocatable, dimension(:,:) :: bufold1, bufold2, bufold3, bufold4
	complex(R_KIND), allocatable, dimension(:,:) :: UL, UR, UN
	complex(R_KIND), allocatable, dimension(:,:) :: KL1, KL2, KR1, KR2, KN1, KN2
	complex(R_KIND), allocatable, dimension(:,:) :: k1, k2, k3, k4, rhobak
	complex(R_KIND), allocatable, dimension(:,:) :: Hnow, Hnxt, Snow, Snxt, Snow_rcd, Snxt_rcd
	complex(R_KIND), allocatable, dimension(:,:) :: Kbuf, Kbuft, vectinv, vect0inv
	complex(R_KIND), allocatable, dimension(:,:,:) :: bufn, bufm, bufn2, bufn3
	complex(R_KIND), allocatable, dimension(:,:,:) :: bufp, bufq
	complex(R_KIND), allocatable, dimension(:,:) :: bugx
	real(R_KIND), allocatable, dimension(:) :: sng_pnti
	character(255) :: fpath, str
	logical :: updateK
	! variables for integral function
	integer :: neval, ier, ieri, last, npts
	integer :: iord(nitgl),  level(nitgl)
	real(R_KIND) :: resultr, resulti, abserr
	real(R_KIND) :: alist(nitgl), blist(nitgl), rlist(nitgl), elist(nitgl)
	integer, allocatable, dimension(:) :: ndin
	real(R_KIND), allocatable, dimension(:) :: points, pts
	! for Gauss–Legendre quadrature
	real(R_KIND) :: ra,rb
	real(R_KIND), allocatable, dimension(:) :: glx, glw

	call mkl_set_dynamic(1)
	call omp_set_nested(1)

	thermaleng = 100.0_R_KIND*kbT
	fermi_apx = 1.0_R_KIND/kbT
	! initialized exp_itp for interpolational exp. function
	n= aint((exp_itp_up-exp_itp_low)/exp_itp_de)+1
	allocate(exp_itp(n))
	do i=1,n
		ee=exp_itp_low+exp_itp_de*(i-1.0_R_KIND)
		exp_itp(i)=exp(ee)
	end do

	n_export_orbit_rho = anint(1.0/dt_std/1000.0)
	write(*,*) "iteration interval to export orbital density: ", n_export_orbit_rho

	ndevatm = size(GF%ccinv,1)
	allocate(Itdata(ntstep+1,3))                                        ! time/I_L/I_R
	allocate(eval(GF%Ddim))
	allocate(heff(GF%Ddim,GF%Ddim),heff0(GF%Ddim,GF%Ddim),sm(GF%Ddim,GF%Ddim),smold(GF%Ddim,GF%Ddim))
	allocate(buf1(GF%Ddim,GF%Ddim),buf2(GF%Ddim,GF%Ddim),buf3(GF%Ddim,GF%Ddim),idty(GF%Ddim,GF%Ddim))
	allocate(UL(GF%Ddim,GF%Ddim),UR(GF%Ddim,GF%Ddim),UN(GF%Ddim,GF%Ddim))
	allocate(heffbuf(GF%Ddim,GF%Ddim),smbuf(GF%Ddim,GF%Ddim))
	allocate(KL1(GF%Ddim,GF%Ddim),KL2(GF%Ddim,GF%Ddim))
	allocate(KR1(GF%Ddim,GF%Ddim),KR2(GF%Ddim,GF%Ddim))
	allocate(KN1(GF%Ddim,GF%Ddim),KN2(GF%Ddim,GF%Ddim))
	allocate(k1(GF%Ddim,GF%Ddim),k2(GF%Ddim,GF%Ddim),k3(GF%Ddim,GF%Ddim),k4(GF%Ddim,GF%Ddim),rhobak(GF%Ddim,GF%Ddim))
	allocate(Hnow(GF%Ddim,GF%Ddim),Hnxt(GF%Ddim,GF%Ddim),Snow(GF%Ddim,GF%Ddim),Snxt(GF%Ddim,GF%Ddim))
	allocate(Snow_rcd(GF%Ddim,GF%Ddim),Snxt_rcd(GF%Ddim,GF%Ddim))
	allocate(sm0(GF%Ddim,GF%Ddim),bugx(GF%Ddim,GF%Ddim))
	allocate(Kbuf(GF%Ddim,GF%Ddim),Kbuft(GF%Ddim,GF%Ddim))
	allocate(bufold1(GF%Ddim,GF%Ddim),bufold2(GF%Ddim,GF%Ddim),bufold3(GF%Ddim,GF%Ddim),bufold4(GF%Ddim,GF%Ddim))
	allocate(bufn(GF%Ddim,GF%Ddim,ncpus),bufm(GF%Ddim,GF%Ddim,ncpus))
	allocate(bufn2(GF%Ddim,GF%Ddim,ncpus))
	allocate(bufn3(GF%Ddim,GF%Ddim,ncpus))
	allocate(bufp(GF%Ddim,GF%Ddim,ncpus),bufq(GF%Ddim,GF%Ddim,ncpus))
	forall(i=1:GF%Ddim,j=1:GF%Ddim) idty(i,j) = c_nil
	forall(i=1:GF%Ddim) idty(i,i) = c_one


	if( (TdHS==.true.) .and. (cv_model==.false.)) then                ! use h/s of device from Gromacs for initial setup, avoid the deviation of QM/MM boundary
		call inputHSt(GF,GF%Ddim,0,.false.)                           ! not orthogonalize h/s here, before orthogonalizing self-energy.
		GF%HD = GF%Hi
		GF%SD = GF%Si
		GF%SD_rcd = GF%Si_rcd
	end if


	mag = 1.0_R_KIND                                                 ! re-scale the imaginary delta energy of retarded/advanced green function
	! check orthbol = .true. for diag_int>0
	if( diag_int>0 .and. (.not. orthbol) ) stop "integrand_fun: set orthbol=.true. for diag_int>0 !"
	if( diag_int==2) write(*,*) "integrand_fun: analytical mode currently supports single-thread only!"

	! read table for Gauss–Legendre quadrature
	allocate(glx(GLpnt),glw(GLpnt))
	do i = 1, GLpnt
		if(GLpnt==1) then
			glx(i) = GL1(i)
			glw(i) = GL1w(i)
		else if (GLpnt==2) then
			glx(i) = GL2(i)
			glw(i) = GL2w(i)
		else if (GLpnt==3) then
			glx(i) = GL3(i)
			glw(i) = GL3w(i)
		else if (GLpnt==4) then
			glx(i) = GL4(i)
			glw(i) = GL4w(i)
		else if (GLpnt==5) then
			glx(i) = GL5(i)
			glw(i) = GL5w(i)
		else
			stop "integrand_fun: error on assigning GLpnt! "
		end if
	end do


	! initialize self energy
	call selfenergyLR(GF,GF%Ef)                                       ! initialize self-energy of L/R leads
	! for orthogonal basis
	if(orthbol) then
		buf1 = GF%SEL                                                 ! left-lead's self-energy
		buf2 = GF%SD
		call orthHS(buf1,buf2,GF%Ddim)
		GF%SEL = buf1
		buf1 = GF%SER                                                 ! righ-lead's self-energy
		buf2 = GF%SD
		call orthHS(buf1,buf2,GF%Ddim)
		GF%SER = buf1
		call orthHS(GF%HD,GF%SD,GF%Ddim)                              ! hamiltonian/overlap matrix
	end if
	! for common cases, real H/S/SE are symmetry
	if(maxval(abs(GF%HD-transpose(GF%HD)))>delta .or. maxval(abs(GF%SD-transpose(GF%SD)))>delta .or. &
		& maxval(abs(GF%SD_rcd-transpose(GF%SD_rcd)))>delta .or. maxval(abs(GF%SEL-transpose(GF%SEL)))>delta &
		& .or. maxval(abs(GF%SER-transpose(GF%SER)))>delta) then
		write (*,*) "Wraning! real-value H/S matrix is not symmetry for present system? "
		write (*,*) "HD: ", maxval(abs(GF%HD-transpose(GF%HD)))
		write (*,*) "SD: ", maxval(abs(GF%SD-transpose(GF%SD)))
		write (*,*) "SD_rcd: ", maxval(abs(GF%SD_rcd-transpose(GF%SD_rcd)))
		write (*,*) "SEL: ", maxval(abs(GF%SEL-transpose(GF%SEL)))
		write (*,*) "SER: ", maxval(abs(GF%SER-transpose(GF%SER)))
	else
		GF%HD = 0.5_R_KIND*(GF%HD+transpose(GF%HD))
	    GF%SD = 0.5_R_KIND*(GF%SD+transpose(GF%SD))
	    GF%SD_rcd = 0.5_R_KIND*(GF%SD_rcd+transpose(GF%SD_rcd))
	    GF%SEL = 0.5_R_KIND*(GF%SEL+transpose(GF%SEL))
	    GF%SER = 0.5_R_KIND*(GF%SER+transpose(GF%SER))
	end if
	! assign variables
	GF%LaL = real(GF%SEL, kind=R_KIND)
	GF%GaL = -1.0_R_KIND*aimag(GF%SEL)
	GF%LaR = real(GF%SER, kind=R_KIND)
	GF%GaR = -1.0_R_KIND*aimag(GF%SER)
	forall(i=1:GF%Ddim,j=1:GF%Ddim) GF%GaN(i,j) = c_nil
	heff = GF%HD+GF%SEL+GF%SER                                        ! no electron-nuclues coupling for initil condition
	sm = GF%SD
	Ef = GF%Ef
	Ddim = GF%Ddim
	! export self-energy
	fpath = trim(fout) // "SEL.txt"
	open(UNIT=199, FILE=trim(fpath), STATUS="REPLACE", IOSTAT=ierr)
	if(ierr /= 0) stop "integrand_fun: openning energy_SEL.txt in default path is error!"
	do i = 1,Ddim
		do j = 1,Ddim
			write(199,'(I5,A,I5,A,F20.15,A,F20.15)') i, tab, j, tab, real(GF%SEL(i,j),kind=R_KIND), tab, aimag(GF%SEL(i,j))
		end do
	end do
	close(199)
	fpath = trim(fout) // "SER.txt"
	open(UNIT=299, FILE=trim(fpath), STATUS="REPLACE", IOSTAT=ierr)
	if(ierr /= 0) stop "integrand_fun: openning energy_SER.txt in default path is error!"
	do i = 1,Ddim
		do j = 1,Ddim
			write(299,'(I5,A,I5,A,F20.15,A,F20.15)') i, tab, j, tab, real(GF%SER(i,j),kind=R_KIND), tab, aimag(GF%SER(i,j))
		end do
	end do
	close(299)

	! analyze eigenvalue
	fpath = trim(fout) // "eigval_hD.txt"                             ! for hD
	open(UNIT=115, FILE=trim(fpath), STATUS="REPLACE", IOSTAT=ierr)
	if(ierr /= 0) stop "integrand_fun: openning eigval_hD.txt in default path is error!"
	heffbuf = GF%HD
	smbuf = GF%SD
	call generalized_eigenvalues(heffbuf, smbuf, eval, buf1, buf2, Ddim)! Note! heffbuf/smbuf are overwritten on exit
	do i = 1, Ddim-1                                                  ! sorting eigenvalue by real part
		do j = i+1, Ddim
			if(real(eval(i))>real(eval(j))) then
				zbuf1 = eval(j)
				eval(j) = eval(i)
				eval(i) = zbuf1
			end if
		end do
	end do
	do i = 1, Ddim
		write(115,*) real(eval(i),kind=R_KIND), tab, aimag(eval(i))
	end do
	close(115)
	fpath = trim(fout) // "eigval_heff.txt"                           ! for heff
	open(UNIT=116, FILE=trim(fpath), STATUS="REPLACE", IOSTAT=ierr)
	if(ierr /= 0) stop "integrand_fun: openning eigval_heff.txt in default path is error!"
	heffbuf = GF%HD+GF%SEL+GF%SER
	smbuf = GF%SD
	call generalized_eigenvalues(heffbuf, smbuf, eval, buf1, buf2, Ddim)    ! Note! heffbuf/smbuf are overwritten on exit
	do i = 1, Ddim-1                                                        ! sorting eigenvalue by real part
		do j = i+1, Ddim
			if(real(eval(i),kind=R_KIND)>real(eval(j),kind=R_KIND)) then
				zbuf1 = eval(j)
				eval(j) = eval(i)
				eval(i) = zbuf1
			end if
		end do
	end do
	do i = 1, Ddim
		write(116,*) real(eval(i),kind=R_KIND), tab, aimag(eval(i))
	end do
	close(116)

	! export fermi distribution function
	if (diag_int==0) then
		elow = GF%Ef - abs(Vb) - abs(VLac) - abs(VRac) - thermaleng
		ehigh = GF%Ef + abs(Vb) + abs(VLac) + abs(VRac) + thermaleng   ! uppper bound for energy integration
	else if(diag_int==1) then
		elow = real(eval(1),kind=R_KIND) - abs(Vb) - abs(VLac) - abs(VRac) - thermaleng
		ehigh = GF%Ef + abs(Vb) + abs(VLac) + abs(VRac) + thermaleng  ! uppper bound for energy integration
	else if(diag_int==2) then
		elow = real(eval(1),kind=R_KIND) - abs(Vb) - abs(VLac) - abs(VRac) - thermaleng
		ehigh = GF%Ef-kbT                                                       ! uppper bound for energy integration
		ehigh_apx = GF%Ef + abs(Vb) + abs(VLac) + abs(VRac) + thermaleng
	end if
	dee = dee
	elow = nint(elow/dee)*dee
	ehigh = nint(ehigh/dee)*dee
	ne = nint(abs(ehigh-elow)/dee)
	fpath = trim(fout) // "fermi_distribution.txt"                    ! for heff
	open(UNIT=112, FILE=trim(fpath), STATUS="REPLACE", IOSTAT=ierr)
	if(ierr /= 0) stop "integrand_fun: openning fermi_distribution.txt in default path is error!"
	do i = 0, ne
		ee = elow+1.0_R_KIND*i*dee
		write(112,*) ee, tab, fermi(ee,Ef), Ef
	end do
	close(112)


	! initialize energy grid
	heffbuf = GF%HD+GF%SEL+GF%SER                                     ! calculate eigenvalue
	smbuf = GF%SD
	call generalized_eigenvalues(heffbuf, smbuf, eval, buf1, buf2, Ddim)! Note! heffbuf/smbuf are overwritten on exit
	vect0 = buf2
	vect0inv = vect0
	call matrixinv(vect0inv,Ddim)
	forall(i=1:Ddim,j=1:Ddim) valt0(i,j) = c_nil
	forall(i=1:Ddim) valt0(i,i) = eval(i)
	do i = 1, Ddim-1                                                  ! sorting eigenvalue by real part
		do j = i+1, Ddim
			if(real(eval(i),kind=R_KIND)>real(eval(j),kind=R_KIND)) then
				zbuf1 = eval(j)
				eval(j) = eval(i)
				eval(i) = zbuf1
			end if
		end do
	end do
	kpnt = 0                                                          ! find eigenvalue within selected energy range
	do i = 1, Ddim
		if(real(eval(i),kind=R_KIND)>=elow .and. real(eval(i),kind=R_KIND)<=ehigh) kpnt = kpnt+1
	end do
	npts = kpnt+2
	allocate(eiibuf(kpnt),pts(npts),ndin(npts),points(npts))
	k = 0
	do i = 1, Ddim
		if(real(eval(i),kind=R_KIND)>=elow .and. real(eval(i),kind=R_KIND)<=ehigh) then
			k = k+1
			eiibuf(k) = real(eval(i), kind=R_KIND)
		end if
	end do
	if(k/=kpnt) stop "integrand_fun: error 1 on generating default energy grid for integration of t0!"
	if(t0_intgl_search .and. diag_int==0) then
		points(1:kpnt) = eiibuf
		deallocate(eiibuf)
		! decide energy grid by math library
		Hglb = GF%HD+GF%SEL+GF%SER                                    ! calculate eigenvalue
		Sglb = GF%SD
		call dqagpe(syst0,elow,ehigh,npts,points,100.0_R_KIND*delta,100.0_R_KIND*delta,nitgl,resultr, &
					& abserr,neval,ier,alist,blist,rlist,elist,pts,iord,level,ndin,last)
		if(ier>0) then
			write(*,*) "integrand_fun: math. function at t0 warning code: ", ier
		else
			write(*,'(A50,F8.4,A,I1)') "integrand_fun: functional integral for rho/ier: ", resultr, "/", ier
		end if
		! sorting alist and blist
		do i = 1,last-1
			do j = i+1,last
				if(alist(i)>alist(j)) then
					db1 = alist(i)
					db2 = blist(i)
					db3 = rlist(i)
					ii = level(i)
					alist(i) = alist(j)
					blist(i) = blist(j)
					rlist(i) = rlist(j)
					level(i) = level(j)
					alist(j) = db1
					blist(j) = db2
					rlist(j) = db3
					level(j) = ii
				end if
			end do
		end do
		! export default energy grid
		fpath = trim(fout) // "energy_grid_t0.txt"
		open(UNIT=199, FILE=trim(fpath), STATUS="REPLACE", IOSTAT=ierr)
		if(ierr /= 0) stop "integrand_fun: openning energy_grid_t0.txt in default path is error!"
		do i = 1, last
			write(199,'(I5,A,F20.15,A,F20.15,A,F20.15,A,I3)') i, tab, alist(i), tab, blist(i), tab, rlist(i), tab, level(i)
		end do
		close(199)
		! assign complete energy grid
		ne = last
		do i = 1, last
			if(level(i) > 0) then
				ne = ne+(-1+2**(level(i)))
			end if
		end do
		ne = ne*GLpnt
		allocate(eii(ne),deii(ne))
		forall(i=1:ne) deii(i) = nil
		k = 0                                                         ! set points ans weighting function for integral
		do i = 1, last
			do j = 1, 2**(level(i))
				ra = alist(i)+(j-1)*(blist(i)-alist(i))/(2**(level(i)))
				rb = alist(i)+(j)*(blist(i)-alist(i))/(2**(level(i)))
				do m = 1, GLpnt
					k = k+1
					eii(k) = glx(m)*(rb-ra)/2.0_R_KIND+(ra+rb)/2.0_R_KIND
					deii(k) = glw(m)*(rb-ra)/2.0_R_KIND
				end do
			end do
		end do
		if(k/=ne) stop "integrand_fun: error 1 on generating energy grid for integration of t0!"
		deallocate(pts,ndin,points)
	else if ( (.not. t0_intgl_search) .and. diag_int==0) then
		ne = kpnt                                                     ! add energy grids
		do i = 1, kpnt-1
			db1 = eiibuf(i+1)-eiibuf(i)
			if(db1>dee) ne = ne+aint(db1/dee)
		end do
		ne = (ne-1)*GLpnt
		allocate(eii(ne),deii(ne))
		forall(i=1:ne) deii(i) = nil
		k=0                                                           ! set points ans weighting function for integral
		ra = eiibuf(1)
		do i =2,kpnt
			db1 = eiibuf(i)-eiibuf(i-1)
			if(db1>dee) then
				do j = 1, aint(db1/dee)+1
					rb = eiibuf(i-1)+db1*j/(aint(db1/dee)+1.0_R_KIND)
					do m = 1, GLpnt
						k = k+1
						eii(k) = glx(m)*(rb-ra)/2.0_R_KIND+(ra+rb)/2.0_R_KIND
						deii(k) = glw(m)*(rb-ra)/2.0_R_KIND
					end do
					ra = rb
				end do
			else
				rb = eiibuf(i)
				do m = 1, GLpnt
					k = k+1
					eii(k) = glx(m)*(rb-ra)/2.0_R_KIND+(ra+rb)/2.0_R_KIND
					deii(k) = glw(m)*(rb-ra)/2.0_R_KIND
				end do
				ra = rb
			end if
		end do
		if(k/=ne) stop "integrand_fun: error 2 on generating energy grid for integration of t0!"
		deallocate(eiibuf)
		deallocate(pts,ndin,points)
	end if
	if(diag_int>0) deallocate(eiibuf,pts,ndin,points)


	if(diag_int==0) then
		! export assigned energy grid
		fpath = trim(fout) // "energy_grid__t0_subdivision.txt"
		open(UNIT=19, FILE=trim(fpath), STATUS="REPLACE", IOSTAT=ierr)
		if(ierr /= 0) stop "integrand_fun: openning energy_grid__t0_subdivision in default path is error!"
		do i = 1, ne
			write(19,'(I5,A,F20.15,A,F20.15)') i, tab, eii(i), tab, deii(i)
		end do
		close(19)
		if(k/=ne) stop "integrand_fun: error on generating assigned energy grid for integration!"
	else
		fpath = trim(fout) // "energy_grid__t0_subdivision.txt"
		open(UNIT=19, FILE=trim(fpath), STATUS="REPLACE", IOSTAT=ierr)
		if(ierr /= 0) stop "integrand_fun: openning energy_grid__t0_subdivision in default path is error!"
		close(19)
	end if


	! guess initial charge density
	!if(diag_int<2) call OMP_SET_NUM_THREADS(ncpus)                                ! call opm for later usage
	call OMP_SET_NUM_THREADS(ncpus)
	if(readrho) then
		! clear old document
		fpath = trim(fout) // "time-current_pre.txt"
		open(UNIT=103, FILE=trim(fpath), STATUS="REPLACE", IOSTAT=ierr)
		if(ierr /= 0) stop "integrand_fun: openning time-current_pre.txt in default path is error!"
		close(103)
		! read rho
		fpath = trim(fi_t0) // "rhoini.txt"
		open(UNIT=173, FILE=trim(fpath), STATUS="OLD", IOSTAT=ierr)
		if(ierr /= 0) stop "integrand_fun: openning rhoini.txt in default path is error!"
		read(173,*) i
		if(i.ne.Ddim) stop "integrand_fun: dimensions of hamiltonian matrix is error!"
		do i = 1, Ddim
			do j = 1, Ddim
				read(173,*) db1, db2
				GF%rho(i,j) = db1*c_one+db2*c_i
			end do
		end do
		close(173)
	else
		if(diag_int==0) then
			! integration
			forall(i=1:GF%Ddim,j=1:GF%Ddim,k=1:ncpus) bufn(i,j,k) = c_nil
			Kbuf = GF%GaL+GF%GaR
			!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ee,tid)
			!$OMP DO
			do i = 1, ne
				tid = OMP_GET_THREAD_NUM()+1
				! calculate retarded green function
				ee = eii(i)
				bufm(:,:,tid) = (ee+c_i*delta*mag)*sm-heff
				call matrixinv(bufm(:,:,tid),Ddim)
				! ====== formula type 1 (disable type 2)
				!bufn(:,:,tid) = bufn(:,:,tid)+fermi(ee,Ef)*aimag(bufm(:,:,tid))*(-2.0_R_KIND*deii(i)/pi)
				! ====== formula type 2 (disable type 1)
				call matrixmul(bufp(:,:,tid), bufm(:,:,tid), Kbuf, Ddim)
				bufq(:,:,tid) = dconjg(transpose(bufm(:,:,tid)))
				call matrixmul(bufm(:,:,tid), bufp(:,:,tid), bufq(:,:,tid), Ddim)
				bufn(:,:,tid) = bufn(:,:,tid)+fermi(ee,Ef)*bufm(:,:,tid)*(2.0_R_KIND*deii(i)/pi)

				! show message
				if(tid .eq. 1) write(*,'(A37,F7.3,A2,A,\)') "initializing charge density...cpu0: ",100.0_R_KIND*i/ne*ncpus," %",creturn
			end do
			!$OMP END DO
			!$OMP END PARALLEL
			! collect data from threads
			forall(i=1:GF%Ddim,j=1:GF%Ddim) GF%rho(i,j) = c_nil
			do i = 1, ncpus
				GF%rho = GF%rho+bufn(:,:,i)
			end do
		else if(diag_int==1) then
			! find number of eval within elow-ehigh
			ii = 0
			do i = 1, Ddim
				if( real(valt0(i,i),kind=R_KIND)>elow .and. real(valt0(i,i),kind=R_KIND)<ehigh) ii=ii+1
			end do
			allocate(sng_pnts(ii+2),ndin_m(ii+2,ncpus),pts_m(ii+2,ncpus))
			ii = 0
			do i = 1, Ddim
				if( real(valt0(i,i),kind=R_KIND)>elow .and. real(valt0(i,i),kind=R_KIND)<ehigh) then
					ii=ii+1
					sng_pnts(ii) = real(valt0(i,i),kind=R_KIND)
				end if
			end do
			! call integrand functions
			forall(i=1:Ddim,j=1:Ddim) Kbuf(i,j) = c_nil
			!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(tid,resultr,resulti,abserr,neval,ier,ieri,last,sng_pnti,pp)
			!$OMP DO
			do i = 1, Ddim
				tid = OMP_GET_THREAD_NUM()+1
				zbufs(1,tid) = valt0(i,i)
				if( (real(zbufs(1,tid),kind=R_KIND)>=elow .and. real(zbufs(1,tid),kind=R_KIND)<=ehigh) .and. (aimag(zbufs(1,tid))>=elow .and. aimag(zbufs(1,tid))<=ehigh) ) then
					pp=4
					allocate(sng_pnti(pp))
					sng_pnti(1)=real(zbufs(1,tid),kind=R_KIND)
					sng_pnti(2)=aimag(zbufs(1,tid))
				else if( (real(zbufs(1,tid),kind=R_KIND)>=elow .and. real(zbufs(1,tid),kind=R_KIND)<=ehigh) .and. (aimag(zbufs(1,tid))<elow .or. aimag(zbufs(1,tid))>ehigh) ) then
					pp=3
					allocate(sng_pnti(pp))
					sng_pnti(1)=real(zbufs(1,tid),kind=R_KIND)
				else if( (real(zbufs(1,tid),kind=R_KIND)<elow .or. real(zbufs(1,tid),kind=R_KIND)>ehigh) .and. (aimag(zbufs(1,tid))>=elow .and. aimag(zbufs(1,tid))<=ehigh) ) then
					pp=3
					allocate(sng_pnti(pp))
					sng_pnti(1)=aimag(zbufs(1,tid))
				else
					pp=2
					allocate(sng_pnti(pp))
				end if
				call dqagpe(integrand1_re,elow,ehigh,pp,sng_pnti,0.01_R_KIND*delta,0.01_R_KIND*delta,nitgl,resultr, &
							& abserr,neval,ier,alist_m(:,tid),blist_m(:,tid),rlist_m(:,tid),elist_m(:,tid), &
							& pts_m(:,tid),iord_m(:,tid),level_m(:,tid),ndin_m(:,tid),last)
				call dqagpe(integrand1_im,elow,ehigh,pp,sng_pnti,0.01_R_KIND*delta,0.01_R_KIND*delta,nitgl,resulti, &
							& abserr,neval,ieri,alist_m(:,tid),blist_m(:,tid),rlist_m(:,tid),elist_m(:,tid), &
							& pts_m(:,tid),iord_m(:,tid),level_m(:,tid),ndin_m(:,tid),last)
				if(ier>0 .and. tid==1) write(*,*) "integrand_fun: warning on integral convergence of RE-rho_ini 1!", ier
				if(ieri>0 .and. tid==1) write(*,*) "integrand_fun: warning on integral convergence of IM-rho_ini 1!", ieri
				Kbuf(i,i) = c_one*resultr+c_i*resulti
				deallocate(sng_pnti)
			end do
			!$OMP END DO
			!$OMP END PARALLEL
			deallocate(sng_pnts,ndin_m,pts_m)
			forall(i=1:Ddim,j=1:Ddim) buf2(j,i) = vect0(j,i)*Kbuf(i,i)
			call matrixmul(buf3,buf2,vect0inv,Ddim)
			GF%rho = aimag(buf3)*(-2.0_R_KIND/pi)
		else if(diag_int==2) then
			forall(i=1:Ddim,j=1:Ddim) Kbuf(i,j) = c_nil
			zbuf1a = exp(fermi_apx*(ehigh-ehigh_apx))
			do i = 1, Ddim
				zbuf1 = ehigh-valt0(i,i)
				zbuf2 = elow-valt0(i,i)
				zbuf3 = ehigh_apx-valt0(i,i)
				Kbuf(i,i) = log(zbuf1)-log(zbuf2)+( cdig_fun(c_nil,fermi_apx*zbuf1)-zbuf1a*cdig_fun(c_nil,fermi_apx*zbuf3) )
			end do
			forall(i=1:Ddim,j=1:Ddim) buf2(j,i) = vect0(j,i)*Kbuf(i,i)
			call matrixmul(buf3,buf2,vect0inv,Ddim)
			GF%rho = aimag(buf3)*(-2.0_R_KIND/pi)
		end if
	end if
	write(*,'(A37,F7.3,A12,F9.3)') "initializing charge density...cpu0: ",100.0_R_KIND," % with rho ", real(ctrace(Ddim,GF%rho), kind=R_KIND)

	! equilibrate charge density deviated from acutal value due to wide-band limit approximation
	! initialize the index/exponent of the propagator factor
	forall(i=1:Ddim,j=1:Ddim) GF%ULi(i,j)=c_nil
	forall(i=1:Ddim,j=1:Ddim) GF%URi(i,j)=c_nil
	forall(i=1:Ddim,j=1:Ddim) GF%UNi(i,j)=c_nil
	! time-loops
	fpath = trim(fout) // "time-current_pre.txt"
	open(UNIT=113, FILE=trim(fpath), STATUS="REPLACE", IOSTAT=ierr)
	if(ierr /= 0) stop "integrand_fun: openning time-current_pre.txt in default path is error!"
	sm0 = GF%SD
	GF%Hi = GF%HD
	GF%Si = GF%SD
	if(BDF) then
		bufold1 = GF%rho
		bufold2 = GF%rho
		bufold3 = GF%rho
		bufold4 = GF%rho
	else
		forall(i=1:Ddim,j=1:Ddim) bufold1(i,j) = c_nil
		forall(i=1:Ddim,j=1:Ddim) bufold2(i,j) = c_nil
		forall(i=1:Ddim,j=1:Ddim) bufold3(i,j) = c_nil
		forall(i=1:Ddim,j=1:Ddim) bufold4(i,j) = c_nil
	end if
	dt = dt_pre
	do n = 0, ntstep_pre
		tt = n*dt
		VL = 0.0_R_KIND
		VR = 0.0_R_KIND

		heff0 = GF%HD
		heff = GF%HD
		smold = GF%SD
		sm = GF%SD

		if(n==0) then
			! update the propagator U
			heff0 = heff0+GF%SEL+GF%SER
			heff = heff+GF%SEL+GF%SER
			call c8mat_expm1 ( Ddim, GF%ULi, UL )
			call c8mat_expm1 ( Ddim, GF%URi, UR )

			! calculate K term
			! ========== general integration for K1 ==========
			if(diag_int==0) then
				forall(i=1:GF%Ddim,j=1:GF%Ddim,k=1:ncpus) bufn(i,j,k) = c_nil
				!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ee,tid)
				!$OMP DO
				do i = 1, ne
					tid = OMP_GET_THREAD_NUM()+1
					! calculate retarded green function
					ee = eii(i)
					bufm(:,:,tid) = (ee+c_i*delta*mag)*sm0-heff0
					call matrixinv(bufm(:,:,tid),Ddim)
					! sum grid data
					bufn(:,:,tid) = bufn(:,:,tid)+bufm(:,:,tid)*(fermi(ee,Ef)*deii(i)*exp(c_i*ee*(tt/(hpa*10.0_R_KIND**12))))
				end do
				!$OMP END DO
				!$OMP END PARALLEL
				forall(i=1:Ddim,j=1:Ddim) Kbuf(i,j) = c_nil           ! collect data
				do i = 1, ncpus
					Kbuf = Kbuf+bufn(:,:,i)
				end do
			else if(diag_int==1) then
				t_glb = tt/(hpa*10.0_R_KIND**12)
				! find number of eval within elow-ehigh
				ii = 0
				do i = 1, Ddim
					if( real(valt0(i,i),kind=R_KIND)>elow .and. real(valt0(i,i),kind=R_KIND)<ehigh) ii=ii+1
				end do
				allocate(sng_pnts(ii+2),ndin_m(ii+2,ncpus),pts_m(ii+2,ncpus))
				ii = 0
				do i = 1, Ddim
					if( real(valt0(i,i),kind=R_KIND)>elow .and. real(valt0(i,i),kind=R_KIND)<ehigh) then
						ii=ii+1
						sng_pnts(ii) = real(valt0(i,i),kind=R_KIND)
					end if
				end do
				! call integrand functions
				forall(i=1:Ddim,j=1:Ddim) Kbuf(i,j) = c_nil
				!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(tid,resultr,resulti,abserr,neval,ier,ieri,last,sng_pnti,pp)
				!$OMP DO
				do i = 1, Ddim
					tid = OMP_GET_THREAD_NUM()+1
					zbufs(1,tid) = valt0(i,i)
					if( (real(zbufs(1,tid),kind=R_KIND)>=elow .and. real(zbufs(1,tid),kind=R_KIND)<=ehigh) .and. (aimag(zbufs(1,tid))>=elow .and. aimag(zbufs(1,tid))<=ehigh) ) then
						pp=4
						allocate(sng_pnti(pp))
						sng_pnti(1)=real(zbufs(1,tid),kind=R_KIND)
						sng_pnti(2)=aimag(zbufs(1,tid))
					else if( (real(zbufs(1,tid),kind=R_KIND)>=elow .and. real(zbufs(1,tid),kind=R_KIND)<=ehigh) .and. (aimag(zbufs(1,tid))<elow .or. aimag(zbufs(1,tid))>ehigh) ) then
						pp=3
						allocate(sng_pnti(pp))
						sng_pnti(1)=real(zbufs(1,tid),kind=R_KIND)
					else if( (real(zbufs(1,tid),kind=R_KIND)<elow .or. real(zbufs(1,tid),kind=R_KIND)>ehigh) .and. (aimag(zbufs(1,tid))>=elow .and. aimag(zbufs(1,tid))<=ehigh) ) then
						pp=3
						allocate(sng_pnti(pp))
						sng_pnti(1)=aimag(zbufs(1,tid))
					else
						pp=2
						allocate(sng_pnti(pp))
					end if
					call dqagpe(integrand2_re,elow,ehigh,pp,sng_pnti,0.01_R_KIND*delta,0.01_R_KIND*delta,nitgl,resultr, &
								& abserr,neval,ier,alist_m(:,tid),blist_m(:,tid),rlist_m(:,tid),elist_m(:,tid), &
								& pts_m(:,tid),iord_m(:,tid),level_m(:,tid),ndin_m(:,tid),last)
					call dqagpe(integrand2_im,elow,ehigh,pp,sng_pnti,0.01_R_KIND*delta,0.01_R_KIND*delta,nitgl,resulti, &
								& abserr,neval,ieri,alist_m(:,tid),blist_m(:,tid),rlist_m(:,tid),elist_m(:,tid), &
								& pts_m(:,tid),iord_m(:,tid),level_m(:,tid),ndin_m(:,tid),last)
					if(ier>0 .and. tid==1) write(*,*) "integrand_fun: warning on integral convergence of RE-rho_ini 2!", ier
					if(ieri>0 .and. tid==1) write(*,*) "integrand_fun: warning on integral convergence of IM-rho_ini 2!", ieri
					Kbuf(i,i) = c_one*resultr+c_i*resulti
					deallocate(sng_pnti)
				end do
				!$OMP END DO
				!$OMP END PARALLEL
				deallocate(sng_pnts,ndin_m,pts_m)
				forall(i=1:Ddim,j=1:Ddim) buf2(j,i) = vect0(j,i)*Kbuf(i,i)
				call matrixmul(Kbuf,buf2,vect0inv,Ddim)
			else if(diag_int==2) then
				ttbuf = tt/(hpa*10.0_R_KIND**12)
				ttbuf2 = ttbuf+c_i*fermi_apx
				forall(i=1:Ddim,j=1:Ddim) Kbuf(i,j) = c_nil
				zbuf1a = exp(fermi_apx*(ehigh-ehigh_apx))
				zbuf2a = exp(c_i*ttbuf*elow)
				zbuf3a = exp(c_i*ttbuf*ehigh)
				zbuf4a = exp(fermi_apx*ehigh)
				!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(zbuf1,zbuf2,zbuf3,zbuf4,zbuf5)
				!$OMP DO
				do i = 1, Ddim
					zbuf1 = ehigh-valt0(i,i)
					zbuf2 = elow-valt0(i,i)
					zbuf3 = ehigh_apx-valt0(i,i)
					if(n>0) then
						zbuf5 = ( log(-c_i*zbuf1*ttbuf2)-log(-c_i*zbuf3*ttbuf2)+log(zbuf3)-log(zbuf1) )
						if(abs(aimag(zbuf5))>5.0e-15) then
							zbuf4 = exp(c_i*(valt0(i,i))*ttbuf2)*zbuf4a*zbuf5
						else
							zbuf4 = c_nil
						end if
						zbuf4 = zbuf4+ zbuf3a*cdig_fun(c_nil,-c_i*zbuf1*ttbuf2)-zbuf3a*zbuf1a*cdig_fun(c_nil,-c_i*zbuf3*ttbuf2)

						zbuf5 = ( log(-c_i*zbuf2*ttbuf)-log(-c_i*zbuf1*ttbuf)+log(zbuf1)-log(zbuf2) )
						if(abs(aimag(zbuf5))>5.0e-15) then
							Kbuf(i,i) = exp(c_i*(valt0(i,i))*ttbuf)*zbuf5
						else
							Kbuf(i,i) = c_nil
						end if
						Kbuf(i,i) = Kbuf(i,i)+zbuf2a*cdig_fun(c_nil,-c_i*zbuf2*ttbuf)-zbuf3a*cdig_fun(c_nil,-c_i*zbuf1*ttbuf)

						Kbuf(i,i) = Kbuf(i,i)+zbuf4
					else
						Kbuf(i,i) = log(zbuf1)-log(zbuf2)+( cdig_fun(c_nil,fermi_apx*zbuf1)-zbuf1a*cdig_fun(c_nil,fermi_apx*zbuf3) )
					end if
				end do
				!$OMP END DO
				!$OMP END PARALLEL
				forall(i=1:Ddim,j=1:Ddim) buf2(j,i) = vect0(j,i)*Kbuf(i,i)
				call matrixmul(Kbuf,buf2,vect0inv,Ddim)
			end if
			call matrixmul(buf1,Kbuf,GF%GaL,Ddim)                     ! KL1 for L electrode
			call matrixmul(buf2,UL,buf1,Ddim)
			KL1 = buf2*(-2.0_R_KIND*c_i/pi)
			call matrixmul(buf1,Kbuf,GF%GaR,Ddim)                     ! KR1 for R electrode
			call matrixmul(buf2,UR,buf1,Ddim)
			KR1 = buf2*(-2.0_R_KIND*c_i/pi)
			! ========== integration for KL2 and KR2 ==========
			if(diag_int==0) then
				forall(i=1:GF%Ddim,j=1:GF%Ddim,k=1:ncpus) bufn(i,j,k) = c_nil
				forall(i=1:GF%Ddim,j=1:GF%Ddim,k=1:ncpus) bufn2(i,j,k) = c_nil
				!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ee,tid)
				!$OMP DO
				do i = 1, ne
					tid = OMP_GET_THREAD_NUM()+1
					! calculate retarded green function
					ee = eii(i)
					bufq(:,:,tid) = (ee+c_i*delta*mag)*sm-heff
					call matrixinv(bufq(:,:,tid),Ddim)
					! for KL2
					bufp(:,:,tid) = idty-UL*exp(c_i*(ee+0.0*VL)*(tt/(hpa*10.0_R_KIND**12)))
					call matrixmul(bufm(:,:,tid),bufp(:,:,tid),bufq(:,:,tid),Ddim)
					bufn(:,:,tid) = bufn(:,:,tid)+bufm(:,:,tid)*(fermi(ee+0.0*VL,Ef)*deii(i))
					! for KR2
					bufp(:,:,tid) = idty-UR*exp(c_i*(ee+0.0*VR)*(tt/(hpa*10.0_R_KIND**12)))
					call matrixmul(bufm(:,:,tid),bufp(:,:,tid),bufq(:,:,tid),Ddim)
					bufn2(:,:,tid) = bufn2(:,:,tid)+bufm(:,:,tid)*(fermi(ee+0.0*VR,Ef)*deii(i))
				end do
				!$OMP END DO
				!$OMP END PARALLEL
				! sum for KL2
				forall(i=1:Ddim,j=1:Ddim) Kbuf(i,j) = c_nil               ! collect data
				do i = 1, ncpus
					Kbuf = Kbuf+bufn(:,:,i)
				end do
				call matrixmul(buf1,Kbuf,GF%GaL,Ddim)
				KL2 = buf1*(-2.0_R_KIND*c_i/pi)
				! sum for KR2
				forall(i=1:Ddim,j=1:Ddim) Kbuf(i,j) = c_nil               ! collect data
				do i = 1, ncpus
					Kbuf = Kbuf+bufn2(:,:,i)
				end do
				call matrixmul(buf1,Kbuf,GF%GaR,Ddim)
				KR2 = buf1*(-2.0_R_KIND*c_i/pi)
			else if(diag_int==1) then
				t_glb = tt/(hpa*10.0_R_KIND**12)
				! find number of eval within elow-ehigh
				ii = 0
				do i = 1, Ddim
					if( real(valt0(i,i),kind=R_KIND)>elow .and. real(valt0(i,i),kind=R_KIND)<ehigh) ii=ii+1
				end do
				allocate(sng_pnts(ii+2),ndin_m(ii+2,ncpus),pts_m(ii+2,ncpus))
				ii = 0
				do i = 1, Ddim
					if( real(valt0(i,i),kind=R_KIND)>elow .and. real(valt0(i,i),kind=R_KIND)<ehigh) then
						ii=ii+1
						sng_pnts(ii) = real(valt0(i,i),kind=R_KIND)
					end if
				end do
				! call integrand functions
				forall(i=1:Ddim,j=1:Ddim) Kbuf(i,j) = c_nil
				forall(i=1:Ddim,j=1:Ddim) Kbuft(i,j) = c_nil
				!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(tid,resultr,resulti,abserr,neval,ier,ieri,last,sng_pnti,pp)
				!$OMP DO
				do i = 1, Ddim
					tid = OMP_GET_THREAD_NUM()+1
					zbufs(1,tid) = valt0(i,i)
					if( (real(zbufs(1,tid),kind=R_KIND)>=elow .and. real(zbufs(1,tid),kind=R_KIND)<=ehigh) .and. (aimag(zbufs(1,tid))>=elow .and. aimag(zbufs(1,tid))<=ehigh) ) then
						pp=4
						allocate(sng_pnti(pp))
						sng_pnti(1)=real(zbufs(1,tid),kind=R_KIND)
						sng_pnti(2)=aimag(zbufs(1,tid))
					else if( (real(zbufs(1,tid),kind=R_KIND)>=elow .and. real(zbufs(1,tid),kind=R_KIND)<=ehigh) .and. (aimag(zbufs(1,tid))<elow .or. aimag(zbufs(1,tid))>ehigh) ) then
						pp=3
						allocate(sng_pnti(pp))
						sng_pnti(1)=real(zbufs(1,tid),kind=R_KIND)
					else if( (real(zbufs(1,tid),kind=R_KIND)<elow .or. real(zbufs(1,tid),kind=R_KIND)>ehigh) .and. (aimag(zbufs(1,tid))>=elow .and. aimag(zbufs(1,tid))<=ehigh) ) then
						pp=3
						allocate(sng_pnti(pp))
						sng_pnti(1)=aimag(zbufs(1,tid))
					else
						pp=2
						allocate(sng_pnti(pp))
					end if
					call dqagpe(integrand1_re,elow,ehigh,pp,sng_pnti,0.01_R_KIND*delta,0.01_R_KIND*delta,nitgl,resultr, &
								& abserr,neval,ier,alist_m(:,tid),blist_m(:,tid),rlist_m(:,tid),elist_m(:,tid), &
								& pts_m(:,tid),iord_m(:,tid),level_m(:,tid),ndin_m(:,tid),last)
					call dqagpe(integrand1_im,elow,ehigh,pp,sng_pnti,0.01_R_KIND*delta,0.01_R_KIND*delta,nitgl,resulti, &
								& abserr,neval,ieri,alist_m(:,tid),blist_m(:,tid),rlist_m(:,tid),elist_m(:,tid), &
								& pts_m(:,tid),iord_m(:,tid),level_m(:,tid),ndin_m(:,tid),last)
					if(ier>0 .and. tid==1) write(*,*) "integrand_fun: warning on integral convergence of RE-rho_ini 3!", ier
					if(ieri>0 .and. tid==1) write(*,*) "integrand_fun: warning on integral convergence of IM-rho_ini 3!", ieri
					Kbuf(i,i) = c_one*resultr+c_i*resulti
					call dqagpe(integrand2_re,elow,ehigh,pp,sng_pnti,0.01_R_KIND*delta,0.01_R_KIND*delta,nitgl,resultr, &
								& abserr,neval,ier,alist_m(:,tid),blist_m(:,tid),rlist_m(:,tid),elist_m(:,tid), &
								& pts_m(:,tid),iord_m(:,tid),level_m(:,tid),ndin_m(:,tid),last)
					call dqagpe(integrand2_im,elow,ehigh,pp,sng_pnti,0.01_R_KIND*delta,0.01_R_KIND*delta,nitgl,resulti, &
								& abserr,neval,ieri,alist_m(:,tid),blist_m(:,tid),rlist_m(:,tid),elist_m(:,tid), &
								& pts_m(:,tid),iord_m(:,tid),level_m(:,tid),ndin_m(:,tid),last)
					if(ier>0 .and. tid==1) write(*,*) "integrand_fun: warning on integral convergence of RE-rho_ini 4!", ier
					if(ieri>0 .and. tid==1) write(*,*) "integrand_fun: warning on integral convergence of IM-rho_ini 4!", ieri
					Kbuft(i,i) = c_one*resultr+c_i*resulti
					deallocate(sng_pnti)
				end do
				!$OMP END DO
				!$OMP END PARALLEL
				deallocate(sng_pnts,ndin_m,pts_m)
				! non-exp(t) term
				forall(i=1:Ddim,j=1:Ddim) buf2(j,i) = vect0(j,i)*Kbuf(i,i)
				call matrixmul(Kbuf,buf2,vect0inv,Ddim)
				! exp(t) term
				forall(i=1:Ddim,j=1:Ddim) buf2(j,i) = vect0(j,i)*Kbuft(i,i)
				call matrixmul(Kbuft,buf2,vect0inv,Ddim)
				! sum for KL2
				call matrixmul(buf1,UL,Kbuft,Ddim)
				buf1 = Kbuf-buf1
				call matrixmul(KL2,buf1,GF%GaL,Ddim)
				KL2 = KL2*(-2.0_R_KIND*c_i/pi)
				! sum for KR2
				call matrixmul(buf1,UR,Kbuft,Ddim)
				buf1 = Kbuf-buf1
				call matrixmul(KR2,buf1,GF%GaR,Ddim)
				KR2 = KR2*(-2.0_R_KIND*c_i/pi)
			else if(diag_int==2) then
				ttbuf = tt/(hpa*10.0_R_KIND**12)
				ttbuf2 = ttbuf+c_i*fermi_apx
				forall(i=1:Ddim,j=1:Ddim) Kbuf(i,j) = c_nil
				forall(i=1:Ddim,j=1:Ddim) Kbuft(i,j) = c_nil
				zbuf1a = exp(fermi_apx*(ehigh-ehigh_apx))
				zbuf2a = exp(c_i*ttbuf*elow)
				zbuf3a = exp(c_i*ttbuf*ehigh)
				zbuf4a = exp(fermi_apx*ehigh)
				!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(zbuf1,zbuf2,zbuf3,zbuf4,zbuf5)
				!$OMP DO
				do i = 1, Ddim
					zbuf1 = ehigh-valt0(i,i)
					zbuf2 = elow-valt0(i,i)
					zbuf3 = ehigh_apx-valt0(i,i)
					! non-exp(t) term
					Kbuf(i,i) = log(zbuf1)-log(zbuf2)+( cdig_fun(c_nil,fermi_apx*zbuf1)-zbuf1a*cdig_fun(c_nil,fermi_apx*zbuf3) )

					! exp(t) term
					if(n>0) then
						zbuf5 = ( log(-c_i*zbuf1*ttbuf2)-log(-c_i*zbuf3*ttbuf2)+log(zbuf3)-log(zbuf1) )
						if(abs(aimag(zbuf5))>5.0e-15) then
							zbuf4 = exp(c_i*(valt0(i,i))*ttbuf2)*zbuf4a*zbuf5
						else
							zbuf4 = c_nil
						end if
						zbuf4 = zbuf4+zbuf3a*cdig_fun(c_nil,-c_i*zbuf1*ttbuf2)-zbuf3a*zbuf1a*cdig_fun(c_nil,-c_i*zbuf3*ttbuf2)

						zbuf5 = (log(-c_i*zbuf2*ttbuf)-log(-c_i*zbuf1*ttbuf)+log(zbuf1)-log(zbuf2))
						if(abs(aimag(zbuf5))>5.0e-15) then
							Kbuft(i,i) = exp(c_i*(valt0(i,i))*ttbuf)*zbuf5
						else
							Kbuft(i,i) = c_nil
						end if
						Kbuft(i,i) = Kbuft(i,i)+zbuf2a*cdig_fun(c_nil,-c_i*zbuf2*ttbuf)-zbuf3a*cdig_fun(c_nil,-c_i*zbuf1*ttbuf)

						Kbuft(i,i) = Kbuft(i,i)+zbuf4
					else
						Kbuft(i,i) = Kbuf(i,i)
					end if
				end do
				!$OMP END DO
				!$OMP END PARALLEL
				! non-exp(t) term
				forall(i=1:Ddim,j=1:Ddim) buf2(j,i) = vect0(j,i)*Kbuf(i,i)
				call matrixmul(Kbuf,buf2,vect0inv,Ddim)
				! exp(t) term
				forall(i=1:Ddim,j=1:Ddim) buf2(j,i) = vect0(j,i)*Kbuft(i,i)
				call matrixmul(Kbuft,buf2,vect0inv,Ddim)
				! sum for KL2
				call matrixmul(buf1,UL,Kbuft,Ddim)
				buf1 = Kbuf-buf1
				call matrixmul(KL2,buf1,GF%GaL,Ddim)
				KL2 = KL2*(-2.0_R_KIND*c_i/pi)
				! sum for KR2
				call matrixmul(buf1,UR,Kbuft,Ddim)
				buf1 = Kbuf-buf1
				call matrixmul(KR2,buf1,GF%GaR,Ddim)
				KR2 = KR2*(-2.0_R_KIND*c_i/pi)
			end if
			! update the total K term
			GF%KL = KL1+KL2+dconjg(transpose(KL1+KL2))
			GF%KR = KR1+KR2+dconjg(transpose(KR1+KR2))
		end if

		! calculate QL/QR, no QN for intial state
		GF%QL = GF%KL                                                 ! QL term
		call matrixmul(buf1,GF%LaL,GF%rho,Ddim)
		call matrixmul(buf2,GF%rho,GF%LaL,Ddim)
		GF%QL = GF%QL+c_i*(buf1-buf2)
		call matrixmul(buf1,GF%GaL,GF%rho,Ddim)
		call matrixmul(buf2,GF%rho,GF%GaL,Ddim)
		GF%QL = GF%QL+(buf1+buf2)
		GF%QR = GF%KR                                                 ! QR term
		call matrixmul(buf1,GF%LaR,GF%rho,Ddim)
		call matrixmul(buf2,GF%rho,GF%LaR,Ddim)
		GF%QR = GF%QR+c_i*(buf1-buf2)
		call matrixmul(buf1,GF%GaR,GF%rho,Ddim)
		call matrixmul(buf2,GF%rho,GF%GaR,Ddim)
		GF%QR = GF%QR+(buf1+buf2)
		! calculate current of L/R electrode, nA  (Note whether factor 2 is double counted in rho )
		IL = real(ctrace(Ddim, GF%QL), kind=R_KIND)*e_over_hpa_pi*pi
		IR = real(ctrace(Ddim, GF%QR), kind=R_KIND)*e_over_hpa_pi*pi
		IN = nil
		write(113,'(I9,A,F20.15,A,F20.15,A,F20.15,A,F20.15,A,F20.9)') n, tab, tt, tab, IL, tab, IR, tab, IN, tab, real(ctrace(Ddim,GF%rho), kind=R_KIND)

		! update rho
		GF%Q = GF%QL+GF%QR
		call matrixmul(buf1,GF%Hi,GF%rho,Ddim)
		call matrixmul(buf2,GF%rho,GF%Hi,Ddim)
		buf3 = -1.0_R_KIND*c_i*(buf1-buf2)-GF%Q
		buf1 = buf3
		if(BDF) then
			if(DForder==1) then                                       ! Backward differentiation method
				GF%rho = GF%rho+buf1*(dt/(hpa*10.0_R_KIND**12))
			else if(DForder==2) then
				GF%rho = bufold1*(4.0_R_KIND/3.0_R_KIND)-bufold2*(1.0_R_KIND/3.0_R_KIND)+ &
						& buf1*(2.0_R_KIND/3.0_R_KIND)*(dt/(hpa*10.0_R_KIND**12))
			else if(DForder==3) then
				GF%rho = bufold1*(18.0_R_KIND/11.0_R_KIND)-bufold2*(9.0_R_KIND/11.0_R_KIND)+ &
						& bufold3*(2.0_R_KIND/11.0_R_KIND)+buf1*(6.0_R_KIND/11.0_R_KIND)*(dt/(hpa*10.0_R_KIND**12))
			else if(DForder==4 .or. DForder==5) then
				GF%rho = bufold1*(48.0_R_KIND/25.0_R_KIND)-bufold2*(36.0_R_KIND/25.0_R_KIND)+ &
						& bufold3*(16.0_R_KIND/25.0_R_KIND)-bufold4*(3.0_R_KIND/25.0_R_KIND)+ &
						& buf1*(12.0_R_KIND/25.0_R_KIND)*(dt/(hpa*10.0_R_KIND**12))
			else
				stop "integrand_fun: error Backward-differentiation order setup!"
			end if
			bufold4 = bufold3
			bufold3 = bufold2
			bufold2 = bufold1
			bufold1 = GF%rho
		else
			if(DForder==1) then                                       ! Adams–Bashforth methods
				buf3 = buf1
			else if(DForder==2) then
				buf3 = buf1*(1.5_R_KIND)-bufold1*(0.5_R_KIND)
			else if(DForder==3) then
				buf3 = buf1*(23.0_R_KIND/12.0_R_KIND)-bufold1*(16.0_R_KIND/12.0_R_KIND)+ &
						& bufold2*(5.0_R_KIND/12.0_R_KIND)
			else if(DForder==4) then
				buf3 = buf1*(55.0_R_KIND/24.0_R_KIND)-bufold1*(59.0_R_KIND/24.0_R_KIND)+ &
						& bufold2*(37.0_R_KIND/24.0_R_KIND)-bufold3*(9.0_R_KIND/24.0_R_KIND)
			else if(DForder==5) then
				buf3 = buf1*(1901.0_R_KIND/720.0_R_KIND)-bufold1*(1387.0_R_KIND/360.0_R_KIND)+ &
						& bufold2*(109.0_R_KIND/30.0_R_KIND)-bufold3*(637.0_R_KIND/360.0_R_KIND)+ &
						& bufold4*(251.0_R_KIND/720.0_R_KIND)
			else
				stop "integrand_fun: error AB order setup!"
			end if
			bufold4 = bufold3
			bufold3 = bufold2
			bufold2 = bufold1
			bufold1 = buf1
			GF%rho = GF%rho+buf3*(dt/(hpa*10.0_R_KIND**12))
		end if

		! show message
		write(*,'(A37,F7.3,A2,A,\)') "equilibrate initial charge density: ",100.0_R_KIND*n/ntstep_pre," %",creturn
	end do
	close(113)
	write(*,*) ""


	! export initial charge density for next usage
	fpath = trim(fout) // "rhoini.txt"
	open(UNIT=99, FILE=trim(fpath), STATUS="REPLACE", IOSTAT=ierr)
	if(ierr /= 0) stop "integrand_fun: openning rhoini.txt in default path is error!"
	write(99,*) Ddim
	do i = 1, Ddim
		do j = 1, Ddim
			write(99,'(F20.15,F20.15)') real(GF%rho(i,j), kind=R_KIND), aimag(GF%rho(i,j))
		end do
	end do
	close(99)


	! initialize the index/exponent of the propagator factor
	forall(i=1:Ddim,j=1:Ddim) GF%ULi(i,j)=c_nil
	forall(i=1:Ddim,j=1:Ddim) GF%URi(i,j)=c_nil
	forall(i=1:Ddim,j=1:Ddim) GF%UNi(i,j)=c_nil
	! time-loops
	fpath = trim(fout) // "time-current.txt"
	open(UNIT=112, FILE=trim(fpath), STATUS="REPLACE", IOSTAT=ierr)
	if(ierr /= 0) stop "integrand_fun: openning time-current.txt in default path is error!"
	close(112)
	fpath = trim(fout) // "orbital_rho.txt"
	open(UNIT=113, FILE=trim(fpath), STATUS="REPLACE", IOSTAT=ierr)
	if(ierr /= 0) stop "integrand_fun: openning orbital_rho.txt in default path is error!"
	close(113)
	! record rhot0
	k=0
	do i=1,size(GF%ccinv,1)
		GF%rhot0(i) = 0.0
		do j=1,GF%numorb(i)
			k=k+1
			GF%rhot0(i)=GF%rhot0(i)+real(GF%rho(k,k))
		end do
	end do

	sm0 = GF%SD
	GF%Hi = GF%HD
	GF%Si = GF%SD
	GF%Si_rcd = GF%SD_rcd
	if(BDF) then
		bufold1 = GF%rho
		bufold2 = GF%rho
		bufold3 = GF%rho
		bufold4 = GF%rho
	else
		forall(i=1:Ddim,j=1:Ddim) bufold1(i,j) = c_nil
		forall(i=1:Ddim,j=1:Ddim) bufold2(i,j) = c_nil
		forall(i=1:Ddim,j=1:Ddim) bufold3(i,j) = c_nil
		forall(i=1:Ddim,j=1:Ddim) bufold4(i,j) = c_nil
	end if
	updateK = .true.
	dt = dt_std
	CALL CPU_TIME(timeold)
	do n = 0, ntstep
		! read hamiltonian/overlap
		if( (TdHS==.true.) .and. (cv_model==.false.) ) then
			if(n==0) then
				smold = GF%SD_rcd
				! read on time h/s matrix
				call  inputHSt(GF,GF%Ddim,n,orthbol)
				! for common case, H/S are symmetry
				GF%Hi = 0.5_R_KIND*(GF%Hi+transpose(GF%Hi))
				GF%Si = 0.5_R_KIND*(GF%Si+transpose(GF%Si))
				GF%Si_rcd = 0.5_R_KIND*(GF%Si_rcd+transpose(GF%Si_rcd))
				! assign variables
				Hnow = GF%Hi
				Snow = GF%Si
				Snow_rcd = GF%Si_rcd
			end if
			if(n>0) then
				smold = Snow_rcd
				! assign variables
				Hnow = Hnxt
				Snow = Snxt
				Snow_rcd = Snxt_rcd
			end if
			heff0 = GF%HD
			! read next time h/s matrix
			call  inputHSt(GF,GF%Ddim,n+1,orthbol)
			Hnxt = 0.5_R_KIND*(GF%Hi+transpose(GF%Hi))
			Snxt = 0.5_R_KIND*(GF%Si+transpose(GF%Si))
			Snxt_rcd = 0.5_R_KIND*(GF%Si_rcd+transpose(GF%Si_rcd))
			GF%Hi = Hnow
			GF%Si = Snow
			GF%Si_rcd = Snow_rcd
		else if( (TdHS==.true.) .and. (cv_model==.true.) ) then
			if(n==0) then
				smold = GF%SD_rcd
				! read on time h/s matrix
				call  inputHSt_cv(GF,GF%Ddim,n,orthbol)
				! for common case, H/S are symmetry
				GF%Hi = 0.5_R_KIND*(GF%Hi+transpose(GF%Hi))
				GF%Si = 0.5_R_KIND*(GF%Si+transpose(GF%Si))
				GF%Si_rcd = 0.5_R_KIND*(GF%Si_rcd+transpose(GF%Si_rcd))
				! assign variables
				Hnow = GF%Hi
				Snow = GF%Si
				Snow_rcd = GF%Si_rcd
			end if
			if(n>0) then
				smold = Snow_rcd
				! assign variables
				Hnow = Hnxt
				Snow = Snxt
				Snow_rcd = Snxt_rcd
			end if
			heff0 = GF%HD
			! read next time h/s matrix
			call  inputHSt_cv_next(GF,GF%Ddim,n+1,orthbol)
			Hnxt = 0.5_R_KIND*(GF%Hi+transpose(GF%Hi))
			Snxt = 0.5_R_KIND*(GF%Si+transpose(GF%Si))
			Snxt_rcd = 0.5_R_KIND*(GF%Si_rcd+transpose(GF%Si_rcd))
			GF%Hi = Hnow
			GF%Si = Snow
			GF%Si_rcd = Snow_rcd
		else if( TdHS==.false. ) then
			! for Runge–Kutta method
			Hnow = GF%HD
			Hnxt = GF%HD
			Snow = GF%SD
			Snxt = GF%SD
			Snow_rcd = GF%SD_rcd
			Snxt_rcd = GF%SD_rcd
			! general variables
			smold = GF%SD_rcd
			heff0 = GF%HD
		end if

		rhobak = GF%rho
		! ! =============== begin Runge–Kutta methods ================
		do u = 1, 4
			if(u .eq. 1) then
				tt = n*dt
				VL = GF%VL(n)
				VR = GF%VR(n)
				heff = Hnow
				sm = Snow
				GF%Hi = heff
				GF%Si = sm
			else if( (u.eq.2) .or. (u.eq.3) ) then
				tt = (n+0.5_R_KIND)*dt
				VL = (GF%VL(n)+GF%VL(n+1))/2.0_R_KIND
				VR = (GF%VR(n)+GF%VR(n+1))/2.0_R_KIND
				heff = (Hnow+Hnxt)/2.0_R_KIND
				sm = (Snow+Snxt)/2.0_R_KIND
				GF%Hi = heff
				GF%Si = sm
			else if(u .eq. 4) then
				tt = (n+1.0_R_KIND)*dt
				VL = GF%VL(n+1)
				VR = GF%VR(n+1)
				heff = Hnxt
				sm = Snxt
				GF%Hi = heff
				GF%Si = sm
			end if
			if(updateK .and. (u.ne.3)) then
				! update the propagator U
				if(TdHS .and. NADi) then
					if(u .eq. 1) then
						buf1 = (Snxt_rcd-smold)/2.0_R_KIND
						buf2 = GF%SD_rcd
						call orthHS(buf1, buf2, Ddim)                 ! orthogonalize self-energy term for non-adiabatic coupling
						GF%GaN = buf1*0.5_R_KIND*fEN/(dt/(hpa*10.0_R_KIND**12))
					else if( (u.eq.2) .or. (u.eq.3) ) then
						buf1 = (Snxt_rcd-Snow_rcd)
						buf2 = GF%SD_rcd
						call orthHS(buf1, buf2, Ddim)                 ! orthogonalize self-energy term for non-adiabatic coupling
						GF%GaN = buf1*0.5_R_KIND*fEN/(dt/(hpa*10.0_R_KIND**12))
					else if(u .eq. 4) then
						! linear extrapolation
						buf1 = 2.0_R_KIND*(Snxt_rcd-Snow_rcd)-(Snxt_rcd-smold)/2.0_R_KIND
						buf2 = GF%SD_rcd
						call orthHS(buf1, buf2, Ddim)                 ! orthogonalize self-energy term for non-adiabatic coupling
						GF%GaN = buf1*0.5_R_KIND*fEN/(dt/(hpa*10.0_R_KIND**12))
					end if
				else
					forall(i=1:Ddim,j=1:Ddim) GF%GaN(i,j)=c_nil
				end if
				heff0 = heff0+GF%SEL+GF%SER-c_i*GF%GaN
				heff = heff+GF%SEL+GF%SER-c_i*GF%GaN
				! integrae the index/exponent of the propagator factor
				if(u .eq. 1) then
					buf1 = GF%ULi
					buf2 = GF%URi
					buf3 = GF%UNi
				else if( (u.eq.2) ) then
					GF%ULi = GF%ULi+(heff-VL*sm)*(-0.5_R_KIND*c_i*dt/(hpa*10.0_R_KIND**12))
					GF%URi = GF%URi+(heff-VR*sm)*(-0.5_R_KIND*c_i*dt/(hpa*10.0_R_KIND**12))
					if(TdHS .and. NADi) GF%UNi = GF%UNi+heff*(-0.5_R_KIND*c_i*dt/(hpa*10.0_R_KIND**12))
					buf1 = GF%ULi
					buf2 = GF%URi
					buf3 = GF%UNi
				else if ( (u.eq.3) ) then
					buf1 = GF%ULi
					buf2 = GF%URi
					buf3 = GF%UNi
				else if(u .eq. 4) then
					GF%ULi = GF%ULi+(heff-VL*sm)*(-0.5_R_KIND*c_i*dt/(hpa*10.0_R_KIND**12))
					GF%URi = GF%URi+(heff-VR*sm)*(-0.5_R_KIND*c_i*dt/(hpa*10.0_R_KIND**12))
					if(TdHS .and. NADi) GF%UNi = GF%UNi+heff*(-0.5_R_KIND*c_i*dt/(hpa*10.0_R_KIND**12))
					buf1 = GF%ULi
					buf2 = GF%URi
					buf3 = GF%UNi
				end if
				if(n==0) then
					call c8mat_expm1 ( Ddim, buf1, UL )
					call c8mat_expm1 ( Ddim, buf2, UR )
					if(TdHS .and. NADi) call c8mat_expm1 ( Ddim, buf3, UN )
				else
					call expmk ( Ddim, buf1, UL )                   ! expmk() is faster than c8mat_expm1() by about 15% for 500*500 array
					call expmk ( Ddim, buf2, UR )                   ! Note! empmk() cannot read null matrix, e.g. matrix at n=0.
					if(TdHS .and. NADi) call expmk ( Ddim, buf3, UN )
				end if
				! initialize energy grid for heff(t)
				heffbuf = heff
				smbuf = sm
				call generalized_eigenvalues(heffbuf, smbuf, eval, buf1, buf2, Ddim)! Note! heffbuf/smbuf are overwritten on exit
				vect = buf2
				vectinv = vect
				call matrixinv(vectinv,Ddim)
				forall(i=1:Ddim,j=1:Ddim) valt(i,j) = c_nil
				forall(i=1:Ddim) valt(i,i) = eval(i)
				do i = 1, Ddim-1                                          ! sorting eigenvalue by real part
					do j = i+1, Ddim
						if(real(eval(i),kind=R_KIND)>real(eval(j),kind=R_KIND)) then
							zbuf1 = eval(j)
							eval(j) = eval(i)
							eval(i) = zbuf1
						end if
					end do
				end do
				elowt = real(eval(1),kind=R_KIND) - abs(Vb) - abs(VLac) - abs(VRac) - thermaleng
				kpnt = 0                                                  ! find eigenvalue within selected energy range
				do i = 1, Ddim
					if(real(eval(i),kind=R_KIND)>(elow-dee) .and. real(eval(i),kind=R_KIND)<(ehigh+dee)) kpnt = kpnt+1
				end do
				npts = kpnt+2
				allocate(eiibuf(kpnt),pts(npts),ndin(npts),points(npts))
				k = 0
				do i = 1, Ddim
					if(real(eval(i),kind=R_KIND)>(elow-dee) .and. real(eval(i),kind=R_KIND)<(ehigh+dee)) then
						k = k+1
						eiibuf(k) = real(eval(i), kind=R_KIND)
					end if
				end do
				if(k/=kpnt) stop "integrand_fun: error 1 on generating energy grid for integration of t!"
				if(TdHS .and. t_intgl_search .and. n>0 .and. diag_int==0) then
					points(1:kpnt) = eiibuf
					deallocate(eiibuf)
					! decide energy grid by math library
					Hglb = GF%Hi+GF%SEL+GF%SER-c_i*GF%GaN
					Sglb = GF%Si
					call dqagpe(syst,elow,ehigh,npts,points,100.0_R_KIND*delta,100.0_R_KIND*delta,nitgl,resultr, &
								& abserr,neval,ier,alist,blist,rlist,elist,pts,iord,level,ndin,last)
					if(ier>0) then
						write(*,*) "integrand_fun: math. function at t warning code: ", ier
					else
						write(*,'(A50,F8.4,A,I1)') "integrand_fun: functional integral for rho/ier: ", resultr, "/", ier
					end if
					! sorting alist and blist
					do i = 1,last-1
						do j = i+1,last
							if(alist(i)>alist(j)) then
								db1 = alist(i)
								db2 = blist(i)
								db3 = rlist(i)
								ii = level(i)
								alist(i) = alist(j)
								blist(i) = blist(j)
								rlist(i) = rlist(j)
								level(i) = level(j)
								alist(j) = db1
								blist(j) = db2
								rlist(j) = db3
								level(j) = ii
							end if
						end do
					end do
					! export default energy grid
					fpath = trim(fout) // "energy_grid_t.txt"
					open(UNIT=199, FILE=trim(fpath), STATUS="REPLACE", IOSTAT=ierr)
					if(ierr /= 0) stop "integrand_fun: openning energy_grid_t.txt in default path is error!"
					do i = 1, last
						write(199,'(I5,A,F20.15,A,F20.15,A,F20.15,A,I3)') i, tab, alist(i), tab, blist(i), tab, rlist(i), tab, level(i)
					end do
					close(199)
					! assign complete energy grid
					net = last
					do i = 1, last
						if(level(i) > 0) then
							net = net+(-1+2**(level(i)))
						end if
					end do
					net = net*GLpnt
					allocate(eiit(net),deiit(net))
					forall(i=1:net) deiit(i) = nil
					k = 0                                                 ! set points ans weighting function for integral
					do i = 1, last
						do j = 1, 2**(level(i))
							ra = alist(i)+(j-1)*(blist(i)-alist(i))/(2**(level(i)))
							rb = alist(i)+(j)*(blist(i)-alist(i))/(2**(level(i)))
							do m = 1, GLpnt
								k = k+1
								eiit(k) = glx(m)*(rb-ra)/2.0_R_KIND+(ra+rb)/2.0_R_KIND
								deiit(k) = glw(m)*(rb-ra)/2.0_R_KIND
							end do
						end do
					end do
					deallocate(points,pts,ndin)
				else if (TdHS .and. (.not. t_intgl_search) .and. n>0 .and. diag_int==0) then
					net = kpnt                                            ! add energy grids
					do i = 1, kpnt-1
						db1 = eiibuf(i+1)-eiibuf(i)
						if(db1>dee) net = net+aint(db1/dee)
					end do
					net = (net-1)*GLpnt
					allocate(eiit(net),deiit(net))
					forall(i=1:net) deiit(i) = 0.0_R_KIND
					k=0
					ra = eiibuf(1)
					do i =2,kpnt
						db1 = eiibuf(i)-eiibuf(i-1)
						if(db1>dee) then
							do j = 1, aint(db1/dee)+1
								rb = eiibuf(i-1)+db1*j/(aint(db1/dee)+1.0_R_KIND)
								do m = 1, GLpnt
									k = k+1
									eiit(k) = glx(m)*(rb-ra)/2.0_R_KIND+(ra+rb)/2.0_R_KIND
									deiit(k) = glw(m)*(rb-ra)/2.0_R_KIND
								end do
								ra = rb
							end do
						else
							rb = eiibuf(i)
							do m = 1, GLpnt
								k = k+1
								eiit(k) = glx(m)*(rb-ra)/2.0_R_KIND+(ra+rb)/2.0_R_KIND
								deiit(k) = glw(m)*(rb-ra)/2.0_R_KIND
							end do
							ra = rb
						end if
					end do
					if(k/=net) stop "integrand_fun: error 3 on generating energy grid for integration of t!"
					deallocate(eiibuf)
					deallocate(pts,ndin,points)
				else if ( ((.not. TdHS) .or. n==0) .and. diag_int==0) then
					net = ne                                              ! add energy grids
					k = net
					allocate(eiit(net),deiit(net))
					eiit = eii
					deiit = deii
					deallocate(eiibuf)
					deallocate(pts,ndin,points)
				end if
				if(diag_int>0) deallocate(eiibuf,pts,ndin,points)

				if(diag_int==0) then
					! export assigned energy grid
					fpath	 = trim(fout) // "energy_grid__t_subdivision.txt"
					open(UNIT=19, FILE=trim(fpath), STATUS="REPLACE", IOSTAT=ierr)
					if(ierr /= 0) stop "integrand_fun: openning energy_grid__t_subdivision in default path is error!"
					do i = 1, net
						write(19,'(I5,A,F20.15,A,F20.15)') i, tab, eiit(i), tab, deiit(i)
					end do
					close(19)
					if(k/=net) stop "integrand_fun: error 4 on generating assigned energy grid for integration of t!"
				else
					fpath	 = trim(fout) // "energy_grid__t_subdivision.txt"
					open(UNIT=19, FILE=trim(fpath), STATUS="REPLACE", IOSTAT=ierr)
					if(ierr /= 0) stop "integrand_fun: openning energy_grid__t_subdivision in default path is error!"
					close(19)
				end if


				! calculate K term
				! ========== general integration for K1 ==========
				if(diag_int==0) then
					forall(i=1:GF%Ddim,j=1:GF%Ddim,k=1:ncpus) bufn(i,j,k) = c_nil
					!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ee,tid)
					!$OMP DO
					do i = 1, ne
						tid = OMP_GET_THREAD_NUM()+1
						! calculate retarded green function
						ee = eii(i)
						bufm(:,:,tid) = (ee+c_i*delta*mag)*sm0-heff0
						call matrixinv(bufm(:,:,tid),Ddim)
						! sum grid data
						bufn(:,:,tid) = bufn(:,:,tid)+bufm(:,:,tid)*(fermi(ee,Ef)*deii(i)*exp(c_i*ee*(tt/(hpa*10.0_R_KIND**12))))
					end do
					!$OMP END DO
					!$OMP END PARALLEL
					forall(i=1:Ddim,j=1:Ddim) Kbuf(i,j) = c_nil                   ! collect data
					do i = 1, ncpus
						Kbuf = Kbuf+bufn(:,:,i)
					end do
				else if(diag_int==1) then
					t_glb = tt/(hpa*10.0_R_KIND**12)
					! find number of eval within elow-ehigh
					ii = 0
					do i = 1, Ddim
						if( real(valt0(i,i),kind=R_KIND)>elow .and. real(valt0(i,i),kind=R_KIND)<ehigh) ii=ii+1
					end do
					allocate(sng_pnts(ii+2),ndin_m(ii+2,ncpus),pts_m(ii+2,ncpus))
					ii = 0
					do i = 1, Ddim
						if( real(valt0(i,i),kind=R_KIND)>elow .and. real(valt0(i,i),kind=R_KIND)<ehigh) then
							ii=ii+1
							sng_pnts(ii) = real(valt0(i,i),kind=R_KIND)
						end if
					end do
					! call integrand functions
					forall(i=1:Ddim,j=1:Ddim) Kbuf(i,j) = c_nil
					!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(tid,resultr,resulti,abserr,neval,ier,ieri,last,sng_pnti,pp)
					!$OMP DO
					do i = 1, Ddim
						tid = OMP_GET_THREAD_NUM()+1
						zbufs(1,tid) = valt0(i,i)
						if( (real(zbufs(1,tid),kind=R_KIND)>=elow .and. real(zbufs(1,tid),kind=R_KIND)<=ehigh) .and. (aimag(zbufs(1,tid))>=elow .and. aimag(zbufs(1,tid))<=ehigh) ) then
							pp=4
							allocate(sng_pnti(pp))
							sng_pnti(1)=real(zbufs(1,tid),kind=R_KIND)
							sng_pnti(2)=aimag(zbufs(1,tid))
						else if( (real(zbufs(1,tid),kind=R_KIND)>=elow .and. real(zbufs(1,tid),kind=R_KIND)<=ehigh) .and. (aimag(zbufs(1,tid))<elow .or. aimag(zbufs(1,tid))>ehigh) ) then
							pp=3
							allocate(sng_pnti(pp))
							sng_pnti(1)=real(zbufs(1,tid),kind=R_KIND)
						else if( (real(zbufs(1,tid),kind=R_KIND)<elow .or. real(zbufs(1,tid),kind=R_KIND)>ehigh) .and. (aimag(zbufs(1,tid))>=elow .and. aimag(zbufs(1,tid))<=ehigh) ) then
							pp=3
							allocate(sng_pnti(pp))
							sng_pnti(1)=aimag(zbufs(1,tid))
						else
							pp=2
							allocate(sng_pnti(pp))
						end if
						call dqagpe(integrand2_re,elow,ehigh,pp,sng_pnti,0.01_R_KIND*delta,0.01_R_KIND*delta,nitgl,resultr, &
									& abserr,neval,ier,alist_m(:,tid),blist_m(:,tid),rlist_m(:,tid),elist_m(:,tid), &
									& pts_m(:,tid),iord_m(:,tid),level_m(:,tid),ndin_m(:,tid),last)
						call dqagpe(integrand2_im,elow,ehigh,pp,sng_pnti,0.01_R_KIND*delta,0.01_R_KIND*delta,nitgl,resulti, &
									& abserr,neval,ieri,alist_m(:,tid),blist_m(:,tid),rlist_m(:,tid),elist_m(:,tid), &
									& pts_m(:,tid),iord_m(:,tid),level_m(:,tid),ndin_m(:,tid),last)
						if(ier>0 .and. tid==1) write(*,*) "integrand_fun: warning on integral convergence of RE-rho_ini 5!", ier
						if(ieri>0 .and. tid==1) write(*,*) "integrand_fun: warning on integral convergence of IM-rho_ini 5!", ieri
						Kbuf(i,i) = c_one*resultr+c_i*resulti
						deallocate(sng_pnti)
					end do
					!$OMP END DO
					!$OMP END PARALLEL
					deallocate(sng_pnts,ndin_m,pts_m)
					forall(i=1:Ddim,j=1:Ddim) buf2(j,i) = vect0(j,i)*Kbuf(i,i)
					call matrixmul(Kbuf,buf2,vect0inv,Ddim)
				else if(diag_int==2) then
					ttbuf = tt/(hpa*10.0_R_KIND**12)
					ttbuf2 = ttbuf+c_i*fermi_apx
					forall(i=1:Ddim,j=1:Ddim) Kbuf(i,j) = c_nil
					zbuf1a = exp(fermi_apx*(ehigh-ehigh_apx))
					zbuf2a = exp(c_i*ttbuf*elow)
					zbuf3a = exp(c_i*ttbuf*ehigh)
					zbuf4a = exp(fermi_apx*ehigh)
					!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(zbuf1,zbuf2,zbuf3,zbuf4,zbuf5)
					!$OMP DO
					do i = 1, Ddim
						zbuf1 = ehigh-valt0(i,i)
						zbuf2 = elow-valt0(i,i)
						zbuf3 = ehigh_apx-valt0(i,i)
						if(n>0 .or. u>1) then
							zbuf5 = ( log(-c_i*zbuf1*ttbuf2)-log(-c_i*zbuf3*ttbuf2)+log(zbuf3)-log(zbuf1))
							if(abs(aimag(zbuf5))>5.0e-15) then
								zbuf4 = exp(c_i*(valt0(i,i))*ttbuf2)*zbuf4a*zbuf5
							else
								zbuf4 = c_nil
							end if
						    zbuf4 = zbuf4+zbuf3a*cdig_fun(c_nil,-c_i*zbuf1*ttbuf2)-zbuf3a*zbuf1a*cdig_fun(c_nil,-c_i*zbuf3*ttbuf2)

							zbuf5 = ( log(-c_i*zbuf2*ttbuf)-log(-c_i*zbuf1*ttbuf)+log(zbuf1)-log(zbuf2))
							if(abs(aimag(zbuf5))>5.0e-15) then
								Kbuf(i,i) = exp(c_i*(valt0(i,i))*ttbuf)*zbuf5
							else
								Kbuf(i,i) = c_nil
							end if
						    Kbuf(i,i) = Kbuf(i,i)+zbuf2a*cdig_fun(c_nil,-c_i*zbuf2*ttbuf)-zbuf3a*cdig_fun(c_nil,-c_i*zbuf1*ttbuf)

							Kbuf(i,i) = Kbuf(i,i)+zbuf4
						else
							Kbuf(i,i) = log(zbuf1)-log(zbuf2)+( cdig_fun(c_nil,fermi_apx*zbuf1)-zbuf1a*cdig_fun(c_nil,fermi_apx*zbuf3) )
						end if
					end do
					!$OMP END DO
					!$OMP END PARALLEL
					forall(i=1:Ddim,j=1:Ddim) buf2(j,i) = vect0(j,i)*Kbuf(i,i)
					call matrixmul(Kbuf,buf2,vect0inv,Ddim)
				end if
				call matrixmul(buf1,Kbuf,GF%GaL,Ddim)                         ! KL1 for L electrode
				call matrixmul(buf2,UL,buf1,Ddim)
				KL1 = buf2*(-2.0_R_KIND*c_i/pi)
				call matrixmul(buf1,Kbuf,GF%GaR,Ddim)                         ! KR1 for R electrode
				call matrixmul(buf2,UR,buf1,Ddim)
				KR1 = buf2*(-2.0_R_KIND*c_i/pi)
				if(TdHS .and. NADi) then
					call matrixmul(buf1,Kbuf,GF%GaN,Ddim)                     ! KN1 for nucleus
					call matrixmul(buf2,UN,buf1,Ddim)
					KN1 = buf2*(-2.0_R_KIND*c_i/pi)
				end if
				! ========== integration for KL2/KR2/KN2 ==========
				if(diag_int==0) then
					forall(i=1:GF%Ddim,j=1:GF%Ddim,k=1:ncpus) bufn(i,j,k) = c_nil
					forall(i=1:GF%Ddim,j=1:GF%Ddim,k=1:ncpus) bufn2(i,j,k) = c_nil
					forall(i=1:GF%Ddim,j=1:GF%Ddim,k=1:ncpus) bufn3(i,j,k) = c_nil
					!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ee,tid)
					!$OMP DO
					do i = 1, net
						tid = OMP_GET_THREAD_NUM()+1
						! calculate retarded green function
						ee = eiit(i)
						! for KL2
						bufq(:,:,tid) = (ee+VL+c_i*delta*mag)*sm-heff
						call matrixinv(bufq(:,:,tid),Ddim)
						bufp(:,:,tid) = idty-UL*exp(c_i*(ee)*(tt/(hpa*10.0_R_KIND**12)))
						call matrixmul(bufm(:,:,tid),bufp(:,:,tid),bufq(:,:,tid),Ddim)
						bufn(:,:,tid) = bufn(:,:,tid)+bufm(:,:,tid)*(fermi(ee,Ef)*deiit(i))
						! for KR2
						bufq(:,:,tid) = (ee+VR+c_i*delta*mag)*sm-heff
						call matrixinv(bufq(:,:,tid),Ddim)
						bufp(:,:,tid) = idty-UR*exp(c_i*(ee)*(tt/(hpa*10.0_R_KIND**12)))
						call matrixmul(bufm(:,:,tid),bufp(:,:,tid),bufq(:,:,tid),Ddim)
						bufn2(:,:,tid) = bufn2(:,:,tid)+bufm(:,:,tid)*(fermi(ee,Ef)*deiit(i))
						if(TdHS .and. NADi) then
							! for KN2
							bufp(:,:,tid) = idty-UN*exp(c_i*ee*(tt/(hpa*10.0_R_KIND**12)))
							call matrixmul(bufm(:,:,tid),bufp(:,:,tid),bufq(:,:,tid),Ddim)
							bufn3(:,:,tid) = bufn3(:,:,tid)+bufm(:,:,tid)*(fermi(ee,Ef)*deiit(i))
						end if
					end do
					!$OMP END DO
					!$OMP END PARALLEL
					! sum for KL2
					forall(i=1:Ddim,j=1:Ddim) Kbuf(i,j) = c_nil                   ! collect data
					do i = 1, ncpus
						Kbuf = Kbuf+bufn(:,:,i)
					end do
					call matrixmul(buf1,Kbuf,GF%GaL,Ddim)
					KL2 = buf1*(-2.0_R_KIND*c_i/pi)
					! sum for KR2
					forall(i=1:Ddim,j=1:Ddim) Kbuf(i,j) = c_nil                   ! collect data
					do i = 1, ncpus
						Kbuf = Kbuf+bufn2(:,:,i)
					end do
					call matrixmul(buf1,Kbuf,GF%GaR,Ddim)
					KR2 = buf1*(-2.0_R_KIND*c_i/pi)
					! sum for KN2
					if(TdHS .and. NADi) then
						forall(i=1:Ddim,j=1:Ddim) Kbuf(i,j) = c_nil                   ! collect data
						do i = 1, ncpus
							Kbuf = Kbuf+bufn3(:,:,i)
						end do
						call matrixmul(buf1,Kbuf,GF%GaN,Ddim)
						KN2 = buf1*(-2.0_R_KIND*c_i/pi)
					end if
				else if(diag_int==1) then
					t_glb = tt/(hpa*10.0_R_KIND**12)
					! === KL2 term ===
					! find number of eval within elow-ehigh
					ii = 0
					do i = 1, Ddim
						if( (real(valt(i,i),kind=R_KIND)-VL)>elow .and. (real(valt(i,i),kind=R_KIND)-VL)<ehigh) ii=ii+1
					end do
					allocate(sng_pnts(ii+2),ndin_m(ii+2,ncpus),pts_m(ii+2,ncpus))
					ii = 0
					do i = 1, Ddim
						if( (real(valt(i,i),kind=R_KIND)-VL)>elow .and. (real(valt(i,i),kind=R_KIND)-VL)<ehigh) then
							ii=ii+1
							sng_pnts(ii) = real(valt(i,i),kind=R_KIND)-VL
						end if
					end do
					! call integrand functions
					forall(i=1:Ddim,j=1:Ddim) Kbuf(i,j) = c_nil
					forall(i=1:Ddim,j=1:Ddim) Kbuft(i,j) = c_nil
					!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(tid,resultr,resulti,abserr,neval,ier,ieri,last,sng_pnti,pp)
					!$OMP DO
					do i = 1, Ddim
						tid = OMP_GET_THREAD_NUM()+1
						zbufs(1,tid) = valt(i,i)-VL
						if( (real(zbufs(1,tid),kind=R_KIND)>=elow .and. real(zbufs(1,tid),kind=R_KIND)<=ehigh) .and. (aimag(zbufs(1,tid))>=elow .and. aimag(zbufs(1,tid))<=ehigh) ) then
							pp=4
							allocate(sng_pnti(pp))
							sng_pnti(1)=real(zbufs(1,tid),kind=R_KIND)
							sng_pnti(2)=aimag(zbufs(1,tid))
						else if( (real(zbufs(1,tid),kind=R_KIND)>=elow .and. real(zbufs(1,tid),kind=R_KIND)<=ehigh) .and. (aimag(zbufs(1,tid))<elow .or. aimag(zbufs(1,tid))>ehigh) ) then
							pp=3
							allocate(sng_pnti(pp))
							sng_pnti(1)=real(zbufs(1,tid),kind=R_KIND)
						else if( (real(zbufs(1,tid),kind=R_KIND)<elow .or. real(zbufs(1,tid),kind=R_KIND)>ehigh) .and. (aimag(zbufs(1,tid))>=elow .and. aimag(zbufs(1,tid))<=ehigh) ) then
							pp=3
							allocate(sng_pnti(pp))
							sng_pnti(1)=aimag(zbufs(1,tid))
						else
							pp=2
							allocate(sng_pnti(pp))
						end if
						call dqagpe(integrand1_re,elow,ehigh,pp,sng_pnti,0.01_R_KIND*delta,0.01_R_KIND*delta,nitgl,resultr, &
									& abserr,neval,ier,alist_m(:,tid),blist_m(:,tid),rlist_m(:,tid),elist_m(:,tid), &
									& pts_m(:,tid),iord_m(:,tid),level_m(:,tid),ndin_m(:,tid),last)
						call dqagpe(integrand1_im,elow,ehigh,pp,sng_pnti,0.01_R_KIND*delta,0.01_R_KIND*delta,nitgl,resulti, &
									& abserr,neval,ieri,alist_m(:,tid),blist_m(:,tid),rlist_m(:,tid),elist_m(:,tid), &
									& pts_m(:,tid),iord_m(:,tid),level_m(:,tid),ndin_m(:,tid),last)
						if(ier>0 .and. tid==1) write(*,*) "integrand_fun: warning on integral convergence of RE-rho_ini 6!", ier
						if(ieri>0 .and. tid==1) write(*,*) "integrand_fun: warning on integral convergence of IM-rho_ini 6!", ieri
						Kbuf(i,i) = c_one*resultr+c_i*resulti
						call dqagpe(integrand2_re,elow,ehigh,pp,sng_pnti,0.01_R_KIND*delta,0.01_R_KIND*delta,nitgl,resultr, &
									& abserr,neval,ier,alist_m(:,tid),blist_m(:,tid),rlist_m(:,tid),elist_m(:,tid), &
									& pts_m(:,tid),iord_m(:,tid),level_m(:,tid),ndin_m(:,tid),last)
						call dqagpe(integrand2_im,elow,ehigh,pp,sng_pnti,0.01_R_KIND*delta,0.01_R_KIND*delta,nitgl,resulti, &
									& abserr,neval,ieri,alist_m(:,tid),blist_m(:,tid),rlist_m(:,tid),elist_m(:,tid), &
									& pts_m(:,tid),iord_m(:,tid),level_m(:,tid),ndin_m(:,tid),last)
						if(ier>0 .and. tid==1) write(*,*) "integrand_fun: warning on integral convergence of RE-rho_ini 7!", ier
						if(ieri>0 .and. tid==1) write(*,*) "integrand_fun: warning on integral convergence of IM-rho_ini 7!", ieri
						Kbuft(i,i) = c_one*resultr+c_i*resulti
						deallocate(sng_pnti)
					end do
					!$OMP END DO
					!$OMP END PARALLEL
					deallocate(sng_pnts,ndin_m,pts_m)
					! non-exp(t) term
					forall(i=1:Ddim,j=1:Ddim) buf2(j,i) = vect(j,i)*Kbuf(i,i)
					call matrixmul(Kbuf,buf2,vectinv,Ddim)
					! exp(t) term
					forall(i=1:Ddim,j=1:Ddim) buf2(j,i) = vect(j,i)*Kbuft(i,i)
					call matrixmul(Kbuft,buf2,vectinv,Ddim)
					! sum
					call matrixmul(buf1,UL,Kbuft,Ddim)
					buf1 = Kbuf-buf1
					call matrixmul(KL2,buf1,GF%GaL,Ddim)
					KL2 = KL2*(-2.0_R_KIND*c_i/pi)
					! === KR2 term ===
					! find number of eval within elow-ehigh
					ii = 0
					do i = 1, Ddim
						if( (real(valt(i,i),kind=R_KIND)-VR)>elow .and. (real(valt(i,i),kind=R_KIND)-VR)<ehigh) ii=ii+1
					end do
					allocate(sng_pnts(ii+2),ndin_m(ii+2,ncpus),pts_m(ii+2,ncpus))
					ii = 0
					do i = 1, Ddim
						if( (real(valt(i,i),kind=R_KIND)-VR)>elow .and. (real(valt(i,i),kind=R_KIND)-VR)<ehigh) then
							ii=ii+1
							sng_pnts(ii) = real(valt(i,i),kind=R_KIND)-VR
						end if
					end do
					! call integrand functions
					forall(i=1:Ddim,j=1:Ddim) Kbuf(i,j) = c_nil
					forall(i=1:Ddim,j=1:Ddim) Kbuft(i,j) = c_nil
					!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(tid,resultr,resulti,abserr,neval,ier,ieri,last,sng_pnti,pp)
					!$OMP DO
					do i = 1, Ddim
						tid = OMP_GET_THREAD_NUM()+1
						zbufs(1,tid) = valt(i,i)-VR
						if( (real(zbufs(1,tid),kind=R_KIND)>=elow .and. real(zbufs(1,tid),kind=R_KIND)<=ehigh) .and. (aimag(zbufs(1,tid))>=elow .and. aimag(zbufs(1,tid))<=ehigh) ) then
							pp=4
							allocate(sng_pnti(pp))
							sng_pnti(1)=real(zbufs(1,tid),kind=R_KIND)
							sng_pnti(2)=aimag(zbufs(1,tid))
						else if( (real(zbufs(1,tid),kind=R_KIND)>=elow .and. real(zbufs(1,tid),kind=R_KIND)<=ehigh) .and. (aimag(zbufs(1,tid))<elow .or. aimag(zbufs(1,tid))>ehigh) ) then
							pp=3
							allocate(sng_pnti(pp))
							sng_pnti(1)=real(zbufs(1,tid),kind=R_KIND)
						else if( (real(zbufs(1,tid),kind=R_KIND)<elow .or. real(zbufs(1,tid),kind=R_KIND)>ehigh) .and. (aimag(zbufs(1,tid))>=elow .and. aimag(zbufs(1,tid))<=ehigh) ) then
							pp=3
							allocate(sng_pnti(pp))
							sng_pnti(1)=aimag(zbufs(1,tid))
						else
							pp=2
							allocate(sng_pnti(pp))
						end if
						call dqagpe(integrand1_re,elow,ehigh,pp,sng_pnti,0.01_R_KIND*delta,0.01_R_KIND*delta,nitgl,resultr, &
									& abserr,neval,ier,alist_m(:,tid),blist_m(:,tid),rlist_m(:,tid),elist_m(:,tid), &
									& pts_m(:,tid),iord_m(:,tid),level_m(:,tid),ndin_m(:,tid),last)
						call dqagpe(integrand1_im,elow,ehigh,pp,sng_pnti,0.01_R_KIND*delta,0.01_R_KIND*delta,nitgl,resulti, &
									& abserr,neval,ieri,alist_m(:,tid),blist_m(:,tid),rlist_m(:,tid),elist_m(:,tid), &
									& pts_m(:,tid),iord_m(:,tid),level_m(:,tid),ndin_m(:,tid),last)
						if(ier>0 .and. tid==1) write(*,*) "integrand_fun: warning on integral convergence of RE-rho_ini 8!", ier
						if(ieri>0 .and. tid==1) write(*,*) "integrand_fun: warning on integral convergence of IM-rho_ini 8!", ieri
						Kbuf(i,i) = c_one*resultr+c_i*resulti
						call dqagpe(integrand2_re,elow,ehigh,pp,sng_pnti,0.01_R_KIND*delta,0.01_R_KIND*delta,nitgl,resultr, &
									& abserr,neval,ier,alist_m(:,tid),blist_m(:,tid),rlist_m(:,tid),elist_m(:,tid), &
									& pts_m(:,tid),iord_m(:,tid),level_m(:,tid),ndin_m(:,tid),last)
						call dqagpe(integrand2_im,elow,ehigh,pp,sng_pnti,0.01_R_KIND*delta,0.01_R_KIND*delta,nitgl,resulti, &
									& abserr,neval,ieri,alist_m(:,tid),blist_m(:,tid),rlist_m(:,tid),elist_m(:,tid), &
									& pts_m(:,tid),iord_m(:,tid),level_m(:,tid),ndin_m(:,tid),last)
						if(ier>0 .and. tid==1) write(*,*) "integrand_fun: warning on integral convergence of RE-rho_ini 9!", ier
						if(ieri>0 .and. tid==1) write(*,*) "integrand_fun: warning on integral convergence of IM-rho_ini 9!", ieri
						Kbuft(i,i) = c_one*resultr+c_i*resulti
						deallocate(sng_pnti)
					end do
					!$OMP END DO
					!$OMP END PARALLEL
					deallocate(sng_pnts,ndin_m,pts_m)
					! non-exp(t) term
					forall(i=1:Ddim,j=1:Ddim) buf2(j,i) = vect(j,i)*Kbuf(i,i)
					call matrixmul(Kbuf,buf2,vectinv,Ddim)
					! exp(t) term
					forall(i=1:Ddim,j=1:Ddim) buf2(j,i) = vect(j,i)*Kbuft(i,i)
					call matrixmul(Kbuft,buf2,vectinv,Ddim)
					! sum
					call matrixmul(buf1,UR,Kbuft,Ddim)
					buf1 = Kbuf-buf1
					call matrixmul(KR2,buf1,GF%GaR,Ddim)
					KR2 = KR2*(-2.0_R_KIND*c_i/pi)
					! === KN2 term ===
					if(TdHS .and. NADi) then
						! find number of eval within elow-ehigh
						ii = 0
						do i = 1, Ddim
							if( (real(valt(i,i),kind=R_KIND))>elow .and. (real(valt(i,i),kind=R_KIND))<ehigh) ii=ii+1
						end do
						allocate(sng_pnts(ii+2),ndin_m(ii+2,ncpus),pts_m(ii+2,ncpus))
						ii = 0
						do i = 1, Ddim
							if( (real(valt(i,i),kind=R_KIND))>elow .and. (real(valt(i,i),kind=R_KIND))<ehigh) then
								ii=ii+1
								sng_pnts(ii) = real(valt(i,i),kind=R_KIND)
							end if
						end do
						! call integrand functions
						forall(i=1:Ddim,j=1:Ddim) Kbuf(i,j) = c_nil
						forall(i=1:Ddim,j=1:Ddim) Kbuft(i,j) = c_nil
						!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(tid,resultr,resulti,abserr,neval,ier,ieri,last,sng_pnti,pp)
						!$OMP DO
						do i = 1, Ddim
							tid = OMP_GET_THREAD_NUM()+1
							zbufs(1,tid) = valt(i,i)
							if( (real(zbufs(1,tid),kind=R_KIND)>=elow .and. real(zbufs(1,tid),kind=R_KIND)<=ehigh) .and. (aimag(zbufs(1,tid))>=elow .and. aimag(zbufs(1,tid))<=ehigh) ) then
								pp=4
								allocate(sng_pnti(pp))
								sng_pnti(1)=real(zbufs(1,tid),kind=R_KIND)
								sng_pnti(2)=aimag(zbufs(1,tid))
							else if( (real(zbufs(1,tid),kind=R_KIND)>=elow .and. real(zbufs(1,tid),kind=R_KIND)<=ehigh) .and. (aimag(zbufs(1,tid))<elow .or. aimag(zbufs(1,tid))>ehigh) ) then
								pp=3
								allocate(sng_pnti(pp))
								sng_pnti(1)=real(zbufs(1,tid),kind=R_KIND)
							else if( (real(zbufs(1,tid),kind=R_KIND)<elow .or. real(zbufs(1,tid),kind=R_KIND)>ehigh) .and. (aimag(zbufs(1,tid))>=elow .and. aimag(zbufs(1,tid))<=ehigh) ) then
								pp=3
								allocate(sng_pnti(pp))
								sng_pnti(1)=aimag(zbufs(1,tid))
							else
								pp=2
								allocate(sng_pnti(pp))
							end if
							call dqagpe(integrand1_re,elow,ehigh,pp,sng_pnti,0.01_R_KIND*delta,0.01_R_KIND*delta,nitgl,resultr, &
										& abserr,neval,ier,alist_m(:,tid),blist_m(:,tid),rlist_m(:,tid),elist_m(:,tid), &
										& pts_m(:,tid),iord_m(:,tid),level_m(:,tid),ndin_m(:,tid),last)
							call dqagpe(integrand1_im,elow,ehigh,pp,sng_pnti,0.01_R_KIND*delta,0.01_R_KIND*delta,nitgl,resulti, &
										& abserr,neval,ieri,alist_m(:,tid),blist_m(:,tid),rlist_m(:,tid),elist_m(:,tid), &
										& pts_m(:,tid),iord_m(:,tid),level_m(:,tid),ndin_m(:,tid),last)
							if(ier>0 .and. tid==1) write(*,*) "integrand_fun: warning on integral convergence of RE-rho_ini 10!", ier
							if(ieri>0 .and. tid==1) write(*,*) "integrand_fun: warning on integral convergence of IM-rho_ini 10!", ieri
							Kbuf(i,i) = c_one*resultr+c_i*resulti
							call dqagpe(integrand2_re,elow,ehigh,pp,sng_pnti,0.01_R_KIND*delta,0.01_R_KIND*delta,nitgl,resultr, &
										& abserr,neval,ier,alist_m(:,tid),blist_m(:,tid),rlist_m(:,tid),elist_m(:,tid), &
										& pts_m(:,tid),iord_m(:,tid),level_m(:,tid),ndin_m(:,tid),last)
							call dqagpe(integrand2_im,elow,ehigh,pp,sng_pnti,0.01_R_KIND*delta,0.01_R_KIND*delta,nitgl,resulti, &
										& abserr,neval,ieri,alist_m(:,tid),blist_m(:,tid),rlist_m(:,tid),elist_m(:,tid), &
										& pts_m(:,tid),iord_m(:,tid),level_m(:,tid),ndin_m(:,tid),last)
							if(ier>0 .and. tid==1) write(*,*) "integrand_fun: warning on integral convergence of RE-rho_ini 11!", ier
							if(ieri>0 .and. tid==1) write(*,*) "integrand_fun: warning on integral convergence of IM-rho_ini 11!", ieri
							Kbuft(i,i) = c_one*resultr+c_i*resulti
							deallocate(sng_pnti)
						end do
						!$OMP END DO
						!$OMP END PARALLEL
						deallocate(sng_pnts,ndin_m,pts_m)
						! non-exp(t) term
						forall(i=1:Ddim,j=1:Ddim) buf2(j,i) = vect(j,i)*Kbuf(i,i)
						call matrixmul(Kbuf,buf2,vectinv,Ddim)
						! exp(t) term
						forall(i=1:Ddim,j=1:Ddim) buf2(j,i) = vect(j,i)*Kbuft(i,i)
						call matrixmul(Kbuft,buf2,vectinv,Ddim)
						! sum
						call matrixmul(buf1,UN,Kbuft,Ddim)
						buf1 = Kbuf-buf1
						call matrixmul(KN2,buf1,GF%GaN,Ddim)
						KN2 = KN2*(-2.0_R_KIND*c_i/pi)
					end if
				else if(diag_int==2) then
					ttbuf = tt/(hpa*10.0_R_KIND**12)
					ttbuf2 = ttbuf+c_i*fermi_apx
					zbuf1a = exp(fermi_apx*(ehigh-ehigh_apx))
					zbuf2a = exp(c_i*ttbuf*elow)
					zbuf3a = exp(c_i*ttbuf*ehigh)
					zbuf4a = exp(fermi_apx*ehigh)
					! === KL2 term ===
					forall(i=1:Ddim,j=1:Ddim) Kbuf(i,j) = c_nil
					forall(i=1:Ddim,j=1:Ddim) Kbuft(i,j) = c_nil
					!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(zbuf1,zbuf2,zbuf3,zbuf4,zbuf5)
					!$OMP DO
					do i = 1, Ddim
						zbuf1 = ehigh-(valt(i,i)-VL)
						zbuf2 = elow-(valt(i,i)-VL)
						zbuf3 = ehigh_apx-(valt(i,i)-VL)
						! non-exp(t) term
						Kbuf(i,i) = log(zbuf1)-log(zbuf2)+( cdig_fun(c_nil,fermi_apx*zbuf1)-zbuf1a*cdig_fun(c_nil,fermi_apx*zbuf3) )
						! exp(t) term
						if(n>0 .or. u>1) then
							zbuf5 = ( log(-c_i*zbuf1*ttbuf2)-log(-c_i*zbuf3*ttbuf2)+log(zbuf3)-log(zbuf1))
							if(abs(aimag(zbuf5))>5.0e-15) then
								zbuf4 = exp(c_i*(valt(i,i)-VL)*ttbuf2)*zbuf4a*zbuf5
							else
								zbuf4 = c_nil
							end if
						    zbuf4 = zbuf4+zbuf3a*cdig_fun(c_nil,-c_i*zbuf1*ttbuf2)-zbuf3a*zbuf1a*cdig_fun(c_nil,-c_i*zbuf3*ttbuf2)

							zbuf5 = ( log(-c_i*zbuf2*ttbuf)-log(-c_i*zbuf1*ttbuf)+log(zbuf1)-log(zbuf2))
							if(abs(aimag(zbuf5))>5.0e-15) then
								Kbuft(i,i) = exp(c_i*(valt(i,i)-VL)*ttbuf)*zbuf5
							else
								Kbuft(i,i) = c_nil
							end if
						    Kbuft(i,i) = Kbuft(i,i)+zbuf2a*cdig_fun(c_nil,-c_i*zbuf2*ttbuf)-zbuf3a*cdig_fun(c_nil,-c_i*zbuf1*ttbuf)

							Kbuft(i,i) = Kbuft(i,i)+zbuf4
						else
							Kbuft(i,i) = Kbuf(i,i)
						end if
					end do
					!$OMP END DO
					!$OMP END PARALLEL
					! non-exp(t) term
					forall(i=1:Ddim,j=1:Ddim) buf2(j,i) = vect(j,i)*Kbuf(i,i)
					call matrixmul(Kbuf,buf2,vectinv,Ddim)
					! exp(t) term
					forall(i=1:Ddim,j=1:Ddim) buf2(j,i) = vect(j,i)*Kbuft(i,i)
					call matrixmul(Kbuft,buf2,vectinv,Ddim)
					! sum
					call matrixmul(buf1,UL,Kbuft,Ddim)
					buf1 = Kbuf-buf1
					call matrixmul(KL2,buf1,GF%GaL,Ddim)
					KL2 = KL2*(-2.0_R_KIND*c_i/pi)
					! === KR2 term ===
					forall(i=1:Ddim,j=1:Ddim) Kbuf(i,j) = c_nil
					forall(i=1:Ddim,j=1:Ddim) Kbuft(i,j) = c_nil
					!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(zbuf1,zbuf2,zbuf3,zbuf4,zbuf5)
					!$OMP DO
					do i = 1, Ddim
						zbuf1 = ehigh-(valt(i,i)-VR)
						zbuf2 = elow-(valt(i,i)-VR)
						zbuf3 = ehigh_apx-(valt(i,i)-VR)
						! non-exp(t) term
						Kbuf(i,i) = log(zbuf1)-log(zbuf2)+( cdig_fun(c_nil,fermi_apx*zbuf1)-zbuf1a*cdig_fun(c_nil,fermi_apx*zbuf3) )
						! exp(t) term
						if(n>0 .or. u>1) then
							zbuf5 = ( log(-c_i*zbuf1*ttbuf2)-log(-c_i*zbuf3*ttbuf2)+log(zbuf3)-log(zbuf1))
							if(abs(aimag(zbuf5))>5.0e-15) then
								zbuf4 = exp(c_i*(valt(i,i)-VR)*ttbuf2)*zbuf4a*zbuf5
							else
								zbuf4 = c_nil
							end if
						    zbuf4 = zbuf4+zbuf3a*cdig_fun(c_nil,-c_i*zbuf1*ttbuf2)-zbuf3a*zbuf1a*cdig_fun(c_nil,-c_i*zbuf3*ttbuf2)

							zbuf5 = ( log(-c_i*zbuf2*ttbuf)-log(-c_i*zbuf1*ttbuf)+log(zbuf1)-log(zbuf2))
							if(abs(aimag(zbuf5))>5.0e-15) then
								Kbuft(i,i) = exp(c_i*(valt(i,i)-VR)*ttbuf)*zbuf5
							else
								Kbuft(i,i) = c_nil
							end if
						    Kbuft(i,i) = Kbuft(i,i)+zbuf2a*cdig_fun(c_nil,-c_i*zbuf2*ttbuf)-zbuf3a*cdig_fun(c_nil,-c_i*zbuf1*ttbuf)

							Kbuft(i,i) = Kbuft(i,i)+zbuf4
						else
							Kbuft(i,i) = Kbuf(i,i)
						end if
					end do
					!$OMP END DO
					!$OMP END PARALLEL
					! non-exp(t) term
					forall(i=1:Ddim,j=1:Ddim) buf2(j,i) = vect(j,i)*Kbuf(i,i)
					call matrixmul(Kbuf,buf2,vectinv,Ddim)
					! exp(t) term
					forall(i=1:Ddim,j=1:Ddim) buf2(j,i) = vect(j,i)*Kbuft(i,i)
					call matrixmul(Kbuft,buf2,vectinv,Ddim)
					! sum
					call matrixmul(buf1,UR,Kbuft,Ddim)
					buf1 = Kbuf-buf1
					call matrixmul(KR2,buf1,GF%GaR,Ddim)
					KR2 = KR2*(-2.0_R_KIND*c_i/pi)
					! === KN2 term ===
					if(TdHS .and. NADi) then
						forall(i=1:Ddim,j=1:Ddim) Kbuf(i,j) = c_nil
						forall(i=1:Ddim,j=1:Ddim) Kbuft(i,j) = c_nil
						!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(zbuf1,zbuf2,zbuf3,zbuf4,zbuf5)
						!$OMP DO
						do i = 1, Ddim
							zbuf1 = ehigh-(valt(i,i))
							zbuf2 = elow-(valt(i,i))
							zbuf3 = ehigh_apx-(valt(i,i))
							! non-exp(t) term
							Kbuf(i,i) = log(zbuf1)-log(zbuf2)+( cdig_fun(c_nil,fermi_apx*zbuf1)-zbuf1a*cdig_fun(c_nil,fermi_apx*zbuf3) )
							! exp(t) term
							if(n>0 .or. u>1) then
								zbuf5 = ( log(-c_i*zbuf1*ttbuf2)-log(-c_i*zbuf3*ttbuf2)+log(zbuf3)-log(zbuf1))
								if(abs(aimag(zbuf5))>5.0e-15) then
									zbuf4 = exp(c_i*(valt(i,i))*ttbuf2)*zbuf4a*zbuf5
								else
									zbuf4 = c_nil
								end if
						        zbuf4 = zbuf4+zbuf3a*cdig_fun(c_nil,-c_i*zbuf1*ttbuf2)-zbuf3a*zbuf1a*cdig_fun(c_nil,-c_i*zbuf3*ttbuf2)

								zbuf5 = ( log(-c_i*zbuf2*ttbuf)-log(-c_i*zbuf1*ttbuf)+log(zbuf1)-log(zbuf2))
								if(abs(aimag(zbuf5))>5.0e-15) then
									Kbuft(i,i) = exp(c_i*(valt(i,i))*ttbuf)*zbuf5
								else
									Kbuft(i,i) = c_nil
								end if
						        Kbuft(i,i) = Kbuft(i,i)+zbuf2a*cdig_fun(c_nil,-c_i*zbuf2*ttbuf)-zbuf3a*cdig_fun(c_nil,-c_i*zbuf1*ttbuf)

								Kbuft(i,i) = Kbuft(i,i)+zbuf4
							else
								Kbuft(i,i) = Kbuf(i,i)
							end if
						end do
						!$OMP END DO
						!$OMP END PARALLEL
						! non-exp(t) term
						forall(i=1:Ddim,j=1:Ddim) buf2(j,i) = vect(j,i)*Kbuf(i,i)
						call matrixmul(Kbuf,buf2,vectinv,Ddim)
						! exp(t) term
						forall(i=1:Ddim,j=1:Ddim) buf2(j,i) = vect(j,i)*Kbuft(i,i)
						call matrixmul(Kbuft,buf2,vectinv,Ddim)
						! sum
						call matrixmul(buf1,UN,Kbuft,Ddim)
						buf1 = Kbuf-buf1
						call matrixmul(KN2,buf1,GF%GaN,Ddim)
						KN2 = KN2*(-2.0_R_KIND*c_i/pi)
					end if
				end if
				! update the total K term
				GF%KL = KL1+KL2+dconjg(transpose(KL1+KL2))
				GF%KR = KR1+KR2+dconjg(transpose(KR1+KR2))
				if(TdHS .and. NADi) GF%KN = KN1+KN2+dconjg(transpose(KN1+KN2))

				if(diag_int==0) deallocate(eiit,deiit)
			end if

			! examine if necessary to update the K terms
			if ((.not. TdHS) .and. (.not. NADi) .and. (abs(Vb)<=delta)) updateK = .false.

			! calculate QL/QR/QN
			GF%QL = GF%KL                                                 ! QL term
			call matrixmul(buf1,GF%LaL,GF%rho,Ddim)
			call matrixmul(buf2,GF%rho,GF%LaL,Ddim)
			GF%QL = GF%QL+c_i*(buf1-buf2)
			call matrixmul(buf1,GF%GaL,GF%rho,Ddim)
			call matrixmul(buf2,GF%rho,GF%GaL,Ddim)
			GF%QL = GF%QL+(buf1+buf2)
			GF%QR = GF%KR                                                 ! QR term
			call matrixmul(buf1,GF%LaR,GF%rho,Ddim)
			call matrixmul(buf2,GF%rho,GF%LaR,Ddim)
			GF%QR = GF%QR+c_i*(buf1-buf2)
			call matrixmul(buf1,GF%GaR,GF%rho,Ddim)
			call matrixmul(buf2,GF%rho,GF%GaR,Ddim)
			GF%QR = GF%QR+(buf1+buf2)
			if(TdHS .and. NADi) then                                      ! QN term
				GF%QN = GF%KN
				call matrixmul(buf1,GF%GaN,GF%rho,Ddim)
				call matrixmul(buf2,GF%rho,GF%GaN,Ddim)
				GF%QN = GF%QN+(buf1+buf2)
			end if

			! calculate current of L/R electrode, nA  (Note whether factor 2 is double counted in rho )
			if(u .eq. 1) then
				IL = -real(ctrace(Ddim, GF%QL), kind=R_KIND)*e_over_hpa_pi*pi
				IR = -real(ctrace(Ddim, GF%QR), kind=R_KIND)*e_over_hpa_pi*pi
				IN = nil
				if(TdHS .and. NADi) IN = real(ctrace(Ddim, GF%QN), kind=R_KIND)*e_over_hpa_pi*pi

				fpath = trim(fout) // "time-current.txt"
				open(UNIT=112, FILE=trim(fpath), ACCESS="APPEND", STATUS="OLD", IOSTAT=ierr)
				if(ierr /= 0) stop "integrand_fun: openning time-current.txt in default path is error!"
				write(112,'(I9,A,F20.15,A,F20.15,A,F20.15,A,F20.15,A,F20.9,A,F20.9)') n, tab, tt, tab, IL, tab, IR, tab, IN, tab, real(ctrace(Ddim,GF%rho), kind=R_KIND), tab, elowt-elow
				close(112)
				! export orbital density
				if(mod(n,n_export_orbit_rho)==0) then
					! diagonal part
					fpath = trim(fout) // "orbital_rho.txt"
				    open(UNIT=113, FILE=trim(fpath), ACCESS="APPEND", STATUS="OLD", IOSTAT=ierr)
				    if(ierr /= 0) stop "integrand_fun: openning orbital_rho.txt in default path is error!"
					write(113,'(I9,A)',advance='no') n,tab
					do uu = 1,Ddim
						write(113,'(F12.7,A)',advance='no') real(GF%rho(uu,uu)),tab
					end do
				    write(113,*) ""
				    close(113)
					! full rho matrix
					fpath = trim(fout) // "rho_matrix.txt"
	                open(UNIT=199, FILE=trim(fpath), STATUS="REPLACE", IOSTAT=ierr)
	                if(ierr /= 0) stop "integrand_fun: openning rho_matrix.txt in default path is error!"
	                write(199,*) Ddim
	                do i = 1, Ddim
		                do j = 1, Ddim
			                write(199,'(F20.15,F20.15)') real(GF%rho(i,j), kind=R_KIND), aimag(GF%rho(i,j))
		                end do
	                end do
	                close(199)
				end if
			end if

			if(TdHS .and. NADi) then
				GF%Q = GF%QL+GF%QR+GF%QN
			else
				GF%Q = GF%QL+GF%QR
			end if
			call matrixmul(buf1,GF%Hi,GF%rho,Ddim)
			call matrixmul(buf2,GF%rho,GF%Hi,Ddim)
			buf3 = -1.0_R_KIND*c_i*(buf1-buf2)-GF%Q

			! record slope function
			if(u .eq. 1) then
				k1 = buf3
				GF%rho = rhobak+k1*(0.5_R_KIND*dt/(hpa*10.0_R_KIND**12))
			else if(u .eq. 2) then
				k2 = buf3
				GF%rho = rhobak+k2*(0.5_R_KIND*dt/(hpa*10.0_R_KIND**12))
			else if(u .eq. 3) then
				k3 = buf3
				GF%rho = rhobak+k3*(dt/(hpa*10.0_R_KIND**12))
			else if(u .eq. 4) then
				k4 = buf3
			end if
		end do
		! ================== end of Runge–Kutta methods ==============


		! update rho
		GF%rho =rhobak
		buf1 = (k1+2.0_R_KIND*k2+2.0_R_KIND*k3+k4)/6.0_R_KIND
		if(BDF) then
			if(DForder==1) then                                       ! Backward differentiation method
				GF%rho = GF%rho+buf1*(dt/(hpa*10.0_R_KIND**12))
			else if(DForder==2) then
				GF%rho = bufold1*(4.0_R_KIND/3.0_R_KIND)-bufold2*(1.0_R_KIND/3.0_R_KIND)+ &
						& buf1*(2.0_R_KIND/3.0_R_KIND)*(dt/(hpa*10.0_R_KIND**12))
			else if(DForder==3) then
				GF%rho = bufold1*(18.0_R_KIND/11.0_R_KIND)-bufold2*(9.0_R_KIND/11.0_R_KIND)+ &
						& bufold3*(2.0_R_KIND/11.0_R_KIND)+buf1*(6.0_R_KIND/11.0_R_KIND)*(dt/(hpa*10.0_R_KIND**12))
			else if(DForder==4 .or. DForder==5) then
				GF%rho = bufold1*(48.0_R_KIND/25.0_R_KIND)-bufold2*(36.0_R_KIND/25.0_R_KIND)+ &
						& bufold3*(16.0_R_KIND/25.0_R_KIND)-bufold4*(3.0_R_KIND/25.0_R_KIND)+ &
						& buf1*(12.0_R_KIND/25.0_R_KIND)*(dt/(hpa*10.0_R_KIND**12))
			else
				stop "integrand_fun: error Backward-differentiation order setup!"
			end if
			bufold4 = bufold3
			bufold3 = bufold2
			bufold2 = bufold1
			bufold1 = GF%rho
		else
			if(DForder==1) then                                       ! Adams–Bashforth methods
				buf3 = buf1
			else if(DForder==2) then
				buf3 = buf1*(1.5_R_KIND)-bufold1*(0.5_R_KIND)
			else if(DForder==3) then
				buf3 = buf1*(23.0_R_KIND/12.0_R_KIND)-bufold1*(16.0_R_KIND/12.0_R_KIND)+bufold2*(5.0_R_KIND/12.0_R_KIND)
			else if(DForder==4) then
				buf3 = buf1*(55.0_R_KIND/24.0_R_KIND)-bufold1*(59.0_R_KIND/24.0_R_KIND)+ &
						& bufold2*(37.0_R_KIND/24.0_R_KIND)-bufold3*(9.0_R_KIND/24.0_R_KIND)
			else if(DForder==5) then
				buf3 = buf1*(1901.0_R_KIND/720.0_R_KIND)-bufold1*(1387.0_R_KIND/360.0_R_KIND)+ &
						& bufold2*(109.0_R_KIND/30.0_R_KIND)-bufold3*(637.0_R_KIND/360.0_R_KIND)+ &
						& bufold4*(251.0_R_KIND/720.0_R_KIND)
			else
				stop "integrand_fun: error AB order setup!"
			end if
			bufold4 = bufold3
			bufold3 = bufold2
			bufold2 = bufold1
			bufold1 = buf1
			GF%rho = GF%rho+buf3*(dt/(hpa*10.0_R_KIND**12))
		end if


		! show message
		CALL CPU_TIME(timenow)
		write(*,'(A30,F7.3,A4,F7.3,A,\)') "solving time-dependent NEGF: ",100.0_R_KIND*n/ntstep," %, ",(timenow-timeold)/ncpus,creturn
		timeold=timenow
	end do
	write(*,*) ""


	! free arrays
	deallocate(Itdata,exp_itp)
	deallocate(heff,heff0,idty,sm,sm0,smold,bugx)
	deallocate(bufold1,bufold2,bufold3,bufold4)
	deallocate(buf1,buf2,buf3,bufn,bufm,bufn2,bufn3,bufp,bufq)
	deallocate(KL1,KL2,KR1,KR2,KN1,KN2,Kbuf,Kbuft)
	deallocate(UL,UR,UN,eval,heffbuf,smbuf,glx,glw)
	deallocate(k1,k2,k3,k4,rhobak)
	deallocate(Hnow,Hnxt,Snow,Snxt,Snow_rcd,Snxt_rcd)
	if(diag_int==0) deallocate(eii,deii)

  end subroutine TDNEGF


  !....................................................................
  ! fermi distribution function
  !....................................................................
  real(R_KIND) function fermi(ex,ef)
		integer :: idx
		real(R_KIND) :: ex, ef, db1
		db1=(ex-ef)/kbT
		if(db1<exp_itp_low) then
			fermi = 1.0_R_KIND
		else if (db1>exp_itp_up) then
		    fermi = 0.0_R_KIND
		else
			db1 = (db1-exp_itp_low)/exp_itp_de
			idx=aint(db1)
			db1=db1-idx
			fermi = exp_itp(idx+1)*(1.0_R_KIND-db1)+exp_itp(idx+2)*db1
			fermi = 1.0_R_KIND/(1.0_R_KIND+fermi)
		end if
  end function fermi


  !....................................................................
  ! read hamiltonian/overlap matrixes h(t)/s(t) from file
  !....................................................................
  subroutine inputHSt(GF,ndim,tstep,orthbol_sw)
	type(GFtype) :: GF
	integer :: i, ndim, tstep,ierr
	real(R_KIND), allocatable, dimension(:,:) :: datain
	character(255) :: fpath,str, strdir
	logical :: orthbol_sw

	write(str,*) mod(tstep,mdtstep)+1
	i = anint(1.0_R_KIND*(tstep-mod(tstep,mdtstep))/mdtstep)
	write(strdir,*) i+1

	! set hamiltonian/overlap matrix at tstep
	! h matrix
	write(fpath,*) trim(fi_hs),trim(adjustl(strdir)),trim("/hamil_step_"),trim(adjustl(str))
	open(UNIT=112, FILE=trim(fpath), STATUS="OLD", IOSTAT=ierr)
	do while(ierr /= 0)
		write(*,*) "integrand_fun: ", trim(fpath), " is not found!"
		pause
		open(UNIT=112, FILE=trim(fpath), STATUS="OLD", IOSTAT=ierr)
	end do
	allocate(datain(1,ndim))
	do i = 1, ndim
		read(112,*) datain
		GF%Hi(i,1:ndim)=datain(1,:)
	end do
	close(112)
	deallocate(datain)

	!GF%Hi = GF%Hi*Hat2eV                           ! Hatree -> eV
	write(*,*) "integrand_fun: ", trim(fpath), " was read!"
	! s matrix
	write(fpath,*) trim(fi_hs),trim(adjustl(strdir)),trim("/overl_step_"),trim(adjustl(str))
	open(UNIT=113, FILE=trim(fpath), STATUS="OLD", IOSTAT=ierr)
	do while(ierr /= 0)
		write(*,*) "integrand_fun: ", trim(fpath), " is not found!"
		pause
		open(UNIT=113, FILE=trim(fpath), STATUS="OLD", IOSTAT=ierr)
	end do
	allocate(datain(1,ndim))
	do i = 1, ndim
		read(113,*) datain
		GF%Si(i,1:ndim)=datain(1,:)
	end do
	close(113)
	deallocate(datain)

	GF%Si_rcd = GF%Si
	write(*,*) "integrand_fun: ", trim(fpath), " was read!"

	if(orthbol_sw) call orthHS(GF%Hi, GF%Si, GF%Ddim)                    ! for orthogonal basis
	return
  end subroutine inputHSt



  !................................................................................
  ! generate hamiltonian/overlap matrixes h(t)/s(t) using capacitance network model
  !................................................................................
  subroutine inputHSt_cv(GF,ndim,tstep,orthbol_sw)
	type(GFtype) :: GF
	integer :: i, j, k, ndim, tstep, ierr, ndevatm
	real(R_KIND), allocatable, dimension(:) :: vm, qm
	logical :: orthbol_sw

	ndevatm = size(GF%ccinv,1)
	allocate(vm(ndevatm),qm(ndevatm))
	GF%Hi = GF%HD
	GF%Si = GF%SD
	GF%Si_rcd = GF%SD
	k=0
	do i=1,size(GF%ccinv,1)
		GF%rhotj(i) = GF%rhoti(i)
		GF%rhoti(i) = -GF%rhot0(i)
		do j=1,GF%numorb(i)
			k=k+1
			GF%rhoti(i)=GF%rhoti(i)+real(GF%rho(k,k))
		end do
	end do

	if(orthbol_sw) call orthHS(GF%Hi, GF%Si, GF%Ddim)                    ! for orthogonal basis
	do i=1,ndevatm
		qm(i)=-GF%rhoti(i)+GF%Lcup(i)*GF%VL(tstep)+GF%Rcup(i)*GF%VR(tstep)
	end do
	do i=1,ndevatm
		vm(i)=0.0
		do j=1,ndevatm
			vm(i)=vm(i)+GF%ccinv(i,j)*qm(j)
		end do
	end do
	! update Hi
	k=0
	do i=1,size(GF%ccinv,1)
		do j=1,GF%numorb(i)
			k=k+1
			GF%Hi(k,k)=GF%Hi(k,k)-vm(i)      ! use negative to coincide with the definition in the propagation U(t)
		end do
	end do
	deallocate(vm,qm)
	return
  end subroutine inputHSt_cv

  subroutine inputHSt_cv_next(GF,ndim,tstep,orthbol_sw)
	type(GFtype) :: GF
	integer :: i, j, k, ndim, tstep, ierr, ndevatm
	real(R_KIND), allocatable, dimension(:) :: vm, qm, rhonext
	logical :: orthbol_sw

	ndevatm = size(GF%ccinv,1)
	allocate(vm(ndevatm),qm(ndevatm), rhonext(ndevatm))
	GF%Hi = GF%HD
	GF%Si = GF%SD
	GF%Si_rcd = GF%SD
	k=0
	do i=1,size(GF%ccinv,1)
		rhonext(i) = 2.0*GF%rhoti(i) - GF%rhotj(i)
	end do

	if(orthbol_sw) call orthHS(GF%Hi, GF%Si, GF%Ddim)                    ! for orthogonal basis
	do i=1,ndevatm
		qm(i)=-rhonext(i)+GF%Lcup(i)*GF%VL(tstep)+GF%Rcup(i)*GF%VR(tstep)
	end do
	do i=1,ndevatm
		vm(i)=0.0
		do j=1,ndevatm
			vm(i)=vm(i)+GF%ccinv(i,j)*qm(j)
		end do
	end do
	! update Hi
	k=0
	do i=1,size(GF%ccinv,1)
		do j=1,GF%numorb(i)
			k=k+1
			GF%Hi(k,k)=GF%Hi(k,k)-vm(i)      ! use negative to coincide with the definition in the propagation U(t)
		end do
	end do
	deallocate(vm,qm,rhonext)
	return
  end subroutine inputHSt_cv_next



  !....................................................................
  ! retarded green function at t0
  !....................................................................
  real(R_KIND) function syst0(x)
	real(R_KIND) :: x
	complex(R_KIND) :: funs(Ddim,Ddim)

	funs = (x+c_i*delta*mag)*Sglb-Hglb
	call matrixinv(funs, Ddim)
	syst0= aimag(ctrace(Ddim,funs))*fermi(x,Ef)*(-2.0_R_KIND/pi)
  end function syst0


  !....................................................................
  ! retarded green function at t
  !....................................................................
  real(R_KIND) function syst(x)
	real(R_KIND) :: x
	complex(R_KIND) :: funs(Ddim,Ddim)

	funs = (x+c_i*delta*mag)*Sglb-Hglb
	call matrixinv(funs, Ddim)
	syst= aimag(ctrace(Ddim,funs))*fermi(x,Ef)*(-2.0_R_KIND/pi)
  end function syst


  !....................................................................
  ! real part and imaginary part of integrand1
  !....................................................................
  real(R_KIND) function integrand1_re(x)
	real(R_KIND) :: x
	complex(R_KIND) :: zval, res
	integer :: tid
	tid = OMP_GET_THREAD_NUM()+1

	zval = zbufs(1,tid)
	res = fermi(x,Ef)*1.0_R_KIND/(x-zval)
	integrand1_re = real(res,kind=R_KIND)
  end function integrand1_re
  real(R_KIND) function integrand1_im(x)
	real(R_KIND) :: x
	complex(R_KIND) :: zval, res
	integer :: tid
	tid = OMP_GET_THREAD_NUM()+1

	zval = zbufs(1,tid)
	res = fermi(x,Ef)*1.0_R_KIND/(x-zval)
	integrand1_im = aimag(res)
  end function integrand1_im


  !....................................................................
  ! real part and imaginary part of integrand2
  !....................................................................
  real(R_KIND) function integrand2_re(x)
	real(R_KIND) :: x
	complex(R_KIND) :: zval, res
	integer :: tid
	tid = OMP_GET_THREAD_NUM()+1

	zval = zbufs(1,tid)
	res = fermi(x,Ef)*exp(c_i*x*t_glb)/(x-zval)
	integrand2_re = real(res,kind=R_KIND)
  end function integrand2_re
  real(R_KIND) function integrand2_im(x)
	real(R_KIND) :: x
	complex(R_KIND) :: zval, res
	integer :: tid
	tid = OMP_GET_THREAD_NUM()+1

	zval = zbufs(1,tid)
	res = fermi(x,Ef)*exp(c_i*x*t_glb)/(x-zval)
	integrand2_im = aimag(res)
  end function integrand2_im

  !....................................................................
  ! exponential integral function to replace incomplete gamma function at alpha=0
  !....................................................................
  function cdig2(alpha, x) result(fn_val)
	complex (R_KIND), intent(in)  :: alpha
	complex (R_KIND), intent(in)  :: x
	complex (R_KIND)              :: fn_val
	if(abs(alpha)>1.0e-15) then
		write(*,*) "Warning! incorrect usage for incomplete gamma function (here, alpha == 0 only)."
	end if
	call e1z ( x, fn_val )
	return
  end function cdig2

  !....................................................................
  ! full-range incomplete gamma function X exp function
  !....................................................................
  function cdig_fun(alpha, x) result(fn_val)
	real (R_KIND) :: rval
	complex (R_KIND)  :: alpha
	complex (R_KIND)  :: x
	complex (R_KIND)  :: fn_val, fn_val1, fn_val2
	rval = abs(real(x,kind=R_KIND))
	if(rval<gamma_thrd) then
		!call e1z ( x, fn_val ) ! exponential integral (replacing incomplete gamma function at alpha=0)
		fn_val = cdig(alpha, x) ! incomplete gamma function
		fn_val = fn_val*exp(x)
	else if(rval>=gamma_thrd .and. rval<(gamma_thrd+5.0_R_KIND)) then
		fn_val1 = cdig(alpha, x) ! incomplete gamma function
		fn_val1 = fn_val1*exp(x)
		fn_val2 = cdig_inf_apprx(x)
		rval=(rval-gamma_thrd)/5.0_R_KIND
		fn_val = (1.0_R_KIND-rval)*fn_val1+rval*fn_val2
	else
		fn_val = cdig_inf_apprx(x)
	end if
	return
  end function cdig_fun
  function cdig_inf_apprx(xi) result(fn_vali)
	complex (R_KIND), intent(in)  :: xi
	integer :: nn
	real (R_KIND) :: cf(17)=(/ 1.0_R_KIND, 2.0_R_KIND, 6.0_R_KIND, 24.0_R_KIND, 120.0_R_KIND, 720.0_R_KIND, 5040.0_R_KIND, &
		& 40320.0_R_KIND, 362880.0_R_KIND, 3628800.0_R_KIND, 39916800.0_R_KIND, 479001600.0_R_KIND, 6227020800.0_R_KIND, &
		& 87178291200.0_R_KIND, 1307674368000.0_R_KIND, 20922789888000.0_R_KIND, 355687428096000.0_R_KIND /)
	complex (R_KIND)              :: z1
	complex (R_KIND)              :: fn_vali
	fn_vali = c_one/xi
	z1 = xi
	do nn = 1,16
		z1 = z1*xi
		if(mod(nn,2)==0) then
			fn_vali = fn_vali+cf(nn)/z1
		else
			fn_vali = fn_vali-cf(nn)/z1
		end if
	end do
	return
  end function cdig_inf_apprx

end module integrand_fun
