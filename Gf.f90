!======================================================================
!    calculations associated with Green's functions 
!======================================================================
	
module gfun
  use globaldef
  use omp_lib
  use matrixoperation
  implicit none  
  
  contains
  
  
  !....................................................................
  ! calculate transmission functino for steady state
  ! ebgn: beginning energy for spectrum
  ! eend: ending energy for spectrum
  ! de:   energy interval for spectrum  
  !....................................................................
  subroutine transmission(GF,ebgn,eend,de)
	type(GFtype) :: GF
	integer :: i, j, ierr, nee
	real(R_KIND) :: ebgn, eend, de
	real(R_KIND) :: Tf, ee
	real(R_KIND), allocatable, dimension(:,:) :: Tfdata(:,:)
	complex(R_KIND), allocatable, dimension(:,:) :: idty, GR, GA
	complex(R_KIND), allocatable, dimension(:,:) :: buf1, buf2, buf3
	complex(R_KIND), allocatable, dimension(:,:) :: GammaL, GammaR
	character(255) :: fpath

	
	nee = nint((eend-ebgn)/de) + 1
	allocate(Tfdata(nee,2))
	allocate(idty(GF%Ddim,GF%Ddim),GR(GF%Ddim,GF%Ddim),GA(GF%Ddim,GF%Ddim))
	allocate(buf1(GF%Ddim,GF%Ddim),buf2(GF%Ddim,GF%Ddim),buf3(GF%Ddim,GF%Ddim))
	allocate(GammaL(GF%Ddim,GF%Ddim),GammaR(GF%Ddim,GF%Ddim))	
	forall(i=1:GF%Ddim,j=1:GF%Ddim) idty(i,j) = c_nil
	forall(i=1:GF%Ddim) idty(i,i) = c_one
	
	
	do j = 1, nee
		ee = ebgn + (j-1.0_R_KIND)*de
		! calculate self-energy of L/R electrodes and broadening energy
		call selfenergyLR(GF,ee)
		GammaL = -2.0_R_KIND*aimag(GF%SEL)
		GammaR = -2.0_R_KIND*aimag(GF%SER)

		! calculate advanced/retarded green function of Ddim
		GR = (ee+c_i*delta)*GF%SD - GF%HD - GF%SEL - GF%SER
		call matrixinv(GR,GF%Ddim)
		GA = dconjg(transpose(GR))
		! check if the imaginary part of surface GF is negative (positive DoS for each orbital)
		i = CheckDoS(GF%Ddim, GR)
		if(i>0) stop "gfun: retarded green function of Device is error!"
		
		! calculate transmission function
		call matrixmul(buf1, GammaL, GA, GF%Ddim)
		call matrixmul(buf2, GR, buf1, GF%Ddim)
		call matrixmul(buf3, GammaR, buf2, GF%Ddim)
		Tf = real(ctrace(GF%Ddim, buf3),kind=R_KIND)
		
		! show message
		write(*,'(A15, F7.4, F7.4, A3, F7.3, A3, A,\)') "transmission: ", ee, Tf, " , ", 100.0_R_KIND*(ee-ebgn)/(eend-ebgn), "% ", creturn
		Tfdata(j,1) = ee
		Tfdata(j,2) = Tf
	end do
	write(*,*) ""
	
	! export transmission function
	fpath = trim(fout) // "transmission.txt"
	open(UNIT=112, FILE=trim(fpath), STATUS="REPLACE", IOSTAT=ierr)
	if(ierr /= 0) stop "gfun: openning transmission.txt is error!"	
	do j = 1, nee
		write(112,*) Tfdata(j,1), tab, Tfdata(j,2)
	end do
	close(112)
	! export steady current
	i = nint((nee+1.0_R_KIND)/2.0_R_KIND)
	fpath = trim(fout) // "steady_current.txt"
	open(UNIT=113, FILE=trim(fpath), STATUS="REPLACE", IOSTAT=ierr)
	if(ierr /= 0) stop "gfun: openning steady_current.txt in default path is error!"	
	do j = 0, nint((nee-1.0_R_KIND)/2.0_R_KIND)
		write(113,*) 2.0_R_KIND*j*de, tab, de*sum(Tfdata(i-j:i+j,2))*e_over_hpa_pi
	end do
	close(113)

	! free arrays
	deallocate(idty,GR,GA,buf1,buf2,buf3,GammaL,GammaR,Tfdata)
  end subroutine transmission
  
  
  
  !....................................................................
  ! calculate selfenergy of L/R electrodes, with frequency/energy representation
  !....................................................................
  subroutine selfenergyLR(GF,Ex)
	type(GFtype) :: GF
	integer :: i, j
	real(R_KIND) :: Ex
	complex(R_KIND), allocatable, dimension(:,:) :: heff
	complex(R_KIND), allocatable, dimension(:,:) :: buf1, buf2, bufn, bufm

	allocate(heff(GF%tdim,GF%tdim))
	heff = GF%Ht - Ex*GF%St                       ! calculate heff defined by orthonormal basis
	forall(i=1:GF%tdim) heff(i,i) = heff(i,i) + Ex*c_one

	! for left electrode
	allocate(buf1(GF%Ldim/2,GF%Ddim),buf2(GF%Ddim,GF%Ddim))
	allocate(bufn(GF%Ddim,GF%Ldim/2),bufm(GF%Ldim/2,GF%Ddim))
	!bufn = heff(1:GF%Ddim,GF%Ddim+1:GF%Ldim/2+GF%Ddim)
	!bufm = heff(GF%Ddim+1:GF%Ldim/2+GF%Ddim,1:GF%Ddim)
	bufn = heff(GF%Ldim/2+1:GF%tdim-GF%Rdim/2,1:GF%Ldim/2)
	bufm = heff(1:GF%Ldim/2,GF%Ldim/2+1:GF%tdim-GF%Rdim/2)
	call surfaceGF(GF ,Ex)
	!call surfaceGF_model(GF,Ex)                  ! other method for surface green function
	call matrixmulnmo(buf1, GF%GLs, bufm, GF%Ldim/2, GF%Ldim/2, GF%Ddim)
	call matrixmulnmo(buf2, bufn, buf1, GF%Ddim, GF%Ldim/2, GF%Ddim)
	GF%SEL = buf2
	deallocate(buf1,buf2,bufn,bufm)

	! for right electrode
	allocate(buf1(GF%Rdim/2,GF%Ddim),buf2(GF%Ddim,GF%Ddim))
	allocate(bufn(GF%Ddim,GF%Rdim/2),bufm(GF%Rdim/2,GF%Ddim))
	!bufn = heff(1:GF%Ddim,GF%Ldim+GF%Ddim+1:GF%Ldim+GF%Ddim+GF%Rdim/2)
	!bufm = heff(GF%Ldim+GF%Ddim+1:GF%Ldim+GF%Ddim+GF%Rdim/2,1:GF%Ddim)
	bufn = heff(GF%Ldim/2+1:GF%tdim-GF%Rdim/2,GF%tdim-GF%Rdim/2+1:GF%tdim)
	bufm = heff(GF%tdim-GF%Rdim/2+1:GF%tdim,GF%Ldim/2+1:GF%tdim-GF%Rdim/2)
	! call surfaceGF(GF,Ex)                       ! done as calculating left electrode
	call matrixmulnmo(buf1, GF%GRs, bufm, GF%Rdim/2, GF%Rdim/2, GF%Ddim)
	call matrixmulnmo(buf2, bufn, buf1, GF%Ddim, GF%Rdim/2, GF%Ddim)
	GF%SER = buf2
	deallocate(buf1,buf2,bufn,bufm)
	
	! free other arrays
	deallocate(heff)
  end subroutine selfenergyLR
  
  
  !....................................................................
  ! calculate surface green function of electrodes by other model
  !....................................................................
  subroutine surfaceGF_model(GF,Ex)
	type(GFtype) :: GF
	integer :: i, j, ndim, maxitr
	real(R_KIND) :: Ex, var
	complex(R_KIND), allocatable, dimension(:,:) :: Tp, Tm, Sp, Sm, Q, pr_u, pr_v
	complex(R_KIND), allocatable, dimension(:,:) ::	idty, buf1, buf2, buf3
	complex(R_KIND), allocatable, dimension(:,:,:) :: u, v
	
	maxitr = 100
	
	! for left electrode
	ndim = GF%Ldim/2
	allocate(idty(ndim,ndim))
	allocate(Tp(ndim,ndim),Tm(ndim,ndim),Sp(ndim,ndim),Sm(ndim,ndim))
	allocate(Q(ndim,ndim),pr_u(ndim,ndim),pr_v(ndim,ndim))
	allocate(u(ndim,ndim,2),v(ndim,ndim,2))
	allocate(buf1(ndim,ndim),buf2(ndim,ndim),buf3(ndim,ndim))
	forall(i=1:ndim,j=1:ndim) idty(i,j) = c_nil   ! define unitary matrix
	forall(i=1:ndim) idty(i,i) = c_one
	! prepare initial parameter
	Q = (Ex+c_i*delta)*GF%SL(1:ndim,1:ndim)-GF%HL(1:ndim,1:ndim)
	call matrixinv(Q, ndim)
	buf1 = GF%HL(ndim+1:2*ndim,1:ndim)-Ex*GF%SL(ndim+1:2*ndim,1:ndim)
	call matrixmul(buf2, Q, buf1, ndim)
	u(:,:,1) = buf2
	buf1 = GF%HL(1:ndim,ndim+1:2*ndim)-Ex*GF%SL(1:ndim,ndim+1:2*ndim)
	call matrixmul(buf2, Q, buf1, ndim)
	v(:,:,1) = buf2
	Tp = u(:,:,1)
	Tm = v(:,:,1)
	pr_u = u(:,:,1)
	pr_v = v(:,:,1)
	do i = 1, maxitr
		! i-th step
		call matrixmul(buf1,u(:,:,1), v(:,:,1),ndim)
		call matrixmul(buf2,v(:,:,1),u(:,:,1),ndim)
		Q = idty-buf1-buf2
		call matrixinv(Q, ndim)                                       ! (1 - uo.vo - vo.uo)^-1
		call matrixmul(buf1,u(:,:,1),u(:,:,1),ndim)
		call matrixmul(buf2,Q,buf1,ndim)
		u(:,:,2) = buf2                                               ! u(i) = Inv(1-u(i-1)*v(i-1)-v(i-1)*u(i-1))*u(i-1)^2
		call matrixmul(buf3,pr_v,buf2,ndim)
		Tp = Tp+buf3                                                  ! T = u0 + v0*u1 + v0*v1*u2 + ...
		call matrixmul(buf1,v(:,:,1),v(:,:,1),ndim)
		call matrixmul(buf2,Q,buf1,ndim)
		v(:,:,2) = buf2                                               ! v(i) = Inv(1-u(i-1)*v(i-1)-v(i-1)*u(i-1))*v(i-1)^2
		call matrixmul(buf3,pr_u,buf2,ndim)
		Tm = Tm+buf3                                                  ! Tbar = v0 + u0*v1 + u0*u1*v2 + ...
		buf1 = u(:,:,1)
		buf2 = u(:,:,2)
		var = abs(ctrace(ndim,buf1)-ctrace(ndim, buf2))
		if(var<delta .and. i>3) exit
		! prepare i+1 step
		u(:,:,1) = u(:,:,2)
		v(:,:,1) = v(:,:,2)
		call matrixmul(buf1,pr_u,u(:,:,1),ndim)
		pr_u = buf1
		call matrixmul(buf1,pr_v,v(:,:,1),ndim)
		pr_v = buf1
	end do
	if (i>= maxitr) write (*,*) "gfun: Warning on convergence for surface green function of L electrode. err= ", var
	! G = Inv(E - H(0) - H(-1)*Tbar)
	buf1 = GF%HL(ndim+1:2*ndim,1:ndim)-Ex*GF%SL(ndim+1:2*ndim,1:ndim)
	call matrixmul(buf2,buf1,Tm,ndim)
	Q = (Ex+c_i*delta)*GF%SL(1:ndim,1:ndim)-GF%HL(1:ndim,1:ndim)-buf2
	call matrixinv(Q, ndim)
	GF%GLs = Q
	! check if the imaginary part of surface GF is negative (positive DoS for each orbital)
	i = CheckDoS(ndim, Q)
	if(i>0) stop "gfun: surface green function for L electrode is error!"
	deallocate(Tp,Tm,Sp,Sm,u,v,pr_u,pr_v,Q,buf1,buf2,buf3,idty)
	
	
	
	! for right electrode
	ndim = GF%Rdim/2
	allocate(idty(ndim,ndim))
	allocate(Tp(ndim,ndim),Tm(ndim,ndim),Sp(ndim,ndim),Sm(ndim,ndim))
	allocate(Q(ndim,ndim),pr_u(ndim,ndim),pr_v(ndim,ndim))
	allocate(u(ndim,ndim,ndim),v(ndim,ndim,ndim))
	allocate(buf1(ndim,ndim),buf2(ndim,ndim),buf3(ndim,ndim))
	forall(i=1:ndim,j=1:ndim) idty(i,j) = c_nil   ! define unitary matrix
	forall(i=1:ndim) idty(i,i) = c_one
	! prepare initial parameter
	Q = (Ex+c_i*delta)*GF%SR(1:ndim,1:ndim)-GF%HR(1:ndim,1:ndim)
	call matrixinv(Q, ndim)
	buf1 = GF%HR(ndim+1:2*ndim,1:ndim)-Ex*GF%SR(ndim+1:2*ndim,1:ndim)
	call matrixmul(buf2, Q, buf1, ndim)
	u(:,:,1) = buf2
	buf1 = GF%HR(1:ndim,ndim+1:2*ndim)-Ex*GF%SR(1:ndim,ndim+1:2*ndim)
	call matrixmul(buf2, Q, buf1, ndim)
	v(:,:,1) = buf2
	Tp = u(:,:,1)
	Tm = v(:,:,1)
	pr_u = u(:,:,1)
	pr_v = v(:,:,1)
	do i = 1, maxitr
		! i-th step
		call matrixmul(buf1,u(:,:,1), v(:,:,1),ndim)
		call matrixmul(buf2,v(:,:,1),u(:,:,1),ndim)
		Q = idty-buf1-buf2
		call matrixinv(Q, ndim)                                       ! (1 - uo.vo - vo.uo)^-1
		call matrixmul(buf1,u(:,:,1),u(:,:,1),ndim)
		call matrixmul(buf2,Q,buf1,ndim)
		u(:,:,2) = buf2                                               ! u(i) = Inv(1-u(i-1)*v(i-1)-v(i-1)*u(i-1))*u(i-1)^2
		call matrixmul(buf3,pr_v,buf2,ndim)
		Tp = Tp+buf3                                                  ! T = u0 + v0*u1 + v0*v1*u2 + ...
		call matrixmul(buf1,v(:,:,1),v(:,:,1),ndim)
		call matrixmul(buf2,Q,buf1,ndim)
		v(:,:,2) = buf2                                               ! v(i) = Inv(1-u(i-1)*v(i-1)-v(i-1)*u(i-1))*v(i-1)^2
		call matrixmul(buf3,pr_u,buf2,ndim)
		Tm = Tm+buf3                                                  ! Tbar = v0 + u0*v1 + u0*u1*v2 + ...
		buf1 = u(:,:,1)
		buf2 = u(:,:,2)
		var = abs(ctrace(ndim,buf1)-ctrace(ndim, buf2))
		if(var<delta .and. i>3) exit
		! prepare i+1 step
		u(:,:,1) = u(:,:,2)
		v(:,:,1) = v(:,:,2)
		call matrixmul(buf1,pr_u,u(:,:,1),ndim)
		pr_u = buf1
		call matrixmul(buf1,pr_v,v(:,:,1),ndim)
		pr_v = buf1
	end do
	if (i>= maxitr) write (*,*) "gfun: Warning on convergence for surface green function of R electrode. err= ",var
	! G = Inv(E - H(0) - H(1)*T)
	buf1 = GF%HR(1:ndim,ndim+1:2*ndim)-Ex*GF%SR(1:ndim,ndim+1:2*ndim)
	call matrixmul(buf2,buf1,Tp,ndim)
	Q = (Ex+c_i*delta)*GF%SR(1:ndim,1:ndim)-GF%HR(1:ndim,1:ndim)-buf2
	call matrixinv(Q, ndim)
	GF%GRs = Q
	! check if the imaginary part of surface GF is negative (positive DoS for each orbital)
	i = CheckDoS(ndim, Q)
	if(i>0) stop "gfun: surface green function for R electrode is error!"
	deallocate(Tp,Tm,Sp,Sm,u,v,pr_u,pr_v,Q,buf1,buf2,buf3,idty)
  end subroutine surfaceGF_model
  

  !....................................................................
  ! calculate surface green function of electrodes, method described in note
  !....................................................................
  subroutine surfaceGF(GF,Ex)
	type(GFtype) :: GF
	integer :: i, j, ndim, maxitr
	real(R_KIND) :: Ex, var
	complex(R_KIND), pointer, dimension(:,:) :: zeta, eta, a, b
	complex(R_KIND), pointer, dimension(:,:) ::	idty, heff, buf1, buf2, bufinv
	maxitr = 100

	! for left electrode
	ndim = GF%Ldim/2
	allocate(idty(ndim,ndim))
	allocate(zeta(ndim,ndim),eta(ndim,ndim),a(ndim,ndim),b(ndim,ndim))
	allocate(heff(GF%Ldim,GF%Ldim),buf1(ndim,ndim),buf2(ndim,ndim),bufinv(ndim,ndim))
	heff = GF%HL-Ex*GF%SL                         ! calculate heff defined by orthonormal basis
	forall(i=1:GF%Ldim) heff(i,i) = heff(i,i) + Ex*c_one
	forall(i=1:ndim,j=1:ndim) idty(i,j) = c_nil   ! define unitary matrix
	forall(i=1:ndim) idty(i,i) = c_one
	zeta = heff(1:ndim,1:ndim)
	eta = zeta
	a = heff(1:ndim,ndim+1:GF%Ldim)
	b = heff(ndim+1:GF%Ldim,1:ndim)
	i = 0
	var = 999.0_R_KIND
	do while (i<=maxitr .and. var>0.1_R_KIND*delta)
		bufinv = (Ex+c_i*delta)*idty - eta
		call matrixinv(bufinv,ndim)
		! update zeta
		call matrixmul(buf1, bufinv, b, ndim)
		call matrixmul(buf2, a, buf1, ndim)
		zeta = zeta + buf2
		var = maxval(abs(buf2))                   ! 
		! update eta
		eta = eta + buf2
		call matrixmul(buf1, bufinv, a, ndim)
		call matrixmul(buf2, b, buf1, ndim)
		eta = eta + buf2
		! update a
		call matrixmul(buf1, bufinv, a, ndim)
		call matrixmul(buf2, a, buf1, ndim)
		a = buf2
		! update b
		call matrixmul(buf1, bufinv, b, ndim)
		call matrixmul(buf2, b, buf1, ndim)
		b = buf2
		! loop
		i = i + 1
	end do
	if (var>delta) write (*,*) "gfun: Warning on convergence for surface green function of L electrode, err = ", var
	buf1 = (Ex+c_i*delta)*idty - zeta
	call matrixinv(buf1, ndim)
	GF%GLs = buf1
	! check if the imaginary part of surface GF is negative (positive DoS for each orbital)
	i = CheckDoS(ndim, buf1)
	if(i>0) stop "gfun: surface green function for L electrode is error!"
	deallocate(zeta,eta,a,b,heff,buf1,buf2,bufinv,idty)
	
	! for right electrode
	ndim = GF%Rdim/2
	allocate(idty(ndim,ndim))
	allocate(zeta(ndim,ndim),eta(ndim,ndim),a(ndim,ndim),b(ndim,ndim))
	allocate(heff(GF%Rdim,GF%Rdim),buf1(ndim,ndim),buf2(ndim,ndim),bufinv(ndim,ndim))
	heff = GF%HR-Ex*GF%SR                         ! calculate heff defined by orthonormal basis
	forall(i=1:GF%Rdim) heff(i,i) = heff(i,i)+Ex*c_one
	forall(i=1:ndim,j=1:ndim) idty(i,j) = c_nil   ! define unitary matrix
	forall(i=1:ndim) idty(i,i) = c_one
	zeta = heff(1:ndim,1:ndim)
	eta = zeta
	a = heff(1:ndim,ndim+1:GF%Rdim)
	b = heff(ndim+1:GF%Rdim,1:ndim)
	i = 0
	var = 999.0_R_KIND
	do while (i<=maxitr .and. var>0.1_R_KIND*delta)
		bufinv = (Ex+c_i*delta)*idty - eta
		call matrixinv(bufinv,ndim)
		! update zeta
		call matrixmul(buf1, bufinv, b, ndim)
		call matrixmul(buf2, a, buf1, ndim)
		zeta = zeta + buf2
		var = maxval(abs(buf2))                   ! 
		! update eta
		eta = eta + buf2
		call matrixmul(buf1, bufinv, a, ndim)
		call matrixmul(buf2, b, buf1, ndim)
		eta = eta + buf2
		! update a
		call matrixmul(buf1, bufinv, a, ndim)
		call matrixmul(buf2, a, buf1, ndim)
		a = buf2
		! update b
		call matrixmul(buf1, bufinv, b, ndim)
		call matrixmul(buf2, b, buf1, ndim)
		b = buf2
		! loop
		i = i + 1
	end do
	if (var>delta) write (*,*) "gfun: Warning on convergence for surface green function of R electrode, err = ", var
	buf1 = (Ex+c_i*delta)*idty - zeta
	call matrixinv(buf1, ndim)
	GF%GRs = buf1
	! check if the imaginary part of surface GF is negative (positive DoS for each orbital)
	i = CheckDoS(ndim, buf1)
	if(i>0) stop "gfun: surface green function for R electrode is error!"
	deallocate(zeta,eta,a,b,heff,buf1,buf2,bufinv,idty)
  end subroutine surfaceGF
  
	
end module gfun