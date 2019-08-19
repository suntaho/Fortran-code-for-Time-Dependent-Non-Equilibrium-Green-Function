module matrixoperation
	use globaldef
	!use mkl95_lapack
	!use mkl95_precision
	implicit none

	contains
	!....................................................................
	! generate hamiltonian having orthogonal basis by given h/s  
	! h will be updated under definition of orthogonal basis, 
	! s will be changed into identity matrix
	!....................................................................
	subroutine orthHS(h, s, nd)
		implicit none

		integer :: nd
		complex(R_KIND) :: h(nd,nd), s(nd,nd)
		integer :: i, j
		complex(R_KIND), allocatable :: u(:,:), buf1(:,:), buf2(:,:), ucjg(:,:)
		
		allocate(buf1(nd,nd),buf2(nd,nd),ucjg(nd,nd),u(nd,nd))
		forall(i=1:nd,j=1:nd) u(i,j) = c_nil
		do j=1,nd                                                     ! upper half of s 
			do i=1,j
				u(i,j) = s(i,j)
			end do
		end do
		call cholesky(u, nd)
		ucjg = conjg(transpose(u))
		call matrixinv(u, nd)
		call matrixinv(ucjg, nd)
		! calculate h with orthogonal basis
		call matrixmul(buf1, h, u, nd)
		call matrixmul(buf2, ucjg, buf1, nd)
		! update h/s
		h = buf2
		forall(i=1:nd,j=1:nd) s(i,j) = c_nil
		forall(i=1:nd) s(i,i) = c_one		

		deallocate(buf1,buf2,u,ucjg)
		return
	end subroutine orthHS


	!....................................................................
	! Cal. cholesky decomposition of a positive definite matrix ia 
	!....................................................................
	subroutine cholesky(ia, na)
		implicit none
		integer :: na, info
		complex(R_KIND)  :: ia(na, na)
		!call zchdc(ia,na,na,work,jpvt,0,info)
		call zpotrf( 'U', na, ia, na, info )
		if (info .ne. 0) write(*,*) 'cholesky decomposition not successful: ', info
		return
	end subroutine cholesky
	
	
	!....................................................................
	! Solve complex-matrix inverse 
	!....................................................................
	subroutine matrixinv(ia, na)
		implicit none

		integer :: na, info, lda
		integer, allocatable :: ipiv(:)
		complex(R_KIND)  :: ia(na, na)
		complex(R_KIND), allocatable :: work(:)
		
		allocate(ipiv(na), work(4*na))
		call zgetrf(na, na, ia, na, ipiv, info) ! indexing
		call zgetri(na, ia, na, ipiv, work, 4*na, info) ! solving
		if (info .ne. 0) write(*,*) 'invrse-matrix not successful: ', info

		deallocate(ipiv, work)
		return
	end subroutine matrixinv
	
	
	!....................................................................
	! Multiply complex-matrix A by complex-matrix B: C=A*B
	!....................................................................
	subroutine matrixmul(ic, ia, ib, n)
		implicit none
		integer :: n
		complex(R_KIND) :: alpha=(1.0_R_KIND,0.0_R_KIND), beta=(0.0_R_KIND,0.0_R_KIND)
		complex(R_KIND)  :: ia(n, n), ib(n, n), ic(n, n)
		
		call zgemm('N','N', n, n, n, alpha, ia, n, ib, n, beta, ic, n)
		return
	end subroutine matrixmul
	
	    !....................................................................
	! Multiply complex-matrix A by complex-matrix B: C_no=A_nm*B_mo
	!....................................................................
	subroutine matrixmulnmo(ic, ia, ib, n, m, o)
		implicit none
		integer :: n, m, o
		complex(R_KIND) :: alpha=(1.0_R_KIND,0.0_R_KIND), beta=(0.0_R_KIND,0.0_R_KIND)
		complex(R_KIND) :: ia(n, m), ib(m, o), ic(n, o)
		
		call zgemm('N','N', n, o, m, alpha, ia, n, ib, m, beta, ic, n)
		return
	end subroutine matrixmulnmo

	!....................................................................
	! Solves A.x = l x for a general complex matrix A
	!....................................................................
	subroutine eigenvalues(A, e_val, e_vec_l, e_vec_r, n)
		complex(R_KIND) :: A(n,n), e_vec_l(n,n), e_vec_r(n,n), e_val(n)
		complex(R_KIND), allocatable :: work(:)
		real(R_KIND), allocatable :: rwork(:)
		integer :: n, error

		ALLOCATE(work(4*n), rwork(4*n), STAT=error)
		if(error /= 0)stop "eigenvalues: allocation problem."

		if(R_KIND==4)then
			call cgeev('N', 'V', n, A, n, e_val, e_vec_l, n, e_vec_r, n, work, 4*n, rwork, error)
			if(error /= 0)stop "eigenvalues: IMKL error."
		else
			call zgeev('N', 'V', n, A, n, e_val, e_vec_l, n, e_vec_r, n, work, 4*n, rwork, error)
			if(error /= 0)stop "eigenvalues: IMKL error."
		end if
		DEALLOCATE(work, rwork)
	end subroutine eigenvalues
  
	
	!....................................................................
	! Solves A.x = l x for a complex Hermitian matrix A, where
	! e_vec is initialized as the upper triangle matrix of A
	!....................................................................
	subroutine heigenvalues_upper(e_vec, e_val, n)
		complex(R_KIND) :: e_vec(n,n)
		real(R_KIND) :: e_val(n)
		complex(R_KIND), allocatable :: work(:)
		real(R_KIND), allocatable :: rwork(:)
		integer :: n, error

		ALLOCATE(work(4*n), rwork(3*n-2), STAT=error)
		if(error /= 0)stop "eigenvalues: allocation problem."

		if(R_KIND==4)then
			call cheev('V', 'U', n, e_vec, n, e_val, work, 4*n, rwork, error)
			if(error /= 0)stop "eigenvalues: IMKL error."
		else
			call zheev('V', 'U', n, e_vec, n, e_val, work, 4*n, rwork, error)
			if(error /= 0)stop "eigenvalues: IMKL error."
		end if
		DEALLOCATE(work, rwork)
	end subroutine heigenvalues_upper
	
	

	!....................................................................
	! Solves A.x = l B.x for general complex matrices A and B. 
	! Note! A and B are overwritten on exit
	!....................................................................
	subroutine generalized_eigenvalues(A, B, e_val, e_vec_l, e_vec_r, n)
		complex(R_KIND) :: A(n,n), B(n,n), e_vec_l(n,n), e_vec_r(n,n), e_val(n)
		complex(R_KIND), allocatable ::  alpha(:), beta(:), work(:)
		real(R_KIND), allocatable :: rwork(:)
		integer :: i, n, error

		ALLOCATE(alpha(n), beta(n), work(8*n), rwork(8*n), STAT=error)
		if(error /= 0)stop "generalized_eigenvalues: allocation problem."

		if(R_KIND==4)then
			call cggev('N', 'V', n, A, n, B, n, alpha, beta, e_vec_l, n, e_vec_r, n, work, 8*n, rwork, error)
			if(error /= 0)stop "generalized_eigenvalues: LaPack error."
		else
			call zggev('N', 'V', n, A, n, B, n, alpha, beta, e_vec_l, n, e_vec_r, n, work, 8*n, rwork, error)
			if(error /= 0)stop "generalized_eigenvalues: LaPack error."
		end if
		forall(i=1:n) e_val(i) = alpha(i)/beta(i)

		DEALLOCATE(alpha, beta, work, rwork)
	end subroutine generalized_eigenvalues
	!....................................................................


	!----------------------------------------------------------------------
	! cal. transpose_conjgate(vec_f)*matrix*vec_i (complex)
	!----------------------------------------------------------------------
	subroutine vec_mat_vec(v1, mi, v2, na, prdt)
		integer :: na
		complex(R_KIND) :: prdt
		complex(R_KIND) :: v1(na), v2(na), mi(na,na)
		complex(R_KIND), allocatable :: vi(:)

		allocate(vi(na))
		vi = c_nil
		call zgemv('N', na, na, c_one, mi, na, v2, 1, c_nil, vi, 1)
		prdt = sum(conjg(v1)*vi)

		deallocate(vi)
	end subroutine vec_mat_vec


	!----------------------------------------------------------------------
	! calulate matrix exponential e = exp(a)
	!----------------------------------------------------------------------
	subroutine expmk ( n, a, e )
		real(R_KIND) :: t
		integer :: ideg, n, lwsp, iexph, ns, iflag, ipiv(n), i, j
		complex(R_KIND) :: a(n,n), e(n,n)
		complex(R_KIND), allocatable :: wsp(:)
		
		t = 1.0_R_KIND
		ideg = 7
		lwsp = 4*n*n+ideg+2
		allocate(wsp(lwsp))
		call ZGPADM(ideg,n,t,a,n,wsp,lwsp,ipiv,iexph,ns,iflag)
		if(iflag <0) write(*,*) "expmk: warning on computing matrix exponential function!"
		do i = 1,n
		  do j = 1,n
			e(j,i) = wsp(iexph)
			iexph = iexph+1
		  end do
		end do
		deallocate(wsp)
	end subroutine expmk

	
	!----------------------------------------------------------------------
	! calulate matrix exponential e = exp(a)
	!----------------------------------------------------------------------
	subroutine c8mat_expm1 ( n, a, e )
		integer :: n, i, j
		complex(R_KIND) :: a(n,n), e(n,n)
		complex(R_KIND), allocatable :: vl(:,:), vr(:,:), ai(:,:), val(:)

		allocate(vl(n,n),vr(n,n),ai(n,n),val(n))
		ai=a
		call eigenvalues(ai, val, vl, vr, n)
		vl=vr
		call matrixinv(vl, n)
		forall(i=1:n,j=1:n) ai(i,j)=c_nil
		forall(i=1:n) ai(i,i)=exp(val(i))
		call matrixmul(e, ai, vl, n)
		call matrixmul(ai, vr, e, n)
		e=ai

		deallocate(vl,vr,ai,val)
		return
	end subroutine c8mat_expm1
	
	

	!....................................................................
	! Checks if the complex matrix A is hermitian.
	!....................................................................
	logical function is_hermitian(r, A)
		integer, intent(in) :: r
		complex(R_KIND), dimension(r,r) :: A
		integer :: i, j
		
		! Default value
		is_hermitian = .TRUE.

		! Check square
		if(r /= size(A,2))then
		  is_hermitian = .FALSE.
		  return
		end if

		! Check Aij = conjg(Aji)
		do i = 1, r
		  do j = 1, r
			if(abs(A(i,j) - conjg(A(j,i))) > delta)then
			  is_hermitian = .FALSE.
			  return
			end if
		  end do
		end do
	end function is_hermitian

	!....................................................................
	! Computes the trace of a complex matrix A
	!....................................................................
	complex(R_KIND) function ctrace(r, A)
		integer, intent(in) :: r
		complex(R_KIND), dimension(r,r) :: A
		integer :: i

		ctrace = (0.0_R_KIND,0.0_R_KIND)
		do i = 1, r
		  ctrace = ctrace + A(i,i)
		end do
	end function ctrace

	
	!....................................................................
	! Checks if the imaginery part of every element of the trace is 
	! negative (positive DoS for each orbital).
	! Input: G - the GF matrix;
	! Output: flag - index for positive element
	!....................................................................
	integer function CheckDoS(r, G)
		integer, intent(in) :: r
		complex(R_KIND), dimension(r,r) :: G
		integer :: i
		
		CheckDoS = 0
		do i = 1, r
			if(aimag(G(i,i)) > nil)then
				CheckDoS = i
				exit
			end if
		end do
	end function CheckDoS

end module matrixoperation

