!======================================================================
!  This program calculates the time-dependent densities and currents 
!  by multiscale (QM/MM; not in windows version) molecular system, emplying wide-band limit approximation.
!
!  Version: 0.11 / 08.06.2018
!======================================================================
    
program TDvs
    use globaldef
    use sys_setup
    use gfun
	use integrand_fun
    implicit none
	
	type(GFtype) :: tGF 

    ! start timing
    call system_clock(time(0),rate)
   
    ! read initial parameters 
	call read_t0(tGF)
    ! calculate transmission function for steady state
    if(genTf) call transmission(tGF,tGF%Ef-erng,tGF%Ef+erng,deng)
   
    ! calculate self-energy of electrode at fermi energy for Lambda and Gamma matrixes
	call TDNEGF(tGF)
   
   
    ! free memory
	call free_array(tGF)
	
    ! end timing and calculate cpu-time
    call system_clock(time(1))
    write(*,'(A45,F12.3)') "runnin-time (sec.) by: ", (time(1)-time(0))/(1.0_R_KIND*rate)
	
end program TDvs

