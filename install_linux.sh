workdir=$(cd "$(dirname "$0")"; pwd)

cd $workdir/src
ifort -c -ipo -O1 -xHost -qopenmp -parallel -autodouble -mp1 -mkl=parallel zgpadm.f90 dqagpe.f Global_Const.f90 cincgam.f90 special_functions.f90 d1mach.f dqelg.f dqk21.f dqpsrt.f MatrixOperation.f90 Gf.f90 Integrand.f90 Sys_Setup.f90 TDvs.f90



ifort  -O1 -ipo -xHost -qopenmp -parallel -autodouble -mp1 -mkl=parallel -o td_qmmm zgpadm.o dqagpe.o cincgam.o special_functions.o Global_Const.o d1mach.o dqelg.o dqk21.o dqpsrt.o MatrixOperation.o Integrand.o Gf.o  Sys_Setup.o TDvs.o 
	   
cd $workdir/
rm -f $workdir/src/*.mod
rm -f $workdir/src/*.o
