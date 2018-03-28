main: sub.o main.o 
	gfortran -o main main.o sub.o
	rm *.mod *.o
sub.o: sub.f90
	gfortran -c -g sub.f90
main.o: main.f90
	gfortran -c -g main.f90
