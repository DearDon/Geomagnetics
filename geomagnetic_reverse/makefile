main: sub1.o sub.o main.o 
	gfortran  -o main sub1.o sub.o main.o
	rm *.o *.mod
main.o: main.f
	gfortran -g -c main.f
sub1.o: sub1.f
	gfortran -g -c sub1.f
sub.o: sub.f
	gfortran -g -c sub.f

