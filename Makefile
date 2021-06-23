all: print_qs.f90
	gfortran ./print_qs.f90 -lm -fdefault-real-8
