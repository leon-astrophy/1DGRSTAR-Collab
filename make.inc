# look for gfortran compiler #
F90=gfortran

# fortran flags #
F90FLAGS= -ffree-line-length-512 -Ofast -O4 
LDFLAGS= -ffree-line-length-512 -Ofast -O4 

# other flags #
MODINC="-I ./"