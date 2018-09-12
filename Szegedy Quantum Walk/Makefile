# Compiler
#FC = mpifort
# Magnus Compiler 
FC = ftn 

# Compiler flags
FCFLAGS = -fopenmp -g #-Wall -O2
#-Wall -Wextra -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan

# Library flags
#LDFLAGS = -llapack -lblas 

# List of executables to be built
PROGRAMS = Szegedy_Changes Szegedy_Multiple

# Build all programs
all: $(PROGRAMS)

# Dependencies
NodeAffinity_Functions.o: sort.o Common_Functions.o
Szegedy_Functions.o: Common_Functions.o
Szegedy_Similarity.o: Szegedy_Functions.o sort.o Common_Functions.o
Szegedy_Changes.o: Szegedy_Functions.o sort.o Common_Functions.o Graph_Functions.o NodeAffinity_Functions.o Szegedy_Similarity.o
Szegedy_Changes: Szegedy_Functions.o sort.o Common_Functions.o Graph_Functions.o NodeAffinity_Functions.o Szegedy_Similarity.o
Szegedy_Multiple.o: Szegedy_Functions.o sort.o Common_Functions.o Graph_Functions.o NodeAffinity_Functions.o Szegedy_Similarity.o
Szegedy_Multiple: Szegedy_Functions.o sort.o Common_Functions.o Graph_Functions.o NodeAffinity_Functions.o Szegedy_Similarity.o

# Automatic compile	
%: %.o 
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

# Clean functions
.PHONY: clean veryclean
clean:
	rm -f *.o *.mod

veryclean:
	rm -f *~ $(PROGRAMS) *.o *.mod 

cleanout:
	rm -f *~ $(PROGRAMS) *.o *.mod *.dat
