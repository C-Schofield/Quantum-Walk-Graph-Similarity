# Fortran Compiler
FC = ftn

# Fortran Compiler flags
FCFLAGS = -fopenmp -g -Wall #-O2#-Wextra
# -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan -fopenmp -Ofast

# Library flags
#LDFLAGS = -llapack -lblas 

# List of executables to be built
PROGRAMS = QW_Multiple QW_Changes
#PROG_OBJS = $(patsubst %.f90,%.o,$(SRCS))

# Find all source files, create a list of corresponding object files
#SRCS=$(wildcard *.f90)
#OBJS=$(patsubst %.f90,%.o,$(SRCS))

# Build all programs
all: $(PROGRAMS) 

# Dependencies
Walk_Functions.o: data_types.o
NodeAffinity_Functions.o: sort.o
QW_Similarity.o: data_types.o Walk_Functions.o Common_Functions.o
QW_Single.o: data_types.o Common_Functions.o Walk_Functions.o QW_Similarity.o NodeAffinity_Functions.o
QW_Single: data_types.o Common_Functions.o Walk_Functions.o QW_Similarity.o NodeAffinity_Functions.o
QW_Multiple.o: data_types.o Common_Functions.o Walk_Functions.o QW_Similarity.o NodeAffinity_Functions.o
QW_Multiple: data_types.o Common_Functions.o Walk_Functions.o QW_Similarity.o NodeAffinity_Functions.o
QW_Changes.o: data_types.o Common_Functions.o Walk_Functions.o QW_Similarity.o NodeAffinity_Functions.o Graph_Functions.o 
QW_Changes: data_types.o Common_Functions.o Walk_Functions.o QW_Similarity.o NodeAffinity_Functions.o Graph_Functions.o 
test.o: Common_Functions.o data_types.o Walk_Functions.o Graph_Functions.o 
test: Common_Functions.o data_types.o Walk_Functions.o Graph_Functions.o 

# Automatic compile	
%: %.o 
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

# Functions
.PHONY: clean veryclean cleanout debug all
debug:
# @echo "SRCS = $(SRCS)"
# @echo "OBJS = $(OBJS)"
# @echo "MODS = $(MODS)"
# @echo "MOD_OBJS = $(MOD_OBJS)"
	@echo "PROGRAMS = $(PROGRAMS)"
# @echo "PRG_OBJS = $(PRG_OBJS)"

clean:
	rm -f *.o *.mod

veryclean:
	rm -f *~ $(PROGRAMS) *.o *.mod 

cleanout:
	rm -f *.dat
