MODE=PRODUCTION

FC_PC=mpif90
GCC=g++

SRC=./src
OBJ=./obj

# tpl
CGAL=./lib/CGAL-5.5.2
BOOST=./lib/boost/1.72.0

# --- Optimization Flags
ifeq ($(MODE), DEBUG)
OPT_FLAGS = -O0 -g
endif

ifeq ($(MODE), PRODUCTION)
OPT_FLAGS = -O5 -pg
endif

# --- Linker Flags
LD_PC_FLAGS=$(OPT_FLAGS)

# --- Compiler Flags
FC_PC_FLAGS=$(OPT_FLAGS) -c -cpp -ffree-form -ffree-line-length-none -fopenmp -DPETSC_AVOID_MPIF_H -fbounds-check

# ---- COMPILE -----------------------------------------------------------------
SRCS=$(wildcard ./src/*.f90)  					# All f90 source files
OBJS=$(patsubst ./src/%.f90,./obj/%.o,$(SRCS))  # All all object files

# --- OUTPUT FILENAME
EXEC=bin/INTERSECTION

# --- Default rule to build all
.PHONY: all
all: MODULES $(OBJS) $(EXEC) # First build extra modules, then object files, then the executable

# --- CGAL interface
.PHONY: CGAL
CGAL:
	g++ -I$(CGAL)/include -I$(BOOST)/include  -c $(SRC)/nefInterface.cpp -o $(OBJ)/nefInterface.o

.PHONY: MODULES
MODULES:
	$(FC_PC) $(FC_PC_FLAGS) $(SRC)/MODULES.f90  -o $(OBJ)/MODULES.o -J$(OBJ)
	$(FC_PC) $(FC_PC_FLAGS) $(SRC)/MOD_POLY.f90 -o $(OBJ)/MOD_POLY.o -J$(OBJ)
	$(FC_PC) $(FC_PC_FLAGS) $(SRC)/MOD_OPENFOAM_GRID.f90 -o $(OBJ)/MOD_OPENFOAM_GRID.o -J$(OBJ)
	$(FC_PC) $(FC_PC_FLAGS) $(SRC)/MOD_INTERSECTION.f90 -o $(OBJ)/MOD_INTERSECTION.o -J$(OBJ)

# --- Compile rules for individual objects
$(OBJS):./obj/%.o: ./src/%.f90 ./obj/MODULES.o
	$(FC_PC) $(FC_PC_FLAGS) $(OBJ)/nefInterface.o  $< -o $@ -J$(OBJ)

# --- Linking all objects
$(EXEC): $(OBJS)
	$(FC_PC) $(OBJ)/nefInterface.o -o $@  $(OBJS) -lstdc++ -lmpfr -lgmp -lgmpxx -L$(BOOST)/lib -lboost_filesystem
 
# --- Delete everything except sources
.PHONY: clean
clean:
	-find ./obj/ -type f -name "*.o" ! -wholename "./obj/nefInterface.o" -exec rm {} \;
	-rm -f obj/*.mod bin/INTERSECTION

cleanCGAL:
	-rm obj/nefInterface.o
