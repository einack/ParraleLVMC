
SRC_DIR=src
OBJ_DIR=objects
BIN_DIR=bin
INCLUDE_DIR=objects
FFLAGS= -std=f2003 -finline-functions -fcheck=all -I$(INCLUDE_DIR) #-Wall

#a_mpi.out: random.f90 module_ising.f90 ising.f90 
#	mpif90 -DDEBUG -o a_mpi.out random.f90 module_ising.f90 ising.f90

#a.out: random.f90 module_ising.f90 ising.f90 
#	gfortran -o a.out random.f90 module_ising.f90 ising.f90
ising_mpi.x: setup $(BIN_DIR)/ising_mpi.x 

$(BIN_DIR)/ising_mpi.x: $(OBJ_DIR)/ising.o $(OBJ_DIR)/functions.o $(OBJ_DIR)/random.o 
	mpif90  -o $@ $^
	
random: setup $(OBJ_DIR)/random.o 

$(OBJ_DIR)/random.o: $(SRC_DIR)/random.f90  
	mpif90 $(FFLAGS)  -c $< -o $@

ising: setup $(OBJ_DIR)/ising.o 

$(OBJ_DIR)/ising.o: $(SRC_DIR)/ising.f90 $(OBJ_DIR)/functions.o 
	mpif90 $(FFLAGS)  -c $< -o $@ -I$(INCLUDE_DIR)

functions: setup $(OBJ_DIR)/functions.o

$(OBJ_DIR)/functions.o: $(SRC_DIR)/module_ising.f90
	mpif90 $(FFLAGS) -J $(INCLUDE_DIR) -c $^ -o $@

setup:
	@mkdir -p $(OBJ_DIR)
	@mkdir -p $(BIN_DIR)
	#@mkdir -p $(INCLUDE_DIR)

test:
	mpirun -np 2 $(BIN_DIR)/ising_mpi.x

clean:
	rm -f *.dat fort*

cleanall:
	rm -f *.o *.out *.mod *.dat fort*
	rm -f -r $(OBJ_DIR) $(BIN_DIR)
