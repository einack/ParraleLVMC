
SRC_DIR=src
OBJ_DIR=objects
BIN_DIR=bin
INCLUDE_DIR=objects
FFLAGS= -std=f2003 -finline-functions -fcheck=all -I$(INCLUDE_DIR) -Wtabs -Wall -Wunused-variable #-DDEBUG 

#a_mpi.out: random.f90 module_ising.f90 ising.f90 
#	mpif90 -DDEBUG -o a_mpi.out random.f90 module_ising.f90 ising.f90

#a.out: random.f90 module_ising.f90 ising.f90 
#	gfortran -o a.out random.f90 module_ising.f90 ising.f90
ising_serial.x: setup $(BIN_DIR)/ising_serial.x 

ising_mpi.x: setup $(BIN_DIR)/ising_mpi.x 

$(BIN_DIR)/ising_serial.x: $(OBJ_DIR)/ising_serial.o $(OBJ_DIR)/functions_serial.o  
	gfortran  -o $@ $^

$(BIN_DIR)/ising_mpi.x: $(OBJ_DIR)/ising_par.o $(OBJ_DIR)/functions.o  
	mpif90  -o $@ $^
	
#random_serial: setup $(OBJ_DIR)/random_serial.o 

#$(OBJ_DIR)/random_serial.o: $(SRC_DIR)/random.f90  
#	gfortran $(FFLAGS)  -c $< -o $@

#random: setup $(OBJ_DIR)/random.o 

#$(OBJ_DIR)/random.o: $(SRC_DIR)/random.f90  
#	mpif90 $(FFLAGS)  -c $< -o $@

ising_serial: setup $(OBJ_DIR)/ising_serial.o 

$(OBJ_DIR)/ising_serial.o: $(SRC_DIR)/ising_serial.f90 $(OBJ_DIR)/functions_serial.o 
	gfortran $(FFLAGS)  -c $< -o $@ -I$(INCLUDE_DIR)

ising_par: setup $(OBJ_DIR)/ising_par.o 

$(OBJ_DIR)/ising_par.o: $(SRC_DIR)/ising_par.f90 $(OBJ_DIR)/functions.o 
	mpif90 $(FFLAGS)  -c $< -o $@ -I$(INCLUDE_DIR)

functions_serial: setup $(OBJ_DIR)/functions_serial.o

$(OBJ_DIR)/functions_serial.o: $(SRC_DIR)/module_ising.f90
	gfortran $(FFLAGS) -J $(INCLUDE_DIR) -c $^ -o $@


functions: setup $(OBJ_DIR)/functions.o

$(OBJ_DIR)/functions.o: $(SRC_DIR)/module_ising.f90
	mpif90 $(FFLAGS) -J $(INCLUDE_DIR) -c $^ -o $@

setup:
	@mkdir -p $(OBJ_DIR)
	@mkdir -p $(BIN_DIR)
	#@mkdir -p $(INCLUDE_DIR)

test:
	./$(BIN_DIR)/ising_serial.x
	#mpirun -np 2 $(BIN_DIR)/ising_mpi.x

test_mpi:
	mpirun -np 2 $(BIN_DIR)/ising_mpi.x


run_ising_serial.x:
	 ./$(BIN_DIR)/ising_serial.x

clean:
	rm -f *.dat fort*

cleanall:
	rm -f *.o *.out *.mod *.dat fort*
	rm -f -r $(OBJ_DIR) $(BIN_DIR)
