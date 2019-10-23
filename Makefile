
SRC_DIR=src
OBJ_DIR=objects
BIN_DIR=bin
INCLUDE_DIR=objects
FFLAGS= -g -std=f2003 -finline-functions -fcheck=all -I$(INCLUDE_DIR) -Wtabs -Wall -Wunused-variable  #-DDEBUG 
FFLAGSOMP= -g -std=f2003 -finline-functions -fcheck=all -I$(INCLUDE_DIR) -Wtabs -Wall -Wunused-variable -fopenmp #-DDEBUG 


ising_omp.x: setup $(BIN_DIR)/ising_omp.x 

ising_serial.x: setup $(BIN_DIR)/ising_serial.x 

ising_par.x: setup $(BIN_DIR)/ising_par.x 

$(BIN_DIR)/ising_omp.x: $(OBJ_DIR)/ising_omp.o $(OBJ_DIR)/functions_omp.o  
	gfortran $(FFLAGSOMP) -o $@ $^

$(BIN_DIR)/ising_serial.x: $(OBJ_DIR)/ising_serial.o $(OBJ_DIR)/functions_serial.o  
	gfortran  -o $@ $^

$(BIN_DIR)/ising_par.x: $(OBJ_DIR)/ising_par.o $(OBJ_DIR)/functions.o  
	mpif90  -o $@ $^
	
##############################################################################

ising_omp: setup $(OBJ_DIR)/ising_omp.o 

$(OBJ_DIR)/ising_omp.o: $(SRC_DIR)/ising_omp.f90 $(OBJ_DIR)/functions_omp.o 
	gfortran $(FFLAGSOMP)  -c $< -o $@ -I$(INCLUDE_DIR)

ising_serial: setup $(OBJ_DIR)/ising_serial.o 

$(OBJ_DIR)/ising_serial.o: $(SRC_DIR)/ising_serial.f90 $(OBJ_DIR)/functions_serial.o 
	gfortran $(FFLAGS)  -c $< -o $@ -I$(INCLUDE_DIR)

ising_par: setup $(OBJ_DIR)/ising_par.o 

$(OBJ_DIR)/ising_par.o: $(SRC_DIR)/ising_par.f90 $(OBJ_DIR)/functions.o 
	mpif90 $(FFLAGS)  -c $< -o $@ -I$(INCLUDE_DIR)

##############################################################################

functions_omp: setup $(OBJ_DIR)/functions_omp.o

$(OBJ_DIR)/functions_omp.o: $(SRC_DIR)/module_ising_omp.f90
	gfortran $(FFLAGSOMP) -J $(INCLUDE_DIR) -c $^ -o $@

functions_serial: setup $(OBJ_DIR)/functions_serial.o

$(OBJ_DIR)/functions_serial.o: $(SRC_DIR)/module_ising_serial.f90
	gfortran $(FFLAGS) -J $(INCLUDE_DIR) -c $^ -o $@

functions: setup $(OBJ_DIR)/functions.o

$(OBJ_DIR)/functions.o: $(SRC_DIR)/module_ising.f90
	mpif90 $(FFLAGS) -J $(INCLUDE_DIR) -c $^ -o $@

setup:
	@mkdir -p $(OBJ_DIR)
	@mkdir -p $(BIN_DIR)
	#@mkdir -p $(INCLUDE_DIR)

test_serial:
	@mkdir -p results_serial
	./$(BIN_DIR)/ising_serial.x < parameters.txt 
	#mpirun -np 2 $(BIN_DIR)/ising_mpi.x

test_omp:
	@mkdir -p results_omp
	./$(BIN_DIR)/ising_omp.x < parameters.txt  

test_par:
	@mkdir -p results_par  
	mpirun -np 2 $(BIN_DIR)/ising_par.x < parameters.txt  


clean:
	rm -f *.dat fort*

cleanall:
	rm -f *.o *.out *.mod *.dat fort*
	rm -f -r $(OBJ_DIR) $(BIN_DIR)
	rm -rf results_*
