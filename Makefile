#!/usr/bin

SRC = ./src/
OBJ = ./obj/
MOD = ./mod/
FCC = /usr/bin/gfortran
GCC = /usr/bin/gcc
FLAGS = -J$(MOD) -O3 -Wall -std=gnu -c
LIBS = /usr/lib64
INCS = /usr/include
SOURCE = mrrmain_tux.F90
OBJECT = $(OBJ)$(SOURCE:.F90=.o)
EXEC = MRR4ADMI
RM = rm -f

$(EXEC): $(OBJECT) $(OBJ)read_RAW_Cfile.o $(OBJ)psg_netcdf.o $(OBJ)extras.o $(OBJ)var_init.o
	$(FCC) $^ -o $@ -L$(LIBS) -lnetcdff

$(OBJECT): $(SRC)$(SOURCE) $(MOD)variables.mod $(MOD)extras.mod $(MOD)psg_nc.mod
	$(FCC) $(FLAGS) -c $(SOURCE) -o $@ 

$(OBJ)read_RAW_Cfile.o: $(SRC)read_RAW_Cfile.c
	$(GCC) -Wall -O3 -c $^ -o $@

$(MOD)extras.mod: $(OBJ)extras.o
$(OBJ)extras.o: $(SRC)extras.F90 $(SRC)var_init.F90
	$(FCC) $(FLAGS) $< -o $@

$(MOD)psg_nc.mod: $(OBJ)psg_netcdf.o
$(OBJ)psg_netcdf.o: $(SRC)psg_netcdf.F90
	$(FCC) $(FLAGS) $^ -o $@ -I$(INCS) -lnetcdf.mod -L$(LIBS) -lnetcdff

$(MOD)variables.mod: $(OBJ)var_init.o
$(OBJ)var_init.o: $(SRC)var_init.F90
	$(FCC) $(FLAGS) $^ -o $@

clean:
	$(RM) $(MOD)*.mod
	$(RM) $(OBJ)*.o
