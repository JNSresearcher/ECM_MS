F90=gfortran
CC = gcc
F=ECM_MS
TARGET =$(F)

#  -Werror=unused-parameter -Werror=unused-variable -Werror=unused-function -Wno-maybe-uninitialized 
FCOPT1 = -O2 -finit-local-zero -fopenmp -ffree-line-length-none
FCOPT2 = -O2
Lsw = -lz -lgfortran -lstdc++ -lgomp -ldl -lpthread
SRC = uncompress m_fparser m_vxc2data vxc2data  $(F)
OBJ = $(addsuffix .o, $(SRC) )
.SUFFIXES: .o .mod .f90 .F90 .f .cpp
%.o: %.f90
	$(F90) -o $@ -c $< $(FCOPT1) 
%.o: %.cpp
	$(CC) -o $@ -c $< $(FCOPT2) 
$(TARGET): $(OBJ)
	$(F90)  -o $@  $^  $(Lsw)
# Cleaning everything
clean:
	del *.o; del *.mod; del *.exe



