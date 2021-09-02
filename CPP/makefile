
# the compiler to use
CC  := g++
CXX := g++-9
XLF := gfortran
#XLF_77 := gfortran -std=legacy -O2 -fdefault-real-8 -fdefault-double-8 -cpp -ffree-form -ffree-line-length-1000 -fPIC

TARGET    := run
#H_FILES   := $(shell find ./ -regex [^\#]*\\.h$)
SRC_FILES := $(shell find ./ -regex [^\#]*\\.cpp$)
OBJ_FILES := $(SRC_FILES:.cpp=.o)
DEP_FILES := $(SRC_FILES:.cpp=.d)

SRCFORT_90_FILES := $(shell find ./ -regex [^\#]*\\.f90$)
SRCFORT_77_FILES := $(shell find ./ -regex [^\#]*\\.f$)
#OBJFORT_90_FILES := chiral-twobody-potentials.o chp-set.o nijmegen_interface.o #$(patsubst %.f90,%.o,$(shell find ./ -regex [^\#]*\\.f90$)) #$(SRCFORT_F90_FILES:./Potentials/chp/%.f90=%.o)
#OBJFORT_77_FILES := pnijm.o #$(patsubst %.f,%.o,$(shell find ./ -regex [^\#]*\\.f$))
OBJFORT_90_FILES := $(patsubst %.f90,%.o,$(shell find ./ -regex [^\#]*\\.f90$)) #$(SRCFORT_F90_FILES:./Potentials/chp/%.f90=%.o)
OBJFORT_77_FILES := $(patsubst %.f,%.o,$(shell find ./ -regex [^\#]*\\.f$))

#CPPFLAGS := -Wall -O2 -std=c++11 
#LDFLAGS := -lblas -llapack -lm


CPPFLAGS := -Wall -std=c++17 -O2 -ggdb -fopenmp
# Stop compilation on warning-encounter (good for neat code but not necessary)
#CPPFLAGS += -Werror
# Use MKL-types
CPPFLAGS += -DMKL_ILP64 -m64 -DMKL_Complex16="std::complex<double>" -DMKL_Complex8="std::complex<float>" -I$(MKLROOT)/include
FORTFLAGS_90 := -O2 -fdefault-real-8 -fdefault-double-8 -cpp -ffree-form -ffree-line-length-1000 -fPIC
FORTFLAGS_77 := -O2

#LDLIBS :=  -Wl,--no-as-needed -lgomp -lgsl -lgslcblas -lpthread -lm -ldl -lgfortran -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_serial_hl  -lhdf5_serial 
LDLIBS :=  -Wl,--no-as-needed -lgomp -lgsl -lpthread -lm -ldl -lgfortran -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_serial_hl  -lhdf5_serial 
# Use MKL-libraries
LDLIBS += -L$(MKLROOT)/lib/intel64 -lmkl_intel_ilp64 -lmkl_gnu_thread -lmkl_core

# Symbol explentations
# @ = recipe name (so $(TARGET) in this case)
# ^ = all dependencies (so ALL $(OBJ_FILES))
# < = first dependency (so first of OBJ_FILES)

# Include the .d files
-include $(DEP_FILES)

$(TARGET): $(OBJ_FILES) $(OBJFORT_90_FILES) $(OBJFORT_77_FILES)
	$(CC) $^ -o $@ $(LDFLAGS) $(LDLIBS) 

%.o: %.cpp Makefile
	$(CC) $(CPPFLAGS) -MMD -MP -c $< -o $@

#$(OBJFORT_90_FILES): $(SRCFORT_90_FILES)
#	$(XLF) $(FORTFLAGS_90) -c $^
#
#$(OBJFORT_77_FILES): $(SRCFORT_77_FILES)
#	$(XLF) $(FORTFLAGS_77) -c $^

%.o: %.f90
	$(XLF) $(FORTFLAGS_90) -c $< -o $@

%.o: %.f
	$(XLF) $(FORTFLAGS_77) -c $< -o $@


.PHONY: clean cleanall

# Clean everything but FORTRAN files
clean:
	rm -f $(TARGET) $(OBJ_FILES) $(DEP_FILES) *.o *.exe

# Clean everything
cleanall:
	rm -f $(TARGET) $(OBJ_FILES) $(OBJFORT_90_FILES) $(OBJFORT_77_FILES) $(DEP_FILES) *.o *.mod *.exe

