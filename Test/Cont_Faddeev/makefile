# the compiler to use
CC  := g++ 
XLF := gfortran -O2 -fdefault-real-8 -fdefault-double-8 -cpp -ffree-form -ffree-line-length-1000 -fPIC

TARGET    := run
#H_FILES   := $(shell find ./ -regex [^\#]*\\.h$)
SRC_FILES := $(shell find ./ -regex [^\#]*\\.cpp$)
OBJ_FILES := $(SRC_FILES:.cpp=.o)
DEP_FILES := $(SRC_FILES:.cpp=.d)

SRCFORT_FILES := $(shell find ./ -regex [^\#]*\\.f90$)
OBJFORT_FILES := $(SRCFORT_FILES:.f90=.o)

#CPPFLAGS := -Wall -O2 -std=c++11 
#LDFLAGS := -lblas -llapack -lm


CPPFLAGS := -Wall -std=c++17 -O2 -ggdb
# Stop compilation on warning-encounter (good for neat code but not necessary)
#CPPFLAGS += -Werror
# Use MKL-types
#CPPFLAGS += -DMKL_ILP64 -m64 -DMKL_Complex16="std::complex<double>" -I$(MKLROOT)/include

LDLIBS :=  -Wl,--no-as-needed -lgomp -lpthread -lm -ldl -lgfortran
# Use MKL-libraries
#LDLIBS += -L$(MKLROOT)/lib/intel64 -lmkl_intel_ilp64 -lmkl_gnu_thread -lmkl_core

# Symbol explentations
# @ = recipe name (so $(TARGET) in this case)
# ^ = all dependencies (so ALL $(OBJ_FILES))
# < = first dependency (so first of OBJ_FILES)

# Include the .d files
-include $(DEP_FILES)

$(TARGET): $(OBJ_FILES) $(OBJFORT_FILES)
	$(CC) $^ -o $@ $(LDFLAGS) $(LDLIBS) 

%.o: %.cpp Makefile
	$(CC) $(CPPFLAGS) -MMD -MP -c $< -o $@

%.o: %.f90
	$(XLF) -c $< -o $@


.PHONY: clean
clean:
	rm -f $(TARGET) $(OBJ_FILES) $(OBJFORT_FILES) $(DEP_FILES) *.o *.mod *.exe
