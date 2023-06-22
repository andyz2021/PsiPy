ifeq ($(TARGET),)
  TARGET = mapfl
endif

ifeq ($(PRECISION),)
  PRECISION = r8
endif

PREFIX = /opt/psi/gnu/corhel/tools/ps_fortran

FC = gfortran-mp-11

FC_OPTS = -O3 -march=native -fcray-pointer -fallow-argument-mismatch -fallow-invalid-boz
LD_OPTS =

FC_OMP_OPTS = -fopenmp

# Change this to whatever it is on local computer to put object/mod files into that directory
OBJ_DIR = /Users/Andy/Desktop/PSI/PsiPy/psipy/fortran/getpb/obj_files

FC_LIB_OPTS = -fPIC

HDF5_INC_FLAGS = -I/opt/psi/gnu/ext_deps/deps/hdf5/include -I/opt/psi/gnu/ext_deps/deps/zlib/include
HDF5_LIB_FLAGS = -L/opt/psi/gnu/ext_deps/deps/hdf5/lib -L/opt/psi/gnu/ext_deps/deps/zlib/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lz

HDF4_INC_FLAGS = -I/opt/psi/gnu/ext_deps/deps/hdf4/include -I/opt/psi/gnu/ext_deps/deps/jpeg/include -I/opt/psi/gnu/ext_deps/deps/zlib/include
HDF4_LIB_FLAGS = -L/opt/psi/gnu/ext_deps/deps/hdf4/lib -L/opt/psi/gnu/ext_deps/deps/jpeg/lib -L/opt/psi/gnu/ext_deps/deps/zlib/lib -lmfhdf -ldf -ljpeg -lz

ZM_NT        = number_types_$(PRECISION)
ZM_PARSE     = zm_parse_$(PRECISION)
ZM_SDS       = zm_sds_$(PRECISION)
ZM_LIB_FLAGS = -L$(PREFIX)/lib -l$(ZM_SDS) -l$(ZM_PARSE)

FCFLAGS = $(FC_OPTS) $(FC_OMP_OPTS) $(HDF5_INC_FLAGS) $(HDF4_INC_FLAGS) $(FC_LIB_OPTS)
LDFLAGS = $(LD_OPTS) $(ZM_LIB_FLAGS) $(HDF5_LIB_FLAGS) $(HDF4_LIB_FLAGS)

OBJS = \
 $(OBJ_DIR)/$(ZM_NT).o \
 $(OBJ_DIR)/zm_parse_modules.o \
 $(OBJ_DIR)/zm_sds_modules.o \
 $(OBJ_DIR)/magnetic_field_function_modules.o \
 $(OBJ_DIR)/magnetic_field_function.o \
 $(OBJ_DIR)/mapfl.o

all: $(OBJS)
	$(FC) $(FCFLAGS) $(OBJS) $(LDFLAGS) -o $(TARGET)

mapflpy: mapflpy.f90 $(OBJS)
	f2py -c mapflpy.f90 $(OBJS) -m mapflpy --fcompiler=$(FC) --f90exec=$(FC) --f77exec=$(FC) --f90flags='$(FCFLAGS)' -lgomp $(LDFLAGS)
	mv mapflpy*.so $(PREFIX)/../ps_python/mapflpy.so

clean:
	rm -f $(TARGET)
	rm -f $(OBJ_DIR)/*.mod $(OBJ_DIR)/*.o $(OBJ_DIR)/*.so

distclean:
	make clean
	rm -f Makefile

install:
	mv $(TARGET) $(PREFIX)/bin/$(TARGET)

$(OBJ_DIR)/$(ZM_NT).o: $(PREFIX)/include/$(ZM_NT).f
	$(FC) -c $(FCFLAGS) $< -o $(OBJ_DIR)/$(ZM_NT).o

$(OBJ_DIR)/zm_parse_modules.o: $(PREFIX)/include/zm_parse_modules.f
	$(FC) -c $(FCFLAGS) $< -o $(OBJ_DIR)/zm_parse_modules.o

$(OBJ_DIR)/zm_sds_modules.o: $(PREFIX)/include/zm_sds_modules.f
	$(FC) -c $(FCFLAGS) $< -o $(OBJ_DIR)/zm_sds_modules.o

$(OBJ_DIR)/magnetic_field_function_modules.o: magnetic_field_function_modules.f90
	$(FC) -c $(FCFLAGS) $< -o $(OBJ_DIR)/magnetic_field_function_modules.o

$(OBJ_DIR)/magnetic_field_function.o: magnetic_field_function.f90
	$(FC) -c $(FCFLAGS) $< -o $(OBJ_DIR)/magnetic_field_function.o

$(OBJ_DIR)/mapfl.o: mapfl.f90
	$(FC) -c $(FCFLAGS) $< -o $(OBJ_DIR)/mapfl.o
