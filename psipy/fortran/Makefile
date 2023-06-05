ifeq ($(TARGET),)
  TARGET = getpb
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
OBJ_DIR = /Users/Andy/Desktop/PSI/PsiPy/psipy/fortran/obj_files


HDF5_INC_FLAGS = -I/opt/psi/gnu/ext_deps/deps/hdf5/include -I/opt/psi/gnu/ext_deps/deps/zlib/include
HDF5_LIB_FLAGS = -L/opt/psi/gnu/ext_deps/deps/hdf5/lib -L/opt/psi/gnu/ext_deps/deps/zlib/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lz

HDF4_INC_FLAGS = -I/opt/psi/gnu/ext_deps/deps/hdf4/include -I/opt/psi/gnu/ext_deps/deps/jpeg/include -I/opt/psi/gnu/ext_deps/deps/zlib/include
HDF4_LIB_FLAGS = -L/opt/psi/gnu/ext_deps/deps/hdf4/lib -L/opt/psi/gnu/ext_deps/deps/jpeg/lib -L/opt/psi/gnu/ext_deps/deps/zlib/lib -lmfhdf -ldf -ljpeg -lz

ZM_NT        = number_types_$(PRECISION)
ZM_PARSE     = zm_parse_$(PRECISION)
ZM_SDS       = zm_sds_$(PRECISION)
ZM_SPLINE    = zm_spline_$(PRECISION)
ZM_LIB_FLAGS = -L$(PREFIX)/lib -l$(ZM_SDS) -l$(ZM_PARSE) -l$(ZM_SPLINE)

FCFLAGS = $(FC_OPTS) $(FC_OMP_OPTS) $(HDF5_INC_FLAGS) $(HDF4_INC_FLAGS)
LDFLAGS = $(LD_OPTS) $(ZM_LIB_FLAGS) $(HDF5_LIB_FLAGS) $(HDF4_LIB_FLAGS)

OBJS = \
 $(OBJ_DIR)/$(ZM_NT).o \
 $(OBJ_DIR)/zm_parse_modules.o \
 $(OBJ_DIR)/zm_sds_modules.o \
 $(OBJ_DIR)/zm_spline_modules.o \
 $(OBJ_DIR)/mas_fields.o
all: $(OBJS)
	f2py -c getpb.f90 -I$(OBJ_DIR) --fcompiler=gfortran --f90flags='-fopenmp' -lgomp $(HDF5_INC_FLAGS) $(HDF4_INC_FLAGS) $(OBJS) $(LDFLAGS) -m getpb

clean:
	rm -f $(TARGET)
	rm -f *.mod *.o
	rm -f $(OBJ_DIR)/*.mod $(OBJ_DIR)/*.o

distclean:
	make clean
	rm -f Makefile

install:
	mv $(TARGET) $(PREFIX)/bin/$(TARGET)

$(OBJ_DIR)/$(ZM_NT).o: $(PREFIX)/include/$(ZM_NT).f
	$(FC) -J $(OBJ_DIR) -c $(FCFLAGS) $< -o $(OBJ_DIR)/$(ZM_NT).o

$(OBJ_DIR)/zm_parse_modules.o: $(PREFIX)/include/zm_parse_modules.f
	$(FC) -J $(OBJ_DIR) -c $(FCFLAGS) $< -o $(OBJ_DIR)/zm_parse_modules.o

$(OBJ_DIR)/zm_sds_modules.o: $(PREFIX)/include/zm_sds_modules.f
	$(FC) -J $(OBJ_DIR) -c $(FCFLAGS) $< -o $(OBJ_DIR)/zm_sds_modules.o

$(OBJ_DIR)/zm_spline_modules.o: $(PREFIX)/include/zm_spline_modules.f
	$(FC) -J $(OBJ_DIR) -c $(FCFLAGS) $< -o $(OBJ_DIR)/zm_spline_modules.o

$(OBJ_DIR)/mas_fields.o: mas_fields.f90
	$(FC) -J $(OBJ_DIR) -c $(FCFLAGS) $< -o $(OBJ_DIR)/mas_fields.o
