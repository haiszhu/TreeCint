SHELL := /bin/zsh

# Resolve project root from this Makefile location (no hardcoded absolute path).
ROOT := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
LIBCINT_SRC := $(ROOT)/external/libcint
LIBCINT_BUILD := $(ROOT)/build/libcint
LIBCINT_INSTALL := $(ROOT)/install/libcint
LIBID_SRC := $(ROOT)/external/libid
DMK_SRC := $(ROOT)/external/dmk
LIBID_MW := $(ROOT)/src/libid/libid.mw
LIBID_MEX_SRC := $(ROOT)/src/libid/libid_mex.f90
TREECINT_WRAPPER := $(ROOT)/src/wrapper
TREECINT_UTILS := $(ROOT)/matlab
TREECINT_BUILD := $(ROOT)/build/treecint-mex
LIBID_MEX_BUILD := $(ROOT)/build/libid-mex
BDMK_BUILD := $(ROOT)/build/bdmk
BDMK_LIB := $(BDMK_BUILD)/libbdmk_ref.a
BDMK_MEX_BUILD := $(ROOT)/build/bdmk-mex
BDMK_MW := $(ROOT)/src/bdmk/bdmk.mw
BDMK_MEX_SRC := $(ROOT)/src/bdmk/bdmk_mex.f90
INT2BDMK_SCF_SRC := $(ROOT)/src/bdmk/bdmk_scf.f90
INT2BDMK_SCF_EXE := $(ROOT)/bin/int2-bdmk-scf
HDF5_PREFIX ?= /opt/homebrew/opt/hdf5
MATLAB_ROOT := $(shell ls -d /Applications/MATLAB_R*.app | sort | tail -n1)
MWRAP ?= ~/mwrap/mwrap

CMAKE ?= cmake
C_COMPILER ?= gcc-14
F_COMPILER ?= gfortran
AR ?= ar
USE_OPENMP ?= 0
OPENMP_FLAG :=
OPENMP_LIB := -lgomp
ifeq ($(USE_OPENMP),1)
OPENMP_FLAG := -fopenmp
endif

# Match the key build flavor from your modified libcint makefile.
LIBCINT_CFLAGS := -fPIC -O2 $(OPENMP_FLAG) -D_Float128=__float128
TREECINT_CFLAGS := -fPIC -O2 $(OPENMP_FLAG) -D_Float128=__float128
LIBID_FFLAGS ?= -O2 -w -framework accelerate -fallow-argument-mismatch
BDMK_FFLAGS ?= -fPIC -O3 -std=legacy -funroll-loops -w -fallow-argument-mismatch

BDMK_SOURCES := \
	$(DMK_SRC)/src/common/specialfunctions/besseljs3d.f \
	$(DMK_SRC)/src/common/specialfunctions/hank103.f \
	$(DMK_SRC)/src/common/specialfunctions/legeexps.f \
	$(DMK_SRC)/src/common/specialfunctions/chebexps.f \
	$(DMK_SRC)/src/common/prini_new.f \
	$(DMK_SRC)/src/common/fmmcommon2d.f \
	$(DMK_SRC)/src/common/lapack_f77.f \
	$(DMK_SRC)/src/common/cumsum.f \
	$(DMK_SRC)/src/common/hkrand.f \
	$(DMK_SRC)/src/common/dlaran.f \
	$(DMK_SRC)/src/common/voltab2d.f \
	$(DMK_SRC)/src/common/voltab3d.f \
	$(DMK_SRC)/src/common/polytens.f \
	$(DMK_SRC)/src/common/tree_data_routs.f \
	$(DMK_SRC)/src/common/tensor_prod_routs.f \
	$(DMK_SRC)/src/common/pts_tree.f \
	$(DMK_SRC)/src/common/tree_routs.f \
	$(DMK_SRC)/src/common/tree_vol_coeffs.f \
	$(DMK_SRC)/src/common/dmk_routs.f \
	$(DMK_SRC)/src/bdmk/sogapproximation/get_sognodes.f \
	$(DMK_SRC)/src/bdmk/sogapproximation/l2dsognodes.f \
	$(DMK_SRC)/src/bdmk/sogapproximation/l3dsognodes.f \
	$(DMK_SRC)/src/bdmk/sogapproximation/sl3dsognodes.f \
	$(DMK_SRC)/src/bdmk/sogapproximation/y2dsognodes.f \
	$(DMK_SRC)/src/bdmk/sogapproximation/y3dsognodes.f \
	$(DMK_SRC)/src/bdmk/bdmk_local_tables.f \
	$(DMK_SRC)/src/bdmk/bdmk_local.f \
	$(DMK_SRC)/src/bdmk/bdmk_pwterms.f \
	$(DMK_SRC)/src/bdmk/bdmk_pwrouts.f \
	$(DMK_SRC)/src/bdmk/boxfgt_md.f \
	$(DMK_SRC)/src/bdmk/bdmk.f

BDMK_OBJECTS := $(patsubst $(DMK_SRC)/src/%.f,$(BDMK_BUILD)/%.o,$(BDMK_SOURCES))

.PHONY: all help libcint-configure libcint-build libcint-install libcint-all libcint-clean libid-build libid-test libid-mex libid-clean bdmk-lib bdmk-mex bdmk-clean int2-bdmk-scf int2-bdmk-scf-clean treecint-mex treecint-clean test-eval-match test-treefun-libcint test-libid-id test-vmunu-parallel clean

all: treecint-mex bdmk-mex int2-bdmk-scf

help:
	@echo "Targets:"
	@echo "  make int2-bdmk-scf     # Build standalone Fortran BDMK executable (no MEX)"
	@echo "  make libcint-all       # configure + build + install upstream libcint"
	@echo "  make libcint-configure # CMake configure for upstream libcint"
	@echo "  make libcint-build     # Build upstream libcint"
	@echo "  make libcint-install   # Install headers/libs to install/libcint"
	@echo "  make libcint-clean     # Remove build/libcint"
	@echo "  make libid-build       # Build libid archive in external/libid"
	@echo "  make libid-test        # Run upstream libid tests"
	@echo "  make libid-mex         # Build iddr_aid_mwrap_mex via mwrap + Fortran API"
	@echo "  make libid-clean       # Clean libid build/test artifacts"
	@echo "  make bdmk-lib          # Build makefile-style BDMK static library from external/dmk Fortran sources"
	@echo "  make bdmk-mex          # Build low-level bdmk wrapper mex linked to libbdmk_ref.a"
	@echo "  make bdmk-clean        # Remove build/bdmk"
	@echo "  make treecint-mex      # Build TreeCint.mexmaca64 + MWrap stubs"
	@echo "  make test-eval-match   # Compare new TreeCint AO vs legacy evaluator"
	@echo "  make test-treefun-libcint # Treefun3 smoke test with libcint AO callback"
	@echo "  make treecint-clean    # Remove build/treecint-mex and TreeCint mex"

libcint-configure:
	@mkdir -p $(LIBCINT_BUILD) $(LIBCINT_INSTALL)
	$(CMAKE) -S $(LIBCINT_SRC) \
		-B $(LIBCINT_BUILD) \
		-DCMAKE_C_COMPILER=$(C_COMPILER) \
		-DCMAKE_BUILD_TYPE=Release \
		-DCMAKE_C_FLAGS="$(LIBCINT_CFLAGS)" \
		-DWITH_CINT2_INTERFACE=ON \
		-DWITH_FORTRAN=ON \
		-DBUILD_SHARED_LIBS=ON \
		-DCMAKE_INSTALL_PREFIX=$(LIBCINT_INSTALL)

libcint-build:
	$(CMAKE) --build $(LIBCINT_BUILD) -j8

libcint-install:
	$(CMAKE) --install $(LIBCINT_BUILD) --prefix $(LIBCINT_INSTALL)

libcint-all: libcint-configure libcint-build libcint-install

libcint-clean:
	rm -rf $(LIBCINT_BUILD)

libid-build:
	$(MAKE) -C $(LIBID_SRC) all FFLAGS="$(LIBID_FFLAGS)"

libid-test: libid-build
	@echo "libid tests already executed by libid-build (upstream all target)."

libid-mex: libid-build
	@mkdir -p $(LIBID_MEX_BUILD) $(TREECINT_UTILS)
	rm -f $(TREECINT_UTILS)/iddr_aid_mwrap_mex.mexmaca64
	cd $(LIBID_MEX_BUILD) && \
		$(MWRAP) -c99complex -mex libid_mex -mb -list $(LIBID_MW) && \
		$(MWRAP) -c99complex -mex libid_mex -c libid_mex.c $(LIBID_MW) && \
		mv -f ./*_mex.m $(TREECINT_UTILS)/
	$(F_COMPILER) -c $(LIBID_FFLAGS) \
		-I$(LIBID_SRC)/src \
		$(LIBID_MEX_SRC) \
		-o $(LIBID_MEX_BUILD)/libid_mex.o
		$(C_COMPILER) -shared $(TREECINT_CFLAGS) \
		-DMATLAB_MEX_FILE -DMWF77_UNDERSCORE1 -DMATLAB_DEFAULT_RELEASE=R2023b -DMX_COMPAT_32=0 \
		-I$(MATLAB_ROOT)/extern/include \
		$(LIBID_MEX_BUILD)/libid_mex.c \
		$(LIBID_MEX_BUILD)/libid_mex.o \
		-L$(MATLAB_ROOT)/bin/maca64 -lmx -lmex -lmat -lm \
		$(LIBID_SRC)/id_lib.a \
		-framework Accelerate \
		-lgfortran -lquadmath \
		-o $(TREECINT_UTILS)/libid_mex.mexmaca64

libid-clean:
	$(MAKE) -C $(LIBID_SRC) clean

bdmk-lib: $(BDMK_LIB)

$(BDMK_LIB): $(BDMK_OBJECTS)
	@mkdir -p $(BDMK_BUILD)
	$(AR) rcs $@ $(BDMK_OBJECTS)

$(BDMK_BUILD)/%.o: $(DMK_SRC)/src/%.f
	@mkdir -p $(dir $@)
	$(F_COMPILER) -c $(BDMK_FFLAGS) $< -o $@

bdmk-mex: bdmk-lib
	@mkdir -p $(BDMK_MEX_BUILD) $(TREECINT_UTILS)
	cd $(BDMK_MEX_BUILD) && \
		$(MWRAP) -c99complex -mex bdmk_mex -mb -list $(BDMK_MW) && \
		$(MWRAP) -c99complex -mex bdmk_mex -c bdmk_mex.c $(BDMK_MW) && \
		mv -f ./*_mex.m $(TREECINT_UTILS)/
	$(F_COMPILER) -c $(BDMK_FFLAGS) $(BDMK_MEX_SRC) -o $(BDMK_MEX_BUILD)/bdmk_mex.o
	$(C_COMPILER) -shared $(TREECINT_CFLAGS) \
		-DMATLAB_MEX_FILE -DMWF77_UNDERSCORE1 -DMATLAB_DEFAULT_RELEASE=R2023b -DMX_COMPAT_32=0 \
		-I$(MATLAB_ROOT)/extern/include \
		$(BDMK_MEX_BUILD)/bdmk_mex.c \
		$(BDMK_MEX_BUILD)/bdmk_mex.o \
		-L$(MATLAB_ROOT)/bin/maca64 -lmx -lmex -lmat -lm \
		$(BDMK_LIB) \
		-framework Accelerate \
		-lgfortran -lquadmath \
		-o $(TREECINT_UTILS)/bdmk_mex.mexmaca64

bdmk-clean:
	rm -rf $(BDMK_BUILD)
	rm -rf $(BDMK_MEX_BUILD)
	rm -f $(TREECINT_UTILS)/bdmk_mex.mexmaca64
	rm -f $(TREECINT_UTILS)/bdmk_wrap_mwrap_mex.m
	rm -f $(TREECINT_UTILS)/bdmk_eval_mex.m

# Standalone Fortran executable: reads treefun+isdf HDF5, runs OMP-parallel BDMK,
# assembles Vmunu via a single DGEMM, writes Vmunu HDF5. No MEX/MATLAB dependency.
int2-bdmk-scf: bdmk-lib
	@mkdir -p $(dir $(INT2BDMK_SCF_EXE))
	$(F_COMPILER) -O3 -fopenmp -fallow-argument-mismatch -w \
		-I$(HDF5_PREFIX)/include \
		$(INT2BDMK_SCF_SRC) \
		$(BDMK_LIB) \
		-L$(HDF5_PREFIX)/lib -lhdf5_fortran -lhdf5 \
		-framework Accelerate \
		-lm \
		-o $(INT2BDMK_SCF_EXE)
	@echo "Built: $(INT2BDMK_SCF_EXE)"

int2-bdmk-scf-clean:
	rm -f $(INT2BDMK_SCF_EXE)

test-vmunu-parallel: int2-bdmk-scf treecint-mex
	KMP_DUPLICATE_LIB_OK=TRUE OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 \
		$(MATLAB_ROOT)/bin/matlab -nojvm -nodesktop -nosplash \
		-r "cd('$(ROOT)/test'); h2o_ccpvdz_vmunu_vijkl_parallel; exit"

treecint-mex: libcint-all libid-test libid-mex
	@mkdir -p $(TREECINT_BUILD) $(TREECINT_UTILS)
	cd $(TREECINT_BUILD) && \
		$(MWRAP) -c99complex -mex TreeCint -mb -list $(TREECINT_WRAPPER)/gateway.mw && \
		$(MWRAP) -c99complex -mex TreeCint -c TreeCint_gateway.c $(TREECINT_WRAPPER)/gateway.mw && \
		mv -f ./*_mex.m $(TREECINT_UTILS)/
	$(C_COMPILER) -shared $(TREECINT_CFLAGS) \
		-Wno-implicit-function-declaration \
		-DMATLAB_MEX_FILE -DMATLAB_DEFAULT_RELEASE=R2023b -DMX_COMPAT_32=0 \
		-I$(MATLAB_ROOT)/extern/include \
		-I$(TREECINT_WRAPPER) -I$(TREECINT_WRAPPER)/gto -I$(TREECINT_WRAPPER)/np_helper \
		-I$(LIBCINT_SRC)/src -I$(LIBCINT_BUILD)/src -I$(LIBCINT_BUILD)/include -I$(LIBCINT_INSTALL)/include \
		$(TREECINT_BUILD)/TreeCint_gateway.c \
		$(TREECINT_WRAPPER)/cgto.c \
		$(TREECINT_WRAPPER)/gto/fill_int2c.c \
		$(TREECINT_WRAPPER)/gto/fill_nr_3c.c \
		$(TREECINT_WRAPPER)/gto/fill_r_3c.c \
		$(TREECINT_WRAPPER)/gto/fill_int2e.c \
		$(TREECINT_WRAPPER)/gto/fill_r_4c.c \
		$(TREECINT_WRAPPER)/gto/ft_ao.c \
		$(TREECINT_WRAPPER)/gto/ft_ao_deriv.c \
		$(TREECINT_WRAPPER)/gto/fill_grids_int2c.c \
		$(TREECINT_WRAPPER)/gto/grid_ao_drv.c \
		$(TREECINT_WRAPPER)/gto/deriv1.c \
		$(TREECINT_WRAPPER)/gto/deriv2.c \
		$(TREECINT_WRAPPER)/gto/nr_ecp.c \
		$(TREECINT_WRAPPER)/gto/nr_ecp_deriv.c \
		$(TREECINT_WRAPPER)/gto/autocode/auto_eval1.c \
		$(TREECINT_WRAPPER)/np_helper/transpose.c \
		$(TREECINT_WRAPPER)/np_helper/pack_tril.c \
		$(TREECINT_WRAPPER)/np_helper/npdot.c \
		$(TREECINT_WRAPPER)/np_helper/condense.c \
		$(TREECINT_WRAPPER)/np_helper/omp_reduce.c \
		$(TREECINT_WRAPPER)/np_helper/np_helper.c \
		-L$(MATLAB_ROOT)/bin/maca64 -lmx -lmex -lmat -lm \
		-L$(LIBCINT_INSTALL)/lib -lcint \
		-Wl,-rpath,$(LIBCINT_INSTALL)/lib \
		-framework Accelerate \
		-lstdc++ $(OPENMP_LIB) -lquadmath \
		-o $(TREECINT_UTILS)/TreeCint.mexmaca64

treecint-clean:
	rm -rf $(TREECINT_BUILD)
	rm -f $(TREECINT_UTILS)/TreeCint.mexmaca64
	rm -rf $(LIBID_MEX_BUILD)
	rm -f $(TREECINT_UTILS)/libid_mex.mexmaca64
	rm -f $(TREECINT_UTILS)/iddr_aid_mwrap_mex.mexmaca64

test-eval-match:
	KMP_DUPLICATE_LIB_OK=TRUE OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
		$(MATLAB_ROOT)/bin/matlab -nojvm -nodesktop -nosplash -r "cd('$(ROOT)/test'); compare_saved_eval; exit"

test-treefun-libcint: treecint-mex
	KMP_DUPLICATE_LIB_OK=TRUE OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
		$(MATLAB_ROOT)/bin/matlab -nojvm -nodesktop -nosplash -r "cd('$(ROOT)/test'); test_treefun_with_libcint; exit"

test-libid-id: treecint-mex
	KMP_DUPLICATE_LIB_OK=TRUE OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
		$(MATLAB_ROOT)/bin/matlab -nojvm -nodesktop -nosplash -r "cd('$(ROOT)/test'); test_libid_iddr_mex; exit"

clean: libcint-clean treecint-clean
clean: libid-clean
clean: bdmk-clean int2-bdmk-scf-clean
