# Build, test and install fortran_shapefuncs without fpm.
#
#   make                  build the static library
#   make test             build and run the unit tests
#   make example          build and run the example
#   make install          install library, modules and pkg-config file
#   make uninstall        remove the installed files
#   make clean            remove the build directory
#
# Variables (override on the command line, e.g. `make FC=ifx FFLAGS=-O3`):
#   FC, FFLAGS            compiler and flags (gfortran, ifx, nvfortran)
#   PREFIX, DESTDIR       install location; DESTDIR stages a package
#   CUBATURES             path to fortran_cubatures' src/cubatures.f90;
#                         cloned from CUBATURES_URL when not given
#
# The library bundles the cubatures module, so programs link only
# -lfortran_shapefuncs.

ifeq ($(origin FC),default)
  FC = gfortran
endif
FFLAGS  ?= -O2
PREFIX  ?= /usr/local
LIBDIR  ?= $(PREFIX)/lib
INCDIR  ?= $(PREFIX)/include
BUILD   ?= build/make
VERSION  = 0.1.0

CUBATURES_URL ?= https://github.com/willklausler/fortran_cubatures
CUBATURES     ?= $(BUILD)/fortran_cubatures/src/cubatures.f90

ifneq (,$(findstring gfortran,$(FC)))
  MODFLAG = -J
else
  MODFLAG = -module
endif

LIB  = $(BUILD)/libfortran_shapefuncs.a
OBJS = $(BUILD)/cubatures.o $(BUILD)/shapefuncs.o
MODS = $(BUILD)/cubatures.mod $(BUILD)/shapefuncs.mod
PC   = $(BUILD)/fortran_shapefuncs.pc

.PHONY: all test example install uninstall clean

all: $(LIB) $(PC)

$(BUILD):
	mkdir -p $@

$(BUILD)/fortran_cubatures/src/cubatures.f90: | $(BUILD)
	git clone --depth 1 $(CUBATURES_URL) $(BUILD)/fortran_cubatures

$(BUILD)/cubatures.o: $(CUBATURES) | $(BUILD)
	$(FC) $(FFLAGS) $(MODFLAG) $(BUILD) -c $< -o $@

$(BUILD)/shapefuncs.o: src/shapefuncs.f90 $(BUILD)/cubatures.o
	$(FC) $(FFLAGS) $(MODFLAG) $(BUILD) -I$(BUILD) -c $< -o $@

$(LIB): $(OBJS)
	ar rcs $@ $^

$(PC): Makefile | $(BUILD)
	printf '%s\n' 'prefix=$(PREFIX)' 'libdir=$(LIBDIR)' 'includedir=$(INCDIR)' '' \
	  'Name: fortran_shapefuncs' \
	  'Description: Finite element shape functions, including mapped infinite elements' \
	  'Version: $(VERSION)' \
	  'Libs: -L$${libdir} -lfortran_shapefuncs' \
	  'Cflags: -I$${includedir}' > $@

$(BUILD)/%: test/%.f90 $(LIB)
	$(FC) $(FFLAGS) -I$(BUILD) $< $(LIB) -o $@

$(BUILD)/%: example/%.f90 $(LIB)
	$(FC) $(FFLAGS) -I$(BUILD) $< $(LIB) -o $@

test: $(BUILD)/shapefuncs_test
	$<

example: $(BUILD)/shapefuncs_example
	$<

install: all
	install -d $(DESTDIR)$(LIBDIR) $(DESTDIR)$(LIBDIR)/pkgconfig $(DESTDIR)$(INCDIR)
	install -m 644 $(LIB) $(DESTDIR)$(LIBDIR)
	install -m 644 $(MODS) $(DESTDIR)$(INCDIR)
	install -m 644 $(PC) $(DESTDIR)$(LIBDIR)/pkgconfig

uninstall:
	rm -f $(DESTDIR)$(LIBDIR)/libfortran_shapefuncs.a \
	      $(DESTDIR)$(LIBDIR)/pkgconfig/fortran_shapefuncs.pc \
	      $(addprefix $(DESTDIR)$(INCDIR)/,$(notdir $(MODS)))

clean:
	rm -rf $(BUILD)
