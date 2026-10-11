# make(1) build, an alternative to fpm.
#   make                  library, tests, examples and app
#   make test             run the tests
#   make FC=ifx FFLAGS=-O2
#   make install PREFIX=$HOME/.local
#   make doc              HTML documentation in build/doc (needs ford)
#   make CUBATURES=../fortran_cubatures/src/cubatures.f90
#                         use a local fortran_cubatures instead of cloning it
#
# The library bundles the cubatures module, so programs link only
# -lfortran_shapefuncs.

ifeq ($(origin FC),default)
FC = gfortran
endif
FFLAGS ?= -O2
PREFIX ?= /usr/local

B := build/make
CUBATURES_TAG ?= v0.2.0
CUBATURES ?= $(B)/fortran_cubatures/src/cubatures.f90
LIB := $(B)/libfortran_shapefuncs.a
DEMOS := $(patsubst example/%.f90,$(B)/%,$(wildcard example/demo_*.f90))
PROGS := $(B)/test_shapefuncs $(DEMOS) $(B)/interpolation_convergence

.PHONY: all test install doc clean

all: $(LIB) $(PROGS)

$(B)/fortran_cubatures/src/cubatures.f90:
	git clone --depth 1 --branch $(CUBATURES_TAG) https://github.com/willklausler/fortran_cubatures $(B)/fortran_cubatures

# Compiling inside $(B) puts the .mod files there for every compiler.
$(LIB): src/shapefuncs.f90 $(CUBATURES)
	mkdir -p $(B)
	cd $(B) && $(FC) $(FFLAGS) -c $(abspath $(CUBATURES)) -o cubatures.o
	cd $(B) && $(FC) $(FFLAGS) -c $(CURDIR)/src/shapefuncs.f90 -o shapefuncs.o
	ar rcs $@ $(B)/cubatures.o $(B)/shapefuncs.o

LINK = $(FC) $(FFLAGS) -I$(B) $< $(LIB) -o $@
$(B)/test_shapefuncs: test/test_shapefuncs.f90 $(LIB)
	$(LINK)
$(B)/demo_%: example/demo_%.f90 $(LIB)
	$(LINK)
$(B)/interpolation_convergence: app/main.f90 $(LIB)
	$(LINK)

test: $(B)/test_shapefuncs
	$<

install: $(LIB)
	install -d $(PREFIX)/lib $(PREFIX)/include $(PREFIX)/share/man/man3
	install -m 644 $(LIB) $(PREFIX)/lib
	install -m 644 $(B)/*.mod $(PREFIX)/include
	install -m 644 man/shapefuncs.3 $(PREFIX)/share/man/man3

# The user manual page is llms.txt, so there is one reference to maintain.
doc:
	mkdir -p build/pages
	{ printf -- '---\ntitle: User manual\n---\n\n'; sed 1d llms.txt; } > build/pages/index.md
	ford ford.md
	cp llms.txt build/doc/

clean:
	rm -rf $(B)
