# Makefile for standalone build of the air EOS, Henry, water-vapour,
# properties, and PFLOTRAN-facing wrapper modules.
#
# Targets:
#   make                 -- build all test drivers
#   make test            -- run basic PR+Peneloux anchor test
#   make sweep           -- run P-T grid sweep (writes CSV)
#   make crosscheck      -- Python cross-validation against sweep CSV
#   make henry           -- run Henry's-law test
#   make watervapour     -- run water-vapour test
#   make properties      -- run viscosity + ideal-gas enthalpy test
#   make wrapper         -- run PFLOTRAN-wrapper module test
#   make all             -- build and run everything
#   make clean           -- remove build artifacts

FC      ?= gfortran
FFLAGS  ?= -O2 -g -cpp -Wall -Wextra -fimplicit-none -std=f2008 \
           -Wno-maybe-uninitialized -Wno-unused-variable \
           -Wno-unused-parameter -Wno-unused-dummy-argument
LDFLAGS ?=

SRC_EOS_MOD       = air_eos_pr_module.F90
SRC_HENRY_MOD     = air_henry_module.F90
SRC_WV_MOD        = air_watervapour_module.F90
SRC_PROP_MOD      = air_properties_module.F90
SRC_WRAP_MOD      = air_eos_module.F90
SRC_TEST          = test_air_eos_pr.F90
SRC_SWEEP         = test_air_eos_pr_sweep.F90
SRC_HENRY_TEST    = test_air_henry.F90
SRC_WV_TEST       = test_air_watervapour.F90
SRC_PROP_TEST     = test_air_properties.F90
SRC_WRAP_TEST     = test_air_eos_module.F90

OBJ_EOS_MOD       = $(SRC_EOS_MOD:.F90=.o)
OBJ_HENRY_MOD     = $(SRC_HENRY_MOD:.F90=.o)
OBJ_WV_MOD        = $(SRC_WV_MOD:.F90=.o)
OBJ_PROP_MOD      = $(SRC_PROP_MOD:.F90=.o)
OBJ_WRAP_MOD      = $(SRC_WRAP_MOD:.F90=.o)
OBJ_TEST          = $(SRC_TEST:.F90=.o)
OBJ_SWEEP         = $(SRC_SWEEP:.F90=.o)
OBJ_HENRY_TEST    = $(SRC_HENRY_TEST:.F90=.o)
OBJ_WV_TEST       = $(SRC_WV_TEST:.F90=.o)
OBJ_PROP_TEST     = $(SRC_PROP_TEST:.F90=.o)
OBJ_WRAP_TEST     = $(SRC_WRAP_TEST:.F90=.o)

EXE_TEST          = test_air_eos_pr
EXE_SWEEP         = test_air_eos_pr_sweep
EXE_HENRY         = test_air_henry
EXE_WV            = test_air_watervapour
EXE_PROP          = test_air_properties
EXE_WRAP          = test_air_eos_module

.PHONY: all test sweep crosscheck henry watervapour properties wrapper clean

all: $(EXE_TEST) $(EXE_SWEEP) $(EXE_HENRY) $(EXE_WV) $(EXE_PROP) $(EXE_WRAP)

$(EXE_TEST): $(OBJ_EOS_MOD) $(OBJ_TEST)
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)

$(EXE_SWEEP): $(OBJ_EOS_MOD) $(OBJ_SWEEP)
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)

$(EXE_HENRY): $(OBJ_EOS_MOD) $(OBJ_HENRY_MOD) $(OBJ_HENRY_TEST)
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)

$(EXE_WV): $(OBJ_EOS_MOD) $(OBJ_WV_MOD) $(OBJ_WV_TEST)
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)

$(EXE_PROP): $(OBJ_EOS_MOD) $(OBJ_PROP_MOD) $(OBJ_PROP_TEST)
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)

$(EXE_WRAP): $(OBJ_EOS_MOD) $(OBJ_HENRY_MOD) $(OBJ_WV_MOD) $(OBJ_PROP_MOD) \
             $(OBJ_WRAP_MOD) $(OBJ_WRAP_TEST)
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)

$(OBJ_EOS_MOD): $(SRC_EOS_MOD)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_HENRY_MOD): $(SRC_HENRY_MOD) $(OBJ_EOS_MOD)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_WV_MOD): $(SRC_WV_MOD) $(OBJ_EOS_MOD)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_PROP_MOD): $(SRC_PROP_MOD) $(OBJ_EOS_MOD)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_WRAP_MOD): $(SRC_WRAP_MOD) $(OBJ_EOS_MOD) $(OBJ_HENRY_MOD) $(OBJ_WV_MOD) $(OBJ_PROP_MOD)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_TEST): $(SRC_TEST) $(OBJ_EOS_MOD)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_SWEEP): $(SRC_SWEEP) $(OBJ_EOS_MOD)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_HENRY_TEST): $(SRC_HENRY_TEST) $(OBJ_EOS_MOD) $(OBJ_HENRY_MOD)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_WV_TEST): $(SRC_WV_TEST) $(OBJ_EOS_MOD) $(OBJ_WV_MOD)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_PROP_TEST): $(SRC_PROP_TEST) $(OBJ_EOS_MOD) $(OBJ_PROP_MOD)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_WRAP_TEST): $(SRC_WRAP_TEST) $(OBJ_WRAP_MOD)
	$(FC) $(FFLAGS) -c $< -o $@

test: $(EXE_TEST)
	./$(EXE_TEST)

sweep: $(EXE_SWEEP)
	./$(EXE_SWEEP)

henry: $(EXE_HENRY)
	./$(EXE_HENRY)

watervapour: $(EXE_WV)
	./$(EXE_WV)

properties: $(EXE_PROP)
	./$(EXE_PROP)

wrapper: $(EXE_WRAP)
	./$(EXE_WRAP)

crosscheck: air_eos_sweep.csv
	python3 cross_validate.py air_eos_sweep.csv

air_eos_sweep.csv: $(EXE_SWEEP)
	./$(EXE_SWEEP)

clean:
	rm -f *.o *.mod $(EXE_TEST) $(EXE_SWEEP) $(EXE_HENRY) $(EXE_WV) \
	      $(EXE_PROP) $(EXE_WRAP) air_eos_sweep.csv
