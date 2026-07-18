# -------------------------------------------------------
# KSROP Makefile (Unix/Linux/macOS; FC=gfortran by default)
# Usage:
#   make            build driver_KS (propagator)
#   make tools      build tle2opm (TLE-to-OPM converter)
#   make tests      build all test binaries
#   make run        run propagator
#   make test       build everything + run test_all.sh (lint + all tests)
#   make clean      remove build artefacts
# -------------------------------------------------------

FC      = gfortran
FFLAGS  = -O2 -Wall

SRC     = driver_KS.F Subrouts.F Legendre.F TLEread.F
TARGET  = driver_KS

TLE2OPM_SRC = tle2opm.F Subrouts.F TLEread.F Legendre.F
TLE2OPM_BIN = tle2opm

TEST_SRC  = test_subrouts.F Subrouts.F Legendre.F
TEST_BIN  = test_subrouts

TEST_BUGS_SRC = test_bugs.F Subrouts.F Legendre.F
TEST_BUGS_BIN = test_bugs

TEST_TLE_SRC  = test_tle.F Subrouts.F TLEread.F Legendre.F
TEST_TLE_BIN  = test_tle

TEST_TLE2SV_SRC = test_tle2sv.F Subrouts.F TLEread.F Legendre.F
TEST_TLE2SV_BIN = test_tle2sv

TEST_TLE2OPM_SRC = test_tle2opm.F Subrouts.F TLEread.F Legendre.F
TEST_TLE2OPM_BIN = test_tle2opm

ALL_TEST_BINS = $(TEST_BIN) $(TEST_BUGS_BIN) $(TEST_TLE_BIN) \
                $(TEST_TLE2SV_BIN) $(TEST_TLE2OPM_BIN)

.PHONY: all tools tests run test clean

all: $(TARGET)

$(TARGET): $(SRC)
	$(FC) $(FFLAGS) $(SRC) -o $(TARGET)

tools: $(TLE2OPM_BIN)

$(TLE2OPM_BIN): $(TLE2OPM_SRC)
	$(FC) $(FFLAGS) $(TLE2OPM_SRC) -o $(TLE2OPM_BIN)

tests: $(ALL_TEST_BINS)

$(TEST_BIN): $(TEST_SRC)
	$(FC) $(FFLAGS) $(TEST_SRC) -o $(TEST_BIN)

$(TEST_BUGS_BIN): $(TEST_BUGS_SRC)
	$(FC) $(FFLAGS) $(TEST_BUGS_SRC) -o $(TEST_BUGS_BIN)

$(TEST_TLE_BIN): $(TEST_TLE_SRC)
	$(FC) $(FFLAGS) $(TEST_TLE_SRC) -o $(TEST_TLE_BIN)

$(TEST_TLE2SV_BIN): $(TEST_TLE2SV_SRC)
	$(FC) $(FFLAGS) $(TEST_TLE2SV_SRC) -o $(TEST_TLE2SV_BIN)

$(TEST_TLE2OPM_BIN): $(TEST_TLE2OPM_SRC)
	$(FC) $(FFLAGS) $(TEST_TLE2OPM_SRC) -o $(TEST_TLE2OPM_BIN)

run: $(TARGET)
	./$(TARGET)

test: $(TARGET) $(TLE2OPM_BIN) $(ALL_TEST_BINS)
	bash test_all.sh

clean:
	rm -f *.o *.obj *.mod *.smod *.ilk *.pdb *.optrpt
	rm -f $(TARGET) $(TLE2OPM_BIN) $(ALL_TEST_BINS)
	rm -f $(TARGET).exe $(TLE2OPM_BIN).exe
	rm -f $(addsuffix .exe,$(ALL_TEST_BINS))
