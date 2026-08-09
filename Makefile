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

SRC     = app/driver_KS.F src/Subrouts.F src/Legendre.F src/TLEread.F \
          src/swx.F src/Cunningham.F src/LegendreTess.F
TARGET  = driver_KS

TLE2OPM_SRC = app/tle2opm.F src/Subrouts.F src/TLEread.F src/Legendre.F
TLE2OPM_BIN = tle2opm

TEST_SRC  = test/test_subrouts.F src/Subrouts.F src/Legendre.F
TEST_BIN  = test_subrouts

TEST_BUGS_SRC = test/test_bugs.F src/Subrouts.F src/Legendre.F
TEST_BUGS_BIN = test_bugs

TEST_TLE_SRC  = test/test_tle.F src/Subrouts.F src/TLEread.F src/Legendre.F
TEST_TLE_BIN  = test_tle

TEST_TLE2SV_SRC = test/test_tle2sv.F src/Subrouts.F src/TLEread.F src/Legendre.F
TEST_TLE2SV_BIN = test_tle2sv

TEST_TLE2OPM_SRC = test/test_tle2opm.F src/Subrouts.F src/TLEread.F src/Legendre.F
TEST_TLE2OPM_BIN = test_tle2opm

TEST_SW_SRC = test/test_sw.F src/swx.F src/Subrouts.F src/Legendre.F
TEST_SW_BIN = test_sw

TEST_CUNNINGHAM_SRC = test/test_cunningham.F src/Cunningham.F \
                      src/Subrouts.F src/Legendre.F
TEST_CUNNINGHAM_BIN = test_cunningham

TEST_LEGENDRE_TESS_SRC = test/test_legendre_tess.F src/LegendreTess.F \
                          src/Cunningham.F src/Subrouts.F src/Legendre.F
TEST_LEGENDRE_TESS_BIN = test_legendre_tess

TEST_DVDT_TESS_SRC = test/test_dvdt_tess.F src/LegendreTess.F \
                      src/Cunningham.F src/Subrouts.F src/Legendre.F
TEST_DVDT_TESS_BIN = test_dvdt_tess

ALL_TEST_BINS = $(TEST_BIN) $(TEST_BUGS_BIN) $(TEST_TLE_BIN) \
                $(TEST_TLE2SV_BIN) $(TEST_TLE2OPM_BIN) $(TEST_SW_BIN) \
                $(TEST_CUNNINGHAM_BIN) $(TEST_LEGENDRE_TESS_BIN) \
                $(TEST_DVDT_TESS_BIN)

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

$(TEST_SW_BIN): $(TEST_SW_SRC)
	$(FC) $(FFLAGS) $(TEST_SW_SRC) -o $(TEST_SW_BIN)

$(TEST_CUNNINGHAM_BIN): $(TEST_CUNNINGHAM_SRC)
	$(FC) $(FFLAGS) $(TEST_CUNNINGHAM_SRC) -o $(TEST_CUNNINGHAM_BIN)

$(TEST_LEGENDRE_TESS_BIN): $(TEST_LEGENDRE_TESS_SRC)
	$(FC) $(FFLAGS) $(TEST_LEGENDRE_TESS_SRC) -o $(TEST_LEGENDRE_TESS_BIN)

$(TEST_DVDT_TESS_BIN): $(TEST_DVDT_TESS_SRC)
	$(FC) $(FFLAGS) $(TEST_DVDT_TESS_SRC) -o $(TEST_DVDT_TESS_BIN)

run: $(TARGET)
	./$(TARGET)

test: $(TARGET) $(TLE2OPM_BIN) $(ALL_TEST_BINS)
	bash test_all.sh

clean:
	rm -f *.o *.obj *.mod *.smod *.ilk *.pdb *.optrpt
	rm -f $(TARGET) $(TLE2OPM_BIN) $(ALL_TEST_BINS)
	rm -f $(TARGET).exe $(TLE2OPM_BIN).exe
	rm -f $(addsuffix .exe,$(ALL_TEST_BINS))
