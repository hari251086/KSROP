# -------------------------------------------------------
# KSROP Makefile
# Usage:
#   make            build propagator
#   make tests      build unit test binary
#   make run        run propagator
#   make test       run all tests (unit + integration)
#   make clean      remove build artefacts
# -------------------------------------------------------

FC      = ifort
FFLAGS  = -O2 -warn all

SRC     = driver_KS.F Subrouts.F Legendre.F
TARGET  = driver_KS

TEST_SRC  = test_subrouts.F Subrouts.F Legendre.F
TEST_BIN  = test_subrouts

TEST_TLE_SRC  = test_tle.F TLEread.F
TEST_TLE_BIN  = test_tle

TEST_TLE2SV_SRC = test_tle2sv.F TLEread.F
TEST_TLE2SV_BIN = test_tle2sv

.PHONY: all tests run test clean

all: $(TARGET)

$(TARGET): $(SRC)
	$(FC) $(FFLAGS) $(SRC) -o $(TARGET)

tests: $(TEST_BIN) $(TEST_TLE_BIN) $(TEST_TLE2SV_BIN)

$(TEST_BIN): $(TEST_SRC)
	$(FC) $(FFLAGS) $(TEST_SRC) -o $(TEST_BIN)

$(TEST_TLE_BIN): $(TEST_TLE_SRC)
	$(FC) $(FFLAGS) $(TEST_TLE_SRC) -o $(TEST_TLE_BIN)

$(TEST_TLE2SV_BIN): $(TEST_TLE2SV_SRC)
	$(FC) $(FFLAGS) $(TEST_TLE2SV_SRC) -o $(TEST_TLE2SV_BIN)

run: $(TARGET)
	./$(TARGET)

test: $(TARGET) $(TEST_BIN)
	@echo "=== Unit tests ==="
	./$(TEST_BIN)
	@echo ""
	@echo "=== Integration test ==="
	python test_driver.py ./$(TARGET)

clean:
	rm -f *.o *.obj *.mod *.smod *.ilk *.pdb *.optrpt
	rm -f $(TARGET) $(TEST_BIN)
	rm -f $(TARGET).exe $(TEST_BIN).exe
