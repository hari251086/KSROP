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

.PHONY: all tests run test clean

all: $(TARGET)

$(TARGET): $(SRC)
	$(FC) $(FFLAGS) $(SRC) -o $(TARGET)

tests: $(TEST_BIN)

$(TEST_BIN): $(TEST_SRC)
	$(FC) $(FFLAGS) $(TEST_SRC) -o $(TEST_BIN)

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
