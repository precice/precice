F90 ?= gfortran

all: solverdummy solverdummy_A solverdummy_B

solverdummy: solverdummy.f90
	$(F90) $^ -o $@ $(shell pkg-config --libs libprecice)

solverdummy_A: solverdummy_A.f90
	$(F90) $^ -o $@ $(shell pkg-config --libs libprecice)

solverdummy_B: solverdummy_B.f90
	$(F90) $^ -o $@ $(shell pkg-config --libs libprecice)

clean:
	rm -f solverdummy solverdummy_A solverdummy_B
