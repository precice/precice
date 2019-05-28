F90 ?= gfortran

all: solverdummy

solverdummy: solverdummy.f90
	$(F90) $^ -o $@ $(shell pkg-config --libs libprecice)

clean:
	rm -f solverdummy
