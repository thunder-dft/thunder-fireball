BINS = fireball.x fireball-MDET.x

all: ${BINS}

.PHONY: clean

clean:
	rm -rf *.a *.o *.mod *genmod.f90 *~

%.o: %.f90
	${F90} ${FFLAGS} -c $< -I../include
