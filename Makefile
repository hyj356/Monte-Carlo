FC = gfortran
TARGET = mc
FFLAGS = -O3 -Wall -fopenmp
OBJS = global_var.o mod_fileIO.o mod_computePE.o mod_random.o main.o
# FILES = ${wildcard *.f90}
# OBJS = ${patsubst %.f90, %.o, $(FILES) }

${TARGET}: ${OBJS}
	${FC} $^ -o $@ -fopenmp	

%.o: %.f90
	${FC} -c $< ${FFLAGS} -o $@

.PHONY: clean

clean:
	rm -rf *.mod *.o ${TARGET}