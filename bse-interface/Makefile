TARGET: libbse.a

BSE_PATH=bse
FC=gfortran
FFLAGS= -g -O0 -Wall

VPATH=$(BSE_PATH)

SRC:= $(shell find $(BSE_PATH) -type f  \( -iname '*.f' ! -iname '.*' ! -iname '[sb]se.f' ! -iname 'popbin.f' \) )
OBJ:= $(SRC:.f=.o)

libbse.a: $(OBJ)
	ar rcs $@ $^

bse_test: libbse.a bse_test.cxx bse_interface.h
	$(CXX) -g -O0 -Wall -std=c++11 -fopenmp bse_test.cxx -o $@ -L./ -lbse -lgfortran

clean:
	rm -f $(BSE_PATH)/*.o *.o bse_test libbse.a
