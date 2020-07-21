TARGET: libgalpy.a petar.galpy

CXX=mpic++
CC=gcc

GALPY_PATH=$(HOME)/.local/lib/python3.6/site-packages/galpy
GALPY_INCLUDE_PATH= $(shell find $(GALPY_PATH) -type d )
GALPY_INCLUDE = $(addprefix -I, $(GALPY_INCLUDE_PATH))

CXXFLAGS +=  -O3 -Wall -std=c++11 
CFLAGS +=  -O3 -Wall $(GALPY_INCLUDE)

SRC:= $(shell find $(GALPY_PATH) -type f \( -iname '*.c' \) )
OBJ:= $(SRC:.c=.o)

libgalpy.a: $(OBJ)
	ar rcs $@ $^

petar.galpy: galpy_test.cxx libgalpy.a galpy_interface.h 
	$(CXX) $(CXXFLAGS) $(GALPY_INCLUDE) $< -o $@ -L./ -lgalpy -lgsl -lgslcblas

clean: 
	rm -f *.o petar.galpy libgalpy.a
