PRGNAME=refflow_utils
CXX=g++
CXX_FLAGS=--std=c++14 -lpthread
INC=
LIB=-lhts

all: $(PRGNAME)

refflow_utils: refflow_utils.cpp refflow_utils.hpp
	$(CXX) -o $@ $< $(LIB) $(CXX_FLAGS)

clean:
	rm *.o $(PRGNAME)
