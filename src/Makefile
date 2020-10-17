PRGNAME=refflow_utils
CXX=g++
CXX_FLAGS=--std=c++14 -lpthread
INC=
LIB=-lhts

DEPS = refflow_utils.hpp split_sam.hpp
OBJ = refflow_utils.o split_sam.o

all: $(PRGNAME)

%.o: %.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(LIB) $(CXX_FLAGS)

refflow_utils: $(OBJ)
	$(CXX) -o $@ $^ $(LIB) $(CXX_FLAGS)

clean:
	rm *.o $(PRGNAME)
