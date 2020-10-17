PRGNAME=refflow_utils
CXX=g++
CXX_FLAGS=--std=c++14 -lpthread
LIB=-lhts

INC = split_sam.hpp merge_sam.hpp refflow_utils.hpp
OBJ = split_sam.o merge_sam.o refflow_utils.o

all: $(PRGNAME)

%.o: %.cpp $(INC)
	$(CXX) -c -o $@ $< $(LIB) $(CXX_FLAGS)

refflow_utils: $(OBJ)
	$(CXX) -o $@ $^ $(LIB) $(CXX_FLAGS)

clean:
	rm *.o $(PRGNAME)
