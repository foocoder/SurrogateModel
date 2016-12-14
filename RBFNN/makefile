OBJECTS = RadialBasisFunction.o testRBF.o
cc = clang++
override CFLAGS += -std=c++11

a.out: $(OBJECTS)
	$(cc) $(CFLAGS) -o a.out $(OBJECTS)

testRBF.o: RadialBasisFunction.h RadialBasisFunction.cpp testRBF.cpp
	$(cc) $(CFLAGS) -c testRBF.cpp

RadialBasisFunction.o: RadialBasisFunction.cpp matrix.h
	$(cc) $(CFLAGS) -c RadialBasisFunction.cpp


.PHONY: clean
clean:
	rm $(OBJECTS) a.out core
