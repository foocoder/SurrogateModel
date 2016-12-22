OBJECTS = mopm.o runmopm.o RadialBasisFunction.o
cc = g++
override CFLAGS += -std=c++11 -msse4.2



a.out: $(OBJECTS)
	$(cc) $(CFLAGS) -o a.out $(OBJECTS)

RadialBasisFunction.o: RadialBasisFunction.cpp matrix.h
	$(cc) $(CFLAGS) -c RadialBasisFunction.cpp
mopm.o: mopm.cpp mopm.h myBitSet.h
	$(cc) $(CFLAGS) -c mopm.cpp
runmopm.o: mopm.h runmopm.cpp myBitSet.h
	$(cc) $(CFLAGS) -c runmopm.cpp

clean:
	rm *.o a.out core
