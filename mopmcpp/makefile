OBJECTS = mopm.o runmopm.o
cc = g++
override CFLAGS += -std=c++11 -msse4.2



runmopm: $(OBJECTS)
	$(cc) $(CFLAGS) -o runmopm $(OBJECTS)

mopm.o: mopm.cpp mopm.h myBitSet.h
	$(cc) $(CFLAGS) -c mopm.cpp
runmopm.o: mopm.h runmopm.cpp myBitSet.h
	$(cc) $(CFLAGS) -c runmopm.cpp

clean:
	rm *.o runmopm
