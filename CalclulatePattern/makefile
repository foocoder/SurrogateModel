obj = CalclulatePattern.o
cc  = g++
override CFLAGS += -std=c++11 -msse4.2

a.out : $(obj)
	$(cc) $(CFLAGS) -o a.out $(obj)

CalclulatePattern.o : myBitSet.h CalclulatePattern.h CalclulatePattern.cpp
	$(cc) $(CFLAGS) -c CalclulatePattern.cpp

.PHONY: clean
clean:
	rm a.out $(obj)
