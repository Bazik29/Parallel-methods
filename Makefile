TARGET = main

COMPILER = g++
FLAGS = -std=c++0x -O2 -Wall -Wextra -Wpedantic -fopenmp

OS := $(shell uname)
ifeq ($(OS),Darwin)
	COMPILER = g++-7
endif

ifdef SYSTEMROOT
    RM = del /Q
else
    RM = rm -f
endif


default: main1 test1 clean

main1: main1.o integral.o utils.o
	$(COMPILER) $(FLAGS) main1.o integral.o utils.o -o $@

main1.o: main1.cpp
	$(COMPILER) -c $(FLAGS) main1.cpp

integral.o: integral.cpp
	$(COMPILER) -c $(FLAGS) integral.cpp

utils.o: utils.cpp
	$(COMPILER) -c $(FLAGS) utils.cpp

test1: test1.o
	$(COMPILER) $(FLAGS) test1.o integral.o utils.o -o $@

test1.o: test1.cpp
	$(COMPILER) -c $(FLAGS) test1.cpp

clean:
	$(RM) *.o

clean_all:
	$(RM) *.o main1 test1
