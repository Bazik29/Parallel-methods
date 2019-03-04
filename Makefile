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


default: $(TARGET)

$(TARGET): main.o integral.o utils.o
	$(COMPILER) $(FLAGS) main.o integral.o utils.o -o $@

main.o: main.cpp
	$(COMPILER) -c $(FLAGS) main.cpp

integral.o: integral.cpp
	$(COMPILER) -c $(FLAGS) integral.cpp

utils.o: utils.cpp
	$(COMPILER) -c $(FLAGS) utils.cpp


clean:
	$(RM) *.o

clean_all:
	$(RM) *.o main