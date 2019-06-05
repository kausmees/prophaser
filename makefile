## Phaser
CC = g++
CFLAGS = -std=c++14 -static -Ofast -g3 -Wall -c -fopenmp -msse2 -mavx
LFLAGS = -static -static-libgcc -static-libstdc++ -fopenmp -g -o 

#CFLAGS = -std=c++14  -Ofast -g3 -Wall -c -fopenmp -msse2 -mavx
#LFLAGS = -fopenmp -g -o 


#CFLAGS = -std=c++14  -O3 -g3 -Wall -c -fopenmp 
#LFLAGS =  -fopenmp -g -o 


SOURCES=$(wildcard *.cpp)
OBJECTS=$(SOURCES:.cpp=.o)
TARGET=phase

all: $(TARGET)


$(TARGET): $(OBJECTS)
	$(CC) $(LFLAGS) $@ $^  -L../../haplotyperProject/libStatGen -lStatGen -lz 
#	$(CC) $(LFLAGS) $@ $^  -L../../haplotyperProject/libStatGen -lStatGen_debug -lz 
%.o: %.cpp %.h
	$(CC) $(CFLAGS) -g $< -I ../../haplotyperProject/libStatGen/include/ -I ../../../Programs/eigen-eigen-5a0156e40feb/
	
	
%.o: %.cpp
	 $(CC) $(CFLAGS) $< -I ../../haplotyperProject/libStatGen/include/  -I ../../../Programs/eigen-eigen-5a0156e40feb/


clean:
	rm -f *.o 