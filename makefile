## Phaser
CC = g++
CFLAGS = -std=c++14 -static -O3 -g3 -Wall -c -fopenmp 
LFLAGS = -static -static-libgcc -static-libstdc++ -fopenmp -g -o 
SOURCES=$(wildcard *.cpp)
OBJECTS=$(SOURCES:.cpp=.o)
TARGET=phase

all: $(TARGET)


$(TARGET): $(OBJECTS)
	$(CC) $(LFLAGS) $@ $^ ../../haplotyperProject/genoUtils/GenoUtils.o -L../../haplotyperProject/libStatGen -lStatGen -lz 
#	$(CC) $(LFLAGS) $@ $^ ../../haplotyperProject/genoUtils/GenoUtils.o -L../../haplotyperProject/libStatGen -lStatGen_debug -lz 
%.o: %.cpp %.h
	$(CC) $(CFLAGS) -g $< -I ../../haplotyperProject/libStatGen/include/ -I ../../haplotyperProject/genoUtils/ 
	
	
%.o: %.cpp
	 $(CC) $(CFLAGS) $< -I ../../haplotyperProject/libStatGen/include/ -I ../../haplotyperProject/genoUtils/ 


clean:
	rm -f *.o 