## Phaser
CC = gcc
#CFLAGS = -std=c++14 -static -Ofast -g3 -Wall -c -fopenmp -msse2 -mavx
CFLAGS = -std=c++14 -static -Ofast -g3 -Wall -c -fopenmp -DNDEBUG 
# -static -static-libgcc -static-libstdc++ 
LFLAGS = -fopenmp -g

#CFLAGS = -std=c++14  -Ofast -g3 -Wall -c -fopenmp -msse2 -mavx
#LFLAGS = -fopenmp -g -o 


#CFLAGS = -std=c++14  -O3 -g3 -Wall -c -fopenmp 
#LFLAGS =  -fopenmp -g -o 


SOURCES=$(wildcard *.cpp)
OBJECTS=$(SOURCES:.cpp=.o)
TARGET=phase
LIBSTATGEN=../libStatGen

all: $(TARGET)


$(TARGET): $(OBJECTS)
	$(CC) $(LFLAGS) $@ $^  -L$(LIBSTATGEN) -lStatGen -lz 
#	$(CC) $(LFLAGS) $@ $^  -L$(LIBSTATGEN) -lStatGen_debug -lz 
%.o: %.cpp %.h
	$(CC) $(CFLAGS) -g $< -I $(LIBSTATGEN)/include/ -I /usr/include/eigen3 #-I ../../../Programs/eigen-eigen-5a0156e40feb/
	
	
%.o: %.cpp
	 $(CC) $(CFLAGS) $< -I $(LIBSTATGEN)/include/ -I /usr/include/eigen3 #-I ../../../Programs/eigen-eigen-5a0156e40feb/


clean:
	rm -f *.o 
