## Phaser
CC = g++-6
CFLAGS = -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_OMP -std=c++14 -static -Ofast -g3 -Wall -c -fopenmp -msse2 -mavx
LFLAGS = -static -static-libgcc -static-libstdc++ -fopenmp -g -o 



#CFLAGS = -std=c++14  -Ofast -g3 -Wall -c -fopenmp -msse2 -mavx
#LFLAGS = -fopenmp -g -o 


#CFLAGS = -std=c++14  -O3 -g3 -Wall -c -fopenmp 
#LFLAGS =  -fopenmp -g -o 


SOURCES=$(wildcard *.cpp)
OBJECTS=$(SOURCES:.cpp=.o)
IGNORES=HaplotypePhaserSym.o 


FILTERED=$(filter-out $(IGNORES),$(OBJECTS))

TARGET=phase_ompmod

all: $(TARGET)

# ../../../Programs/eigen-eigen-5a0156e40feb/


$(TARGET): $(FILTERED)
	$(CC) $(LFLAGS) $@ $^  -L../../haplotyperProject/libStatGen -lStatGen -lz 
#	$(CC) $(LFLAGS) $@ $^  -L../../haplotyperProject/libStatGen -lStatGen_debug -lz 
%.o: %.cpp %.h
	$(CC) $(CFLAGS) -g $< -I ../../haplotyperProject/libStatGen/include/ -I /home/kristiina/Projects/thrust/ -I ../../../Programs/eigen-3.3.9/ 
	
	
%.o: %.cpp
	 $(CC) $(CFLAGS) $< -I ../../haplotyperProject/libStatGen/include/ -I /home/kristiina/Projects/thrust/  -I ../../../Programs/eigen-3.3.9/   


clean:
	rm -f *.o 