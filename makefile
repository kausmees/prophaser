############################ GPU version ############################

#CC=nvc++
#CFLAGS = -std=c++17 -mp=gpu -Minfo=mp -Wall -O4 -c -fast  -Mipa=fast -v -stdpar   -acc --remarks  -Mlre=assoc -Mdse -DTARGET_GPU
#LFLAGS = -gpu=cc80,managed -mp=gpu -acc -o 
#TARGET=phase_ompmod_gpu

############################ CPU version############################

CC = g++


CFLAGS = -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_OMP -std=c++14 -static -Ofast -g3 -Wall -c -fopenmp -msse2 -mavx
LFLAGS =  -static -static-libgcc -static-libstdc++ -fopenmp -g -o 

TARGET=phase_ompmod_cpu
####################################################################





SOURCES=$(wildcard *.cpp)
OBJECTS=$(SOURCES:.cpp=.o)
IGNORES=HaplotypePhaserSym.o testmain.o
FILTERED=$(filter-out $(IGNORES),$(OBJECTS))
all: $(TARGET)


############################ GPU version ############################
#
#$(TARGET): $(FILTERED)
#	$(CC) $(LFLAGS) $@ $^ -L../../haplotyperProject/libStatGen -lStatGen -lz -cuda 
#
#HaplotypPhaser.o: HaplotypePhaser.cpp HaplotypePhaser.h
#	$(CC) $(CFLAGS)  -gpu=cc80,managed -cuda $< -I ../../haplotyperProject/libStatGen/include/ -I ./../../Programs/eigen-3.3.9/  -I ../libdivide -cuda
#
#
#%.o: %.cpp %.h
#	$(CC) $(CFLAGS)  $< -I ../../haplotyperProject/libStatGen/include/ -I ./../../Programs/eigen-3.3.9/  -I ../libdivide 
#
#
#%.o: %.cpp
#	 $(CC) $(CFLAGS) $< -I ../../haplotyperProject/libStatGen/include/  -I ./../../Programs/eigen-3.3.9/  -I ../libdivide 


############################ CPU version ############################


$(TARGET): $(FILTERED)
	$(CC) $(LFLAGS) $@ $^  -L../../haplotyperProject/libStatGen -lStatGen -lz
#       $(CC) $(LFLAGS) $@ $^  -L../../haplotyperProject/libStatGen -lStatGen_debug -lz

%.o: %.cpp %.h
	$(CC) $(CFLAGS) -g $< -I ../../haplotyperProject/libStatGen/include/ -I /home/kristiina/Projects/thrust/ -I ../../../Programs/eigen-3.3.9/


%.o: %.cpp
	$(CC) $(CFLAGS) $< -I ../../haplotyperProject/libStatGen/include/   -I ../../../Programs/eigen-3.3.9/


clean:
	rm -f *.o
