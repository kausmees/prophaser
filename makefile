############################## nvc++ ################################

# module load buildenv-nvhpc/21.5-cu11a 


#CC = hpcsdk213/Linux_x86_64/2021/compilers/bin/nvc++
CC=nvc++
CFLAGS = -std=c++17 -mp=gpu -Minfo=mp -Wall -O4 -c -fast  -Mipa=fast -v -stdpar   -acc --remarks  -Mlre=assoc -Mdse -DTARGET_GPU
#--remarks
LFLAGS = -gpu=cc80,managed -mp=gpu -acc -o 
TARGET=phase_ompmod_berzelius

#CFLAGS = -std=c++17 -Wall -c -v  -fast -O4 --remarks -Mlre=assoc -Mipa=fast -Mdse
#TARGET=phase_ompmod_berzelius_cpu


# -mp=gpu
#  -tp=sandybridge
#  -fno-exceptions
# -Mipa -Mkeepasm
#-gopt -Mdwarf3 -Memit-dwarf-inlined
# -Mneginfo=all -Minfo=all

#LFLAGS =  -g -mp=gpu -gpu=cc80,managed -cuda -o

############################## gcc++ ################################

#module load buildenv-gcccuda/11.2-8.3.1-bare

#CC = g++


#CFLAGS = -std=c++17  -Ofast -g3 -Wall -c -fopenmp -msse2 -mavx -fsanitize=address
#LFLAGS = -fopenmp -g -lasan -L /gcc/x86_64-redhat-linux/8/ -o 

#-fsanitize=address

#CFLAGS = -std=c++14  -O3 -g3 -Wall -c -fopenmp 
#LFLAGS =  -fopenmp -g -o 
#TARGET=phase_ompmod_berzelius_cpu_gccasan
###########################################################

############################## gcc++ ################################


#module load gcc/system

#CC = g++


#CFLAGS = -std=c++14 -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_OMP  -Ofast -g3 -Wall -c -fopenmp -msse2 -mavx
#-fsanitize=address
#LFLAGS = -fopenmp -g -o

#-fsanitize=address

#CFLAGS = -std=c++14  -O3 -g3 -Wall -c -fopenmp 
#LFLAGS =  -fopenmp -g -o 
#TARGET=phase_ompmod_berzelius_cpu_gcc
###########################################################





SOURCES=$(wildcard *.cpp)
OBJECTS=$(SOURCES:.cpp=.o)

IGNORES=HaplotypePhaserSym.o testmain.o

FILTERED=$(filter-out $(IGNORES),$(OBJECTS))





all: $(TARGET)


$(TARGET): $(FILTERED)
	$(CC) $(LFLAGS) $@ $^ -L/proj/gcae_berzelius/users/carl/libStatGen -lStatGen -lz -cuda 
#	$(CC) $(LFLAGS) $@ $^  -L../../haplotyperProject/libStatGen -lStatGen_debug -lz 
HaplotypPhaser.o: HaplotypePhaser.cpp HaplotypePhaser.h
	$(CC) $(CFLAGS)  -gpu=cc80,managed -cuda $< -I /proj/gcae_berzelius/users/carl/libStatGen/include/ -I /home/x_kriau/prophaser_libs/eigen-3.3.9/  -I ../libdivide -cuda
#-I /home/x_kriau/prophaser_libs/thrust/ -I /home/x_kriau/prophaser_libs/thrust/dependencies/cub/ -cuda


%.o: %.cpp %.h
	$(CC) $(CFLAGS)  $< -I /proj/gcae_berzelius/users/carl/libStatGen/include/ -I /home/x_kriau/prophaser_libs/eigen-3.3.9/  -I ../libdivide 
#-I /home/x_kriau/prophaser_libs/thrust/ -I /home/x_kriau/prophaser_libs/thrust/dependencies/cub/


%.o: %.cpp
	 $(CC) $(CFLAGS) $< -I /home/x_kriau/prophaser_libs/libStatGen/include/  -I /home/x_kriau/prophaser_libs/eigen-3.3.9/  -I ../libdivide 
#-I /home/x_kriau/prophaser_libs/thrust/ -I /home/x_kriau/prophaser_libs/thrust/dependencies/cub/


clean:
	rm -f *.o 
