CC = g++ 
CFLAGS = -Wall -gstabs -std=c++11

all : mpi user
clean : 
	rm -f user mpi /dev/shm/simpi_shared_mem
mpi : mpi.cpp user simpi_temp.h
	$(CC) $(CFLAGS) mpi.cpp -o mpi -lrt
user : user.cpp simpi_temp.h
	$(CC) $(CFLAGS) user.cpp -o user -lrt