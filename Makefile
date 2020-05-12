CC = g++ 
CFLAGS = -Wall -std=c++11

all : mpi user
clean : 
	rm -f user mpi /dev/shm/simpi_shared_mem
mpi : mpi.cpp user simpi.h
	$(CC) $(CFLAGS) mpi.cpp -o mpi 
user : user.cpp simpi.h
	$(CC) $(CFLAGS) user.cpp -o user 

#If a mac user, use the below makefile.
# CC = g++ 
# CFLAGS = -Wall -std=c++11

# all : mpi user
# clean : 
# 	rm -f user mpi /dev/shm/simpi_shared_mem
# mpi : mpi.cpp user simpi.h
# 	$(CC) $(CFLAGS) mpi.cpp -o mpi 
# user : user.cpp simpi.h
# 	$(CC) $(CFLAGS) user.cpp -o user 