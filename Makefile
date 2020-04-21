mpi.o : mpi.cpp user.o
	g++ mpi.cpp -o mpi.o -std=c++11 -lrt
user.o : user.cpp
	g++ user.cpp -o user.o -std=c++11 -lrt