mpi : mpi.cpp user
	g++ mpi.cpp -o mpi -std=c++11
user : user.cpp
	g++ user.cpp -o user -std=c++11