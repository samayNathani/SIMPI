#include <sys/mman.h>
#include <sys/stat.h>       
#include <fcntl.h>
#include <unistd.h>
#include <sys/time.h>
#include <vector>
#include <string>
#include <sstream>
class simpi
{
private:
	int id;
	std::vector<int> shared_mem_fds;
	std::vector<std::string> shared_mem_names;
	std::string  get_shared_mem_name() { 
		timeval curTime;
		gettimeofday(&curTime, NULL);
		unsigned long micro = curTime.tv_sec*(uint64_t)1000000+curTime.tv_usec;
		std::ostringstream ss;
		ss << "simpi_" << micro;
		return ss.str();
	}
	cleanup(std::string unique_id) {
		metadata = get_matrix_metadata(unique_id);
		// metadtaa.fd, metadtaa.name, metadtaa.mmap_addr
		closed(fd);
		munmap();
		shm_unlink();
	}

public:
	simpi(int id_)
	{
		id = id_;
	}
	~simpi()
	{
		//clean up shared mem
		for(int fd : shared_mem_fds) {
			close(fd);
		}
		for(std::string shared_mem_name : shared_mem_names) {
			shm_unlink(shared_mem_name.c_str());
		}

	}
	int get_id() { return id; }
	//
	double* create_matrix(int x, int y)
	{
    double size = x*y*sizeof(double);
		
		if (id == 0)
		{
		//add name to list of shared mems
		std::string shared_mem_name = get_shared_mem_name();
		shared_mem_names.push_back(shared_mem_name);

		//create a shared mem object
		int fd = shm_open(shared_mem_name.c_str(), O_RDWR | O_CREAT | O_EXCL, 0777);
		ftruncate(fd, size);
		double * matrix = (double*)mmap(NULL, sizeof(double)*size, PROT_READ|PROT_WRITE, MAP_SHARED, fd, 0);
		shared_mem_fds.push_back(fd);
    
		return matrix;
		}
		else
		{
			shared_mem_fds.back = shm_open(shared_mem_names.back.c_str(),O_RDWR, 0777);
			double * matrix = (double*)mmap(NULL, sizeof(double)*size, PROT_READ|PROT_WRITE, MAP_SHARED, shared_mem_fds.back, 0);
			return matrix;
		}
	}

	void free_matrix(std::string shared_mem_name) { 
		//get matrix metadata
		//close fd, shmunmap, munmap
	}

class matrix //similar stuff for vector
{
private:
	int xdim, ydim;
	std::string unique_id;
	simpi *psimpi = NULL; //for later reference
public:
	int getx() { return xdim; }
    int gety() {return ydim;}
	int getid() { return psimpi->get_id(); }
	double *arr;
	matrix(simpi &simp, int x, int y) //constructor
	{
		//use simp and init the matrix for all processes. The id is also in simp
		psimpi = &simp;
    arr = simp.create_matrix(x,y);
		xdim = x;
		ydim = y;
		// unique_id = id_from_simpi;
	}
	~matrix() //destructor
	{
		//use psimpi for getting rid of the mem and unlink stuff
		psimpi->free_matrix(unique_id);
	}
	double &get(int x, int y)
	{
		return arr[x + x*y];
	}
	void invert()
	{
	}
};

class vector //similar stuff for vector
{
private:
	int dim;
	simpi *psimpi = NULL; //for later reference
public:
	float *arr;
	vector(simpi *simp, int a) //constructor
	{
		//use simp and init the matrix for all processes. The id is also in simp
		psimpi = simp;
	}
	~matrix() //destructor
	{
		//use psimpi for getting rid of the mem and unlink stuff
	}
};



void main(int argc, char *argv[])
{

	int a = 5;
	fct(a);
	cout << a;

	int id = ...;
	simpi simp(id, ...);

	...

			matrix A(&simp, 10, 15);
	vector v(simp, 10);
	... A.get(4, 5) = 7;
	float f = A.get(1, 1);

	C = A * B;
}