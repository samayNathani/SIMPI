#include <sys/mman.h>
#include <sys/stat.h>       
#include <fcntl.h>
#include <unistd.h>
#include <sys/time.h>
#include <vector>
#include <string>
#include <unordered_map> 
#include <sstream>
#include <string>

#define SYNCH_OBJECT_MEM_NAME "simpi_shared_mem"

typedef struct matrix_metadata {
	std::string unique_id;
	int file_descriptor;
	size_t size;
	double* matrix_data;
} matrix_metadata;

typedef struct synch_object {
	int par_count;
	std::string last_matrix_id;
	int ready[];
} synch_object;


class simpi
{
public:
	simpi(int id_)
	{
		id = id_;
		int fd = shm_open(SYNCH_OBJECT_MEM_NAME, O_RDWR, 0777);
		synch_info = (synch_object*)mmap(NULL, sizeof(synch_object), PROT_READ|PROT_WRITE, MAP_SHARED,fd, 0);
	}
	~simpi()	
	{
		
		//clean up shared mem
		for(std::pair<std::string, matrix_metadata> matrix : matrix_info)
 		{
			free_matrix(matrix.second.unique_id);
		}
		shm_unlink(SYNCH_OBJECT_MEM_NAME);
	
	}
	int get_id() { return id; }

	std::pair<std::string, double*> create_matrix(int x, int y)
	{
    size_t size = x*y*sizeof(double);
		
		if (id == 0)
		{
		//generate a uniqueid for matrix
		std::string unique_id = get_shared_mem_name();
		//create a shared mem object
		int fd = shm_open(unique_id.c_str(), O_RDWR | O_CREAT | O_EXCL, 0777);
		ftruncate(fd, size);
		//allocate matrix
		double * matrix = (double*)mmap(NULL, sizeof(double)*size, PROT_READ|PROT_WRITE, MAP_SHARED, fd, 0);

		matrix_metadata metadata;
		metadata.size = size;
		metadata.file_descriptor = fd;
		metadata.unique_id = unique_id;
		metadata.matrix_data = matrix;
		matrix_info[unique_id] = metadata;

		//write name to synch_object so that other processes can get the uniqie id
		synch_info->last_matrix_id = unique_id;

		synch(id, synch_info->par_count, synch_info->ready);
		std::pair<std::string, double*> pass_back;
		pass_back = std::make_pair(unique_id, matrix);
		return pass_back;
		}
		else
		{
			//wait for id 0 to create the shared memory for matrix
			synch(id, synch_info->par_count, synch_info->ready);
			//get the unique id from the synch object
			std::string unique_id = synch_info->last_matrix_id;
			//open and allocate the shared memory
			int fd = shm_open(unique_id.c_str(), O_RDWR, 0777);
			double * matrix = (double*)mmap(NULL, sizeof(double)*size, PROT_READ|PROT_WRITE, MAP_SHARED, fd, 0);

			//create a metadata object
			matrix_metadata metadata;
			metadata.unique_id = unique_id;
			metadata.file_descriptor = fd;
			metadata.matrix_data = matrix;
			matrix_info[unique_id] = metadata;
			std::pair<std::string, double*> pass_back = std::make_pair(unique_id, matrix);
			return pass_back;
		}
	}

	void free_matrix(std::string unique_id) { 
		//get matrix metadata
		//close fd, shmunmap, munmap
		matrix_metadata metadata = matrix_info[unique_id];
		close(metadata.file_descriptor);
		shm_unlink(metadata.unique_id.c_str());
		munmap(metadata.matrix_data, metadata.size);
	}

private:
	int id;
	synch_object* synch_info;
	std::string sync_shared_mem_name; //this shared meory holds: shared mem name for last allocated matrix id, synch 
	std::unordered_map<std::string, matrix_metadata> matrix_info;
	std::string  get_shared_mem_name() { 
		//gets a unique name for each shared memory based on the time that each was made
		timeval curTime;
		gettimeofday(&curTime, NULL);
		unsigned long micro = curTime.tv_sec*(uint64_t)1000000+curTime.tv_usec;
		std::ostringstream ss;
		ss << "simpi_" << micro;
		return ss.str();
	}
	void  synch(int par_id, int par_count, int* ready)
	{
    int synchid = ready[par_count] + 1;  // synchid is progressing with each synch - I know, one from section 1 and 3 had a similar idea.
    ready[par_id] = synchid;
    int breakout = 0;
    while (1) {
        breakout = 1;
        for (int i = 0; i < par_count; i++) {
            if (ready[i] < synchid)  //"less than" is super important here. Do you know why? Because one process could run ahead and
                                     // increment its ready[par_is] var already.
            {
                breakout = 0;
                break;
            }
        }
        if (breakout == 1) {
            ready[par_count] = synchid;  // and here we increment the additional variable
            break;
        }
    }
 }
};
class matrix //similar stuff for vector
{
private:
	int xdim, ydim;
	std::string unique_id;
	simpi *psimpi = NULL; //for later reference
public:
	int get_x() { return xdim; }
  int get_y() {return ydim;}
	double *arr;
	matrix(simpi &simp, int x, int y) //constructor
	{
		//use simp and init the matrix for all processes. The id is also in simp
		psimpi = &simp;
    std::pair<std::string, double*> pass_back ( psimpi->create_matrix(x,y));
		unique_id = pass_back.first;
		arr = pass_back.second;
		xdim = x;
		ydim = y;
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
};

class vector //similar stuff for vector
{
private:
	int dim;
	simpi *psimpi = NULL; //for later reference
	std::string unique_id;
public:
	double *arr;

	vector(simpi &simp, int a) //constructor
	{
		//use simp and init the matrix for all processes. The id is also in simp
		psimpi = &simp;
    std::pair<std::string, double*> pass_back ( psimpi->create_matrix(1,a));
		unique_id = pass_back.first;
		arr = pass_back.second;
		dim = a;
	}
	~vector() //destructor
	{
		//use psimpi for getting rid of the mem and unlink stuff
		psimpi->free_matrix(unique_id);
	}
	double &get(int y)
	{
		return arr[y];
	}
};