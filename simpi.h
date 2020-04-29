#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <unistd.h>

#include <iostream>
#include <map>
#include <string>
#include <vector>

#define SYNCH_OBJECT_MEM_NAME "simpi_shared_mem"
#define UNIQUE_ID_SIZE 23
typedef struct matrix_metadata {
  char unique_id[UNIQUE_ID_SIZE];
  int file_descriptor;
  size_t size;
  double* matrix_data;
} matrix_metadata;

typedef struct synch_object {
  int par_count;
  char last_matrix_id[UNIQUE_ID_SIZE];
  int ready[];
} synch_object;

class simpi {
 public:
  simpi(int id_, size_t synch_size)
  {
    id = id_;
    int fd = shm_open(SYNCH_OBJECT_MEM_NAME, O_RDWR, 0777);
    synch_info = (synch_object*)mmap(NULL, synch_size, PROT_READ | PROT_WRITE,
                                     MAP_SHARED, fd, 0);
  }
  ~simpi()
  {
    for (std::pair<std::string, matrix_metadata> matrix : matrix_info) {
      free_matrix(matrix.second.unique_id);
    }
    shm_unlink(SYNCH_OBJECT_MEM_NAME);
  }
  int get_id() { return id; }

  std::pair<std::string, double*> create_matrix(int x, int y);

  void free_matrix(std::string unique_id);

 private:
  int id;
  synch_object* synch_info;
  std::map<std::string, matrix_metadata> matrix_info;
  std::string sync_shared_mem_name;
  std::string get_shared_mem_name();
  void synch(int par_id, int par_count, int* ready);
};
class matrix  // similar stuff for vector
{
 private:
  int xdim, ydim;
  std::string unique_id;
  simpi* mysimpi = NULL;  // for later reference
 public:
  int get_x() { return xdim; }
  int get_y() { return ydim; }
  double* arr;
  matrix(simpi& simp, int x, int y)  // constructor
  {
    // use simp and init the matrix for all processes. The id is also in simp
    mysimpi = &simp;
    std::pair<std::string, double*> pass_back(mysimpi->create_matrix(x, y));
    unique_id = pass_back.first;
    arr = pass_back.second;
    xdim = x;
    ydim = y;
  }
  ~matrix()  // destructor
  {
    // use mysimpi for getting rid of the mem and unlink stuff
    mysimpi->free_matrix(unique_id);
  }
  double& get(int x, int y) { return arr[x + y * xdim]; }
};

class vector  // similar stuff for vector
{
 private:
  int dim;
  simpi* mysimpi = NULL;  // for later reference
  std::string unique_id;

 public:
  double* arr;

  vector(simpi& simp, int a)  // constructor
  {
    // use simp and init the matrix for all processes. The id is also in simp
    mysimpi = &simp;
    std::pair<std::string, double*> pass_back(mysimpi->create_matrix(1, a));
    unique_id = pass_back.first;
    arr = pass_back.second;
    dim = a;
  }
  ~vector()  // destructor
  {
    // use mysimpi for getting rid of the mem and unlink stuff
    mysimpi->free_matrix(unique_id);
  }
  double& get(int y) { return arr[y]; }
};

std::string simpi::get_shared_mem_name()
{
  // gets a unique name for each shared memory based on the time that each was
  // made
  timeval curTime;
  gettimeofday(&curTime, NULL);
  unsigned long micro = curTime.tv_sec * (uint64_t)1000000 + curTime.tv_usec;
  std::string name = "simpi_" + std::to_string(micro);
  return name;
}

void simpi::synch(int par_id, int par_count, int* ready)
{
  int synchid = ready[par_count] + 1;
  ready[par_id] = synchid;
  int breakout = 0;
  while (1) {
    breakout = 1;
    for (int i = 0; i < par_count; i++) {
      if (ready[i] < synchid) {
        breakout = 0;
        break;
      }
    }
    if (breakout == 1) {
      ready[par_count] = synchid;
      // and here we increment the additional variable
      break;
    }
  }
}

std::pair<std::string, double*> simpi::create_matrix(int x, int y)
{
  size_t size = x * y * sizeof(double);

  if (id == 0) {
    // generate a uniqueid for matrix
    std::string unique_id = get_shared_mem_name();
    // create a shared mem object
    int fd = shm_open(unique_id.c_str(), O_RDWR | O_CREAT | O_EXCL, 0777);
    ftruncate(fd, size);
    // allocate matrix
    double* matrix = (double*)mmap(NULL, sizeof(double) * size,
                                   PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);

    matrix_metadata metadata;
    metadata.size = size;
    metadata.file_descriptor = fd;
    strcpy(metadata.unique_id, unique_id.c_str());
    metadata.matrix_data = matrix;
    matrix_info[unique_id] = metadata;

    // write name to synch_object so that other processes can get the uniqie
    // id
    strcpy(synch_info->last_matrix_id, unique_id.c_str());
    synch(id, synch_info->par_count, synch_info->ready);
    std::pair<std::string, double*> pass_back;
    pass_back = std::make_pair(unique_id, matrix);
    return pass_back;
  }
  else {
    // wait for id 0 to create the shared memory for matrix
    synch(id, synch_info->par_count, synch_info->ready);
    // get the unique id from the synch object
    std::string unique_id = synch_info->last_matrix_id;
    // open and allocate the shared memory
    int fd = shm_open(unique_id.c_str(), O_RDWR, 0777);
    double* matrix = (double*)mmap(NULL, sizeof(double) * size,
                                   PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);

    // create a metadata object
    matrix_metadata metadata;
    strcpy(metadata.unique_id, unique_id.c_str());
    metadata.file_descriptor = fd;
    metadata.matrix_data = matrix;
    matrix_info[unique_id] = metadata;
    std::pair<std::string, double*> pass_back =
        std::make_pair(unique_id, matrix);
    return pass_back;
  }
}

void simpi::free_matrix(std::string unique_id)
{
  // get matrix metadata
  // close fd, shmunmap, munmap
  matrix_metadata metadata = matrix_info[unique_id];
  close(metadata.file_descriptor);
  shm_unlink(metadata.unique_id);
  munmap(metadata.matrix_data, metadata.size);
}
