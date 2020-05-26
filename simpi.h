#include <fcntl.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <unistd.h>

#include <iostream>
#include <map>
#include <string>
#include <vector>

#define SYNCH_OBJECT_MEM_NAME "/simpi_shared_mem"
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

// static methods
void INIT_SIMPI(int par_id, size_t synch_size);

class simpi {
 public:
  simpi(int id_, size_t synch_size)
  {
    id = id_;
    int fd = shm_open(SYNCH_OBJECT_MEM_NAME, O_RDWR, 0777);
    if (fd == -1) {
      perror("Unable to shm_open synch_object: ");
      exit(1);
    }
    synch_info = (synch_object*)mmap(NULL, synch_size, PROT_READ | PROT_WRITE,
                                     MAP_SHARED, fd, 0);
    if (synch_info == MAP_FAILED) {
      perror("Unable to mmap synch_info: ");
      exit(1);
    }
  }
  ~simpi()
  {
    for (std::pair<std::string, matrix_metadata> matrix : matrix_info) {
      free_matrix(matrix.second.unique_id);
    }
    shm_unlink(SYNCH_OBJECT_MEM_NAME);
  }
  int get_id() { return id; }
  synch_object* get_synch_info() { return synch_info; }

  std::pair<std::string, double*> create_matrix(int x, int y);

  void free_matrix(std::string unique_id);
  void synch();

 private:
  int id;
  synch_object* synch_info;
  std::map<std::string, matrix_metadata> matrix_info;
  std::string sync_shared_mem_name;
  std::string get_shared_mem_name();
};

class vector  // similar stuff for vector
{
 private:
  int dim;
  simpi* mysimpi = NULL;  // for later reference
  std::string unique_id;

 public:
  double* arr;
  int get_size() { return dim; }
  void set(int pos, double val) { arr[pos] = val; }
  double& get(int pos) { return arr[pos]; }
  void print();
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
};

class matrix  // similar stuff for vector
{
 private:
  int xdim, ydim;
  std::string unique_id;
  simpi* mysimpi = NULL;  // for later reference
 public:
  double* arr;
  void print();
  int get_x() { return xdim; }
  int get_y() { return ydim; }
  int determinant(double* A, int n, int order);
  void adjoint(double* A, double* adj, int order, int par_id, int par_count);
  void getCofactor(double* A, double* temp, int p, int q, int n, int order);
  double get_algbera(int pos) { return arr[pos]; }
  void set(int pos, int val) { arr[pos] = val; }
  matrix& inverse();
  void solveSystem(vector* constants, vector* solution);
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
