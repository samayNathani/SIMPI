#include <fcntl.h>
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

class matrix  // similar stuff for vector
{
 private:
  int xdim, ydim;
  std::string unique_id;
  simpi* mysimpi = NULL;  // for later reference
 public:
  int get_x() { return xdim; }
  int get_y() { return ydim; }
  int determinant(double* A, int n, int order);
  void adjoint(double* A, double* adj, int order, int par_id, int par_count);
  void getCofactor(double* A, double* temp, int p, int q, int n, int order);
  double get_algbera(int pos) { return arr[pos]; }
  void set(int pos, int val) { arr[pos] = val; }
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
  matrix& inverse()
  {
    matrix* inverse = new matrix(*mysimpi, xdim, ydim);
    matrix* adj = new matrix(*mysimpi, xdim, ydim);

    // Find determinant of A[][]
    int det = determinant(arr, xdim, xdim);
    if (det == 0) {
      std::cout << "Singular matrix, can't find its inverse";
      return *inverse;
    }
    std::cout << "Determinant is: " << det;
    // std::chrono::steady_clock::time_point end1 =
    // std::chrono::steady_clock::now();

    // std::chrono::steady_clock::time_point begin2 =
    // std::chrono::steady_clock::now();

    // Find adjoint
    // float adj[order*order];
    mysimpi->synch();

    adjoint(arr, adj->arr, xdim, mysimpi->get_id(),
            mysimpi->get_synch_info()->par_count);
    mysimpi->synch();

    // std::chrono::steady_clock::time_point end2 =
    // std::chrono::steady_clock::now();

    // std::chrono::steady_clock::time_point begin3 =
    // std::chrono::steady_clock::now(); Find Inverse using formula "inverse(A)
    // = adj(A)/det(A)"
    for (int i = 0; i < xdim; i++)
      for (int j = 0; j < xdim; j++)
        inverse->get(i, j) = adj->get(i, j) / double(det);

    return *inverse;
  }
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

void simpi::synch()
{
  int par_count = synch_info->par_count;
  int* ready = synch_info->ready;
  int synchid = ready[par_count] + 1;
  ready[id] = synchid;
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
    int fd = shm_open(unique_id.c_str(), O_RDWR | O_CREAT, 0777);
    if (fd == -1) {
      std::string msg = std::to_string(id) + ": Unable to shm_open matrix (" +
                        unique_id + "):  ";
      perror(msg.c_str());
      exit(1);
    }
    ftruncate(fd, size);

    // allocate matrix
    double* matrix =
        (double*)mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    if (matrix == MAP_FAILED) {
      perror("Unable to mmap matrix: ");
      exit(1);
    }

    matrix_metadata metadata;
    metadata.size = size;
    metadata.file_descriptor = fd;
    strcpy(metadata.unique_id, unique_id.c_str());
    metadata.matrix_data = matrix;
    matrix_info[unique_id] = metadata;

    // write name to synch_object so that other processes can get the uniqie
    // id
    strcpy(synch_info->last_matrix_id, unique_id.c_str());
    synch();
    std::pair<std::string, double*> pass_back;
    pass_back = std::make_pair(unique_id, matrix);
    return pass_back;
  }
  else {
    // wait for id 0 to create the shared memory for matrix
    synch();
    // get the unique id from the synch object
    std::string unique_id = synch_info->last_matrix_id;
    // open and allocate the shared memory
    int fd = shm_open(unique_id.c_str(), O_RDWR, 0777);
    if (fd == -1) {
      std::string msg = std::to_string(id) + ": Unable to shm_open matrix (" +
                        unique_id + "):  ";
      perror(msg.c_str());
      exit(1);
    }
    double* matrix =
        (double*)mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    if (matrix == MAP_FAILED) {
      std::string msg = std::to_string(id) + ": Unable to mmap matrix: ";
      perror(msg.c_str());
      exit(1);
    }

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

// Functions for Inverse of Matrix
int matrix::determinant(double* A, int n, int order)
{
  int D = 0;  // Initialize result

  //  Base case : if matrix contains single element
  if (n == 1)
    return A[0];

  double temp[order * order];  // To store cofactors

  int sign = 1;  // To store sign multiplier

  // Iterate for each element of first row
  for (int f = 0; f < n; f++) {
    // Getting Cofactor of A[0][f]
    getCofactor(A, temp, 0, f, n, order);
    D += sign * A[0 + f * order] * determinant(temp, n - 1, order);

    // terms are to be added with alternate sign
    sign = -sign;
  }

  return D;
}

void matrix::adjoint(double* A,
                     double* adj,
                     int order,
                     int par_id,
                     int par_count)
{
  if (order == 1) {
    adj[0] = 1;
    return;
  }

  int rpp = order / par_count;
  int start = par_id * rpp;
  int end = start + rpp;

  // temp is used to store cofactors of A[][]
  int sign = 1;
  double temp[order * order];

  for (int i = 0; i < order; i++) {
    for (int j = start; j < end; j++) {
      // Get cofactor of A[i][j]
      getCofactor(A, temp, i, j, order, order);

      // sign of adj[j][i] positive if sum of row
      // and column indexes is even.
      sign = ((i + j) % 2 == 0) ? 1 : -1;

      // Interchanging rows and columns to get the
      // transpose of the cofactor matrix
      adj[j + i * order] = (sign) * (determinant(temp, order - 1, order));
    }
  }
}

void matrix::getCofactor(double* A,
                         double* temp,
                         int p,
                         int q,
                         int n,
                         int order)
{
  int i = 0, j = 0;

  // Looping for each element of the matrix
  for (int row = 0; row < n; row++) {
    for (int col = 0; col < n; col++) {
      //  Copying into temporary matrix only those element
      //  which are not in given row and column
      if (row != p && col != q) {
        temp[(i) + (j++) * order] = A[row + col * order];

        // Row is filled, so increase row index and
        // reset col index
        if (j == n - 1) {
          j = 0;
          i++;
        }
      }
    }
  }
}