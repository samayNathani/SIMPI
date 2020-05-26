#include "simpi.h"
// init global simpi
static simpi* main_simpi;

// simpi init function

// simpi shutdown?

/******************Simpi Functions*************************/
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
/******************Matrix Functions*************************/

matrix& matrix::inverse()
{
  matrix* inverse = new matrix(*main_simpi, xdim, ydim);
  matrix* adj = new matrix(*main_simpi, xdim, ydim);

  // Find determinant of A[][]
  int det = determinant(arr, xdim, xdim);
  if (det == 0) {
    std::cout << "Singular matrix, can't find its inverse";
    return *inverse;
  }
  std::cout << "Determinant is: " << det;
  main_simpi->synch();

  adjoint(arr, adj->arr, xdim, main_simpi->get_id(),
          main_simpi->get_synch_info()->par_count);
  main_simpi->synch();

  for (int i = 0; i < xdim; i++)
    for (int j = 0; j < xdim; j++)
      inverse->get(i, j) = adj->get(i, j) / double(det);

  return *inverse;
}

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

void matrix::solveSystem(vector* constants, vector* solution)
{
  int processCount = main_simpi->get_synch_info()->par_count;
  int id = main_simpi->get_id();

  // TODO: initialize shared memory if id is 0
  vector* prev = new vector(
      *main_simpi,
      constants->get_size());  // shared mem containing a copy of values
  // vector *solution = new vector(*main_simpi, constants->get_size()); //
  // shared mem containing actual calculated values

  matrix* saveEq = new matrix(*main_simpi, get_x(),
                              get_y());  // save equations from modification
  vector* saveConst =
      new vector(*main_simpi, constants->get_size());  // saves input vector

  // divide up work
  int n = constants->get_size();

  int work = n / processCount;
  /*
   * implement remainder?
   */
  int i, j, k;
  int start = id * work;
  int end = start + work;

  // Save Matrix and Vector
  for (i = start; i < end; i++) {
    for (j = 0; j < get_y(); j++) {
      saveEq->get(i, j) = get(i, j);
    }
    saveConst->get(i) = constants->get(i);
  }

  // synch, wait for all process before solving
  main_simpi->synch();

  // setup, switch var coefficient with row solution, and divide by coefficient
  for (i = start; i < end; i++) {
    double temp = get(i, i);
    get(i, i) = constants->get(i);
    constants->get(i) = temp;
    for (j = 0; j < get_y(); j++) {
      if (j != i) {
        get(i, j) *= -1;
      }
      get(i, j) /= constants->get(i);
      main_simpi->synch();
    }
    prev->get(i) = 0.0;
    solution->get(i) = constants->get(i);
    main_simpi->synch();
  }

  main_simpi->synch();
  // first iteration by trying substituting 1
  for (i = start; i < end; i++) {
    double rowSum = 0;
    for (j = 0; j < get_y(); j++) {
      if (j == i) {
        rowSum += get(i, j);
      }
      else {
        rowSum += (get(i, j) * prev->get(j));
      }
      main_simpi->synch();
    }
    solution->get(i) = rowSum;
    main_simpi->synch();
  }

  // wait for all processes before repeating iterations with calculated results
  main_simpi->synch();

  for (k = 0; k < 1000; k++) {
    for (i = start; i < end; i++) {
      // save prev value for comparision
      prev->get(i) = solution->get(i);
      main_simpi->synch();
    }
    for (i = start; i < end; i++) {
      // save prev value for comparision
      // prev->get(i) = solution->get(i);
      double rowSum = 0;
      for (j = 0; j < get_y(); j++) {
        if (j == i) {
          rowSum += get(i, j);
        }
        else {
          rowSum += (get(i, j) * prev->get(j));
        }
      }
      solution->get(i) = rowSum;
      main_simpi->synch();
    }
    // wait at end of each loop for all processes before beginning next
    // iteration
    main_simpi->synch();
  }
  main_simpi->synch();
  // restore original matrix and vector
  for (i = start; i < end; i++) {
    for (j = 0; j < get_y(); j++) {
      get(i, j) = saveEq->get(i, j);
    }
    constants->get(i) = saveConst->get(i);
  }
  // wait for all processes before returning solution vector
  main_simpi->synch();
  // return solution;
}