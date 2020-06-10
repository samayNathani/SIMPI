#include "simpi.h"
// init global simpi
static simpi* main_simpi;

// simpi init function
void SIMPI_INIT(int par_id, size_t synch_size)
{
  main_simpi = new simpi(par_id, synch_size);
}
// simpi synch
void SIMPI_SYNCH()
{
  main_simpi->synch();
}

void SIMPI_FINALIZE()
{
  main_simpi->synch();
  delete main_simpi;
}

/******************Simpi Functions*************************/
simpi::simpi(int _id, int _num_workers)
{
  id = _id;
  num_workers = _num_workers;
  size_t synchObjectSize =
      sizeof(synch_object) + sizeof(int) * (num_workers + 1);
  int fd = shm_open(SYNCH_OBJECT_MEM_NAME, O_RDWR, 0777);
  if (fd == -1) {
    perror("Unable to shm_open synch_object: ");
    exit(1);
  }
  shm_fd = fd;
  synch_info = (synch_object*)mmap(NULL, synchObjectSize,
                                   PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
  if (synch_info == MAP_FAILED) {
    perror("Unable to mmap synch_info: ");
    exit(1);
  }
}
simpi::~simpi()
{
  for (std::pair<std::string, matrix_metadata> matrix : matrix_info) {
    free_matrix(matrix.second.unique_id);
  }
  shm_unlink(SYNCH_OBJECT_MEM_NAME);
  close(shm_fd);
}
void simpi::synch()
{
  int* ready = synch_info->ready;
  int synchid = ready[num_workers] + 1;
  ready[id] = synchid;
  int breakout = 0;
  while (1) {
    breakout = 1;
    for (int i = 0; i < num_workers; i++) {
      if (ready[i] < synchid) {
        breakout = 0;
        break;
      }
    }
    if (breakout == 1) {
      ready[num_workers] = synchid;
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
matrix::matrix(int x, int y)  // constructor
{
  // use main_simp init the matrix for all processes. The id is also in simp
  std::pair<std::string, double*> pass_back(main_simpi->create_matrix(x, y));
  unique_id = pass_back.first;
  arr = pass_back.second;
  xdim = x;
  ydim = y;
}
matrix::~matrix()  // destructor
{
  // use main_simpi for getting rid of the mem and unlink stuff
  main_simpi->free_matrix(unique_id);
}

std::ostream& operator<<(std::ostream& out, const matrix& m)
{
  if (main_simpi->get_id() == 0) {
    for (int i = 0; i < m.xdim; i++) {
      out << "\n";
      for (int j = 0; j < m.ydim; j++) {
        out << m.arr[i + j * m.xdim];
        out << ", ";
      }
    }
    out << "\n";
    return out;
  }
  else {
    return out;
  }
}

void matrix::inverse(matrix* inverse)
  {
    if (get_x() != get_y()) {
      std::cout << "Invalid Matrix";
      exit(1);
    }


    //matrix* inverse = new matrix(get_x(), get_y());
    matrix* adj = new matrix(get_x(), get_y());

    // Find determinant of A[][]
    int det = determinant(arr, get_x(), get_y());
    if (det == 0) {
      std::cout << "Singular matrix, can't find its inverse";
      return;
    }
    std::cout << "Determinant is: " << det << '\n';

    main_simpi->synch();

    adjoint(arr, adj->arr, xdim, main_simpi->get_id(),
            main_simpi->get_num_workers());
    main_simpi->synch();


    int num_processes = main_simpi->get_num_workers();
    int parID = main_simpi->get_id();
    int n = get_x();
    if (num_processes > n) {
        num_processes = n;
    }

    int rpp = n/num_processes;
    int start = main_simpi->get_id() * rpp;
    int end = start + rpp;

    for (int i = start; i < end; i++)
      for (int j = 0; j < get_x(); j++) {
        inverse->get(i, j) = (adj->get(i, j)) / double(det);
      }
    if (n%num_processes != 0) {
        int leftover = n%num_processes;
        if (parID < leftover) {
            parID += (n-leftover);
            int start = parID;
            int end = start + 1;
            for (int i = start; i< end; i++){
                for (int j = 0; j<n; j++) {
                    inverse->get(i,j) = (adj->get(i,j))/double(det);
                }
            }
        }
    }
    main_simpi->synch();
    return;
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

  int num_processes = par_count;
  int parID = par_id;
  int n = order;
  if (num_processes > n) {
    num_processes = n;
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
  if (n%num_processes != 0) {
    int leftover = n%num_processes;
    if (parID < leftover) {
        parID += (n-leftover);
        int start = parID;
        int end = start + 1;
        for (int i = 0; i <order; i++) {
            for (int j = start; j<end; j++) {
                getCofactor(A, temp, i, j, order, order);
                sign = ((i+j)%2==0)? 1: -1; 
                adj[j+i*order] = (sign)*(determinant(temp, order-1, order));


            }
        }
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

void matrix::luDecomposition(matrix* lower, matrix* upper) {
  std::cout << "Code is running!" << '\n';

  if (get_x() != get_y()) {
    std::cout << "Invalid Matrix";
    exit(1);
  }

  for (int i = 0; i < get_x(); i++) {
    int num_processes = main_simpi->get_num_workers();
    int parID = main_simpi->get_id();
    int total = xdim - i;
    if (num_processes > total) {
      num_processes = total;
    }
    int rpp = total/num_processes;
    int start = rpp*main_simpi->get_id() + i;
    int end = start + rpp;
    // std::cout << "start: " << start << '\n';
    // std::cout << "end: " << end << '\n';
    // std::cout << "rpp: " << rpp << '\n' << '\n';

    // Upper Triangular 
    for (int k = start; k< end; k++) {
      if (k>=get_x()) {
        break;
      }
      // Summation of L(i, j) * U(j, k) 
      float sum = 0;
      for (int j = 0; j < i; j++) 
        sum += (lower->get(i,j)* upper->get(j,k)); 

      // Evaluating U(i, k) 
      upper->get(i,k) = get(i,k) - sum; 

    }
    if (total%num_processes != 0) {
      int leftover = total%num_processes;
      if (parID < leftover) {
        parID += (xdim-leftover);
        int start = parID;
        std::cout << "leftover row being worked on: " << start << '\n'<<'\n';

        int end = start + 1;
        for (int a = start; a<end; a++) {
          // Summation of L(i, j) * U(j, k) 
          float sum = 0;
          for (int j= 0; j<i; j++) 
            sum+= (lower->get(i,j)*upper->get(j,a));
          // Evaluating U(i, k) 
          upper->get(i,a) = get(i,a) - sum;
        }
      }
    }

    main_simpi->synch();

    total = get_x() - i;
    parID = main_simpi->get_id();
    num_processes = main_simpi->get_num_workers();

    for (int k = start; k<end; k++) {
      if (k>=get_x()) {
        break;
      }
      if (i == k) {
        lower->get(i,i) = 1; // Diagonal as 1 
      }
      else {
        // Summation of L(k, j) * U(j, i)
        float sum = 0;
        for (int j = 0; j<i; j++) 
          sum += (lower->get(k,j) * upper->get(j,i));
        // Evaluating L(k, i)
        lower->get(k,i) = ((get(k,i)-sum) / upper->get(i,i));
      }
    }
    if (total%num_processes != 0) {
      int leftover = total%num_processes;
      if (parID < leftover) {
        parID += (get_x()-leftover);
        int start = parID;
        int end = start + 1;
        for (int a = start; a<end; a++) {
          if (i == a)
            lower->get(i,i) = 1;
          else {
            // Summation of L(k, j) * U(j, i) 
            float sum = 0;
            for (int j = 0; j<i; j++) 
              sum+= (lower->get(a,j) * upper->get(j,i));

            // Evaluating L(k, i) 
            lower->get(a,i) = (get(a,i)-sum) / upper->get(i,i);
          }
        }
      }
    }
    main_simpi->synch();




  }
return;
}

void matrix::inverseLU(matrix* inv) {
  //Check if matrix is square
  if (get_x() != get_y()) {
    std::cout << "Invalid Matrix";
    exit(1);
  }
  // int det = determinant(arr, get_x(), get_y());
  // if (det == 0) {
  //   std::cout << "Singular matrix, can't find its inverse";
  //   exit(1);
  // }
  // std::cout<<"Stage 1";


  //Solve for lower and upper matrices
  matrix* upper = new matrix(get_x(), get_y());
  matrix* lower = new matrix(get_x(), get_y());
  luDecomposition(lower,upper);
  std::cout<<"Stage 2";
  main_simpi->synch();

  //Create Identity nxn Matrix
  matrix* identity = new matrix(get_x(), get_y());
  if (main_simpi->get_id() == 0) {
    for (int i = 0; i<get_x(); i++) {
      for (int j = 0; j<get_x(); j++) {
        if (i == j) {
          identity->get(i,j) = 1;
        }
        else {
          identity->get(i,j) = 0;                
        }
      }
    }
  }

  main_simpi->synch();
  std::cout<<"Stage 3";
  // Calculate columns per parallel process
  int cpp = get_x()/main_simpi->get_num_workers();
  int start = cpp*main_simpi->get_id();
  int end = start + cpp;

  // Initialize necessary arrays for future calculations
  // Each array is local to its own process
  float identity_col[get_x()];
  float z_col[get_x()];
  float soln_col[get_x()];

  for (int a = start; a<end; a++) {

    //Get individual columns of identity matrix
    for (int b = 0; b<get_x(); b++) {
      identity_col[b] = identity->get(b,a);
    }

    //Reset Z column to solve for again
    for (int d = 0; d<get_x(); d++) {
      z_col[d] = 0;
    }
  
    //Solve LZ = I
    (*lower).forward_substitution(identity_col, z_col);

    //Reset X column to solve for again
    for (int d = 0; d<get_x(); d++) {
      soln_col[d] = 0;
    }

    //Solve UX = Z
    (*upper).backward_substitution(z_col, soln_col);

    //Input X column to corresponding columnn in final inverse matrix
    for (int c = 0; c<get_x(); c++) {
      inv->get(c,a) = soln_col[c];
    }
  }
  int parID = main_simpi->get_id();
  if (get_x()%main_simpi->get_num_workers() != 0) {
    int leftover = get_x() % main_simpi->get_num_workers();
    if (parID < leftover) {
      parID += (get_x()-leftover);
      int start = parID;
      int end = start + 1;
      for (int a = start; a<end; a++) {
        //Get individual columns of identity matrix
        for (int b = 0; b<get_x(); b++) {
          identity_col[b] = identity->get(b,a);
        }

        //Reset Z column to solve for again
        for (int d = 0; d<get_x(); d++) {
          z_col[d] = 0;
        }
      
        //Solve LZ = I
        (*lower).forward_substitution(identity_col, z_col);

        //Reset X column to solve for again
        for (int d = 0; d<get_x(); d++) {
          soln_col[d] = 0;
        }

        //Solve UX = Z
        (*upper).backward_substitution(z_col, soln_col);

        //Input X column to corresponding columnn in final inverse matrix
        for (int c = 0; c<get_x(); c++) {
          inv->get(c,a) = soln_col[c];
        }
      }
    }
  }  

  main_simpi->synch();
  return;

}

void matrix::forward_substitution(float *b, float* x)
{
    double suma;
    for(int i=0; i < get_x(); i=i+1)
    {
        suma = 0;
        for(int j=0;j<i;j=j+1)
            suma = suma+ get(i,j) * x[j];

        x[i] = (b[i]-suma)/get(i,i);
    }
}

void matrix::backward_substitution(float* b, float* x)
{
    double suma;
    for(int i= get_x()-1; i>=0; i=i-1)
    {
        suma=0;
        for(int j = get_x() - 1; j > i; j =j-1)
            suma = suma + get(i,j) * x[j];
 
        x[i]=(b[i] - suma)/get(i,i);
    }
}


void matrix::jacobi(vector* constants, vector* solution) {

    int processCount = main_simpi->get_num_workers();
    int id = main_simpi->get_id();

    vector *prev = new vector(constants->get_size()); // shared mem containing a copy of values
    //vector *solution = new vector(*main_simpi, constants->get_size()); // shared mem containing actual calculated values

    matrix* saveEq = new matrix(get_x(), get_y() + 1); // save equations from modification
    vector* saveConst = new vector(constants->get_size()); // saves input vector


    int work = constants->get_size() / processCount;

    int i, j, k;
    int start = id * work;
    int end = start + work;

    main_simpi->synch();

    //Save Matrix and Vector
    for (i = start; i < end; i++) {
        for (j = 0; j < get_y(); j++) {
            saveEq->get(i,j) = get(i, j);
        }
        saveEq->get(i,get_y() + 1) = constants->get(i);
        saveConst->get(i) = constants->get(i);
    }

    //synch, wait for all process before solving
    main_simpi->synch();

    //setup, switch var coefficient with row solution, and divide by coefficient
    for (i = start; i < end; i++) {
        double temp = get(i, i);
        get(i, i) = constants->get(i);
        constants->get(i) = temp;
        for (j = 0; j < get_y(); j++) {
            if (j != i) {
                get(i, j) *= -1;
            }
            get(i, j) /= constants->get(i);
            //main_simpi->synch();
        }
        prev->get(i) = 1.0;
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
            } else {
                rowSum += (get(i, j) * prev->get(j));
            }
            //main_simpi->synch();
        }
        solution->get(i) = rowSum;
        main_simpi->synch();
    }

    //wait for all processes before repeating iterations with calculated results
    //main_simpi->synch();

    for (k = 0; k < 1000; k++)
    {
        for (i = start; i < end; i++)
        {
            //save prev value for comparision
            prev->get(i) = solution->get(i);
            main_simpi->synch();
        }
        for (i = start; i < end; i++)
        {
            double rowSum = 0;
            for (j = 0; j < get_y(); j++)
            {
                if (j == i) {
                    rowSum += get(i, j);
                } else {
                    rowSum += (get(i, j) * prev->get(j));
                }
            }
            solution->get(i) = rowSum;
            //main_simpi->synch();
        }
        //wait at end of each loop for all processes before beginning next iteration
        main_simpi->synch();
    }
    main_simpi->synch();
    //restore original matrix and vector
    for (i = start; i < end; i++) {
        for (j = 0; j < get_y(); j++) {
            get(i, j) = saveEq->get(i,j);
        }
        constants->get(i) = saveConst->get(i);
    }
    //wait for all processes before returning solution vector
    main_simpi->synch();
    //return solution;
    return;
}

bool matrix::isDiagonallyDominant()
{
        for(int i = 0; i < get_x(); i ++)
        {
            double sq;
            double rest = 0;
            for(int j = 0; j < get_y(); j++)
            {
                if(i==j)
                {
                    sq = get(i,j);
                }
                else
                {
                    rest+=get(i,j);
                }
            }
            if (sq < rest)
            {
                return false;
            }
        }
        return true;
}

void matrix::solveSystem(vector *constants, vector* solution)
{
    bool dd = isDiagonallyDominant();
    main_simpi->synch();
    // hi
    if (dd)
    {
        jacobi(constants, solution);
    }
    else
    {
        failSafe(constants, solution);
    }
    main_simpi->synch();
}

void matrix::failSafe(vector* constants, vector* solution)
{
    matrix* inv = new matrix(get_x(), get_y());
    inverse(inv);
    main_simpi->synch();
    if (main_simpi->get_id() == 0)
    {
      int n = constants->get_size();
      double sol;
      for(int i = 0; i < n; i++)
      {
          sol = 0;
          for(int j = 0; j < n; j++)
          {
              sol += (inv->get(i, j)*constants->get(j));
          }
          solution->get(i) = sol;
      }
    }
    main_simpi->synch();
    return;
}


/******************Vector Functions*************************/
vector::vector(int a)
{
  // use simp and init the matrix for all processes. The id is also in simp
  std::pair<std::string, double*> pass_back(main_simpi->create_matrix(1, a));
  unique_id = pass_back.first;
  arr = pass_back.second;
  dim = a;
}
vector::~vector()  // destructor
{
  // use main_simpi for getting rid of the mem and unlink stuff
  main_simpi->free_matrix(unique_id);
}


  // if (main_simpi->get_id() == 0) {
  //   for (int i = 0; i < xdim; i++) {
  //     std::cout << "\n";
  //     for (int j = 0; j < ydim; j++) {
  //       std::cout << std::fixed << std::setprecision(7) << identity->get(i,j);
  //       std::cout << ", ";
  //     }
  //   }
  // }