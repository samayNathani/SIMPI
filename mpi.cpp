#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <unistd.h>

#include <cstring>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "simpi.h"

int main(int argc, char* argv[])
{
  if (argc != 3) {
    printf("Usage: ./prog2 prog1 <num workers>");
    exit(2);
  }

  char progname[100];
  int fd;
  // strcpy(progname, "./");
  strcpy(progname, argv[1]);

  int numWorkers = atoi(argv[2]);

  // create shared mem for workers
  size_t synchObjectSize =
      sizeof(synch_object) + sizeof(int) * (numWorkers + 1);
  fd = shm_open(SYNCH_OBJECT_MEM_NAME, O_RDWR | O_CREAT | O_EXCL, 0777);
  ftruncate(fd, synchObjectSize);
  synch_object* shared_mem = (synch_object*)mmap(
      NULL, synchObjectSize, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);

  // initialize ready to zero
  for (int i = 0; i <= numWorkers; i++) {
    shared_mem->ready[i] = 0;
  }
  shared_mem->par_count = numWorkers;

  for (int i = 0; i < numWorkers; i++) {
    char workerID[3];
    std::string size = std::to_string(synchObjectSize);
    char size_cstring[size.size() + 1];
    strcpy(size_cstring, size.c_str());
    sprintf(workerID, "%d", i);
    char* args[] = {progname, workerID, size_cstring, NULL};
    if (fork() == 0) {
      execv(progname, args);
    }
  }
  exit(0);
}
