#include <sys/mman.h>
#include <sys/stat.h>       
#include <fcntl.h>
#include <unistd.h>
#include <sys/time.h>
#include <vector>
#include <string>
#include <unordered_map> 
#include <sstream>
#include "simpi.h"

int main(int argc, char* argv[]){
    
  if (argc != 3) {
        printf("Usage: ./prog2 prog1 <num workers>");
        exit(2);
    }

 char progname[100];
    //strcpy(progname, "./");
    strcpy(progname, argv[1]);

    int numWorkers = atoi(argv[2]);
    int fd = shm_open(SYNCH_OBJECT_MEM_NAME, O_RDWR | O_CREAT | O_EXCL, 0777);
	  ftruncate(fd, sizeof(synch_object));
	  synch_object* shared_mem = (synch_object*)mmap(NULL, sizeof(synch_object), PROT_READ|PROT_WRITE, MAP_SHARED, fd, 0);
    int ready[numWorkers+1];
    for(int i = 0; i < numWorkers;i++) {
      ready[i] = 0;
    } 
    shared_mem->ready = ready;

    for (int i = 0; i < numWorkers; i++) {
        char workerID[3];
        sprintf(workerID, "%d", i);
        char* args[] = {progname, workerID, NULL};
        if (fork() == 0) {
            execv(progname, args);
        }
    }
    exit(0);
}
