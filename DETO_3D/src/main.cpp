#include "deto.h"
#include "mpi.h"
#include <unistd.h>   //just for the sleep() function
#include "universe.h"
#include "inputdeto.h"


using namespace DETO_NS;

int main(int argc, char **argv)
{
    
    MPI_Init(&argc,&argv);
    
#if defined(DETO_WAIT_ATTACH_RANK)
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == DETO_WAIT_ATTACH_RANK) {
#endif
#if defined(DETO_WAIT_ATTACH) || defined(DETO_WAIT_ATTACH_RANK)
    int i = 0;
    char hostname[256];
    gethostname(hostname, sizeof(hostname));
    printf("PID %d on %s ready for attach\n", getpid(), hostname);
    fflush(stdout);
    while (0 == i)
      sleep(5);
#endif
#if defined(DETO_WAIT_ATTACH_RANK)
    }
#endif
    
    
    
    // Create an object from the fundamental class (this object is essentially your whole simulation)
    DETO *deto = new DETO(argc,argv);
    
    
    
    
    

    // Read the input file (or, automatically, the restart file if there is one that is active)
    deto->inputdeto->file();
    
    
	deto -> printall();
    
    delete deto;
     
    MPI_Finalize();

}

