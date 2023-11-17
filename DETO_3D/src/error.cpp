#include "error.h"

using namespace DETO_NS;

// ---------------------------------------------------------------
// Initialize class
Error::Error(DETO *deto) : Pointers(deto)
{
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    
    if (me == MASTER) fprintf(screen,"Generating error class\n");   
}


// ---------------------------------------------------------------
// Class destructor
Error::~Error()
{
    
}


// ---------------------------------------------------------------
// This produces a simple error taking care of turning off all the parallel processes
void Error::errsimple(string msg)
{
    if (me==MASTER) fprintf(screen,"\n---------------------------------------------------------------\n");
    if (me==MASTER) fprintf(screen,"%s",msg.c_str());
    if (me==MASTER) fprintf(screen,"\n---------------------------------------------------------------\n");
    MPI_Finalize();
    exit(1);
}
