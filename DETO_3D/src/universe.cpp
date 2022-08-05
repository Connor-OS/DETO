#include "universe.h"
#include "output.h"
#include <string.h>

using namespace DETO_NS;


Universe::Universe(DETO *deto) : Pointers(deto)
{
    // The inputcprs class is run by all processors in COMM_WORLD. This reads the id of the processor
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

    
    if (me == MASTER) fprintf(screen,"Generating universe class\n");

    SCnames.clear();
    SCnp.clear();
    nsc = 0;
    color = 0;
    key = 0;


}


// ---------------------------------------------------------------
// Class destructor
Universe::~Universe()
{
    /*
     delete [] color_each;
     */
}


// ---------------------------------------------------------------
// Creates the universe of sub-communicators
void Universe::create()
{   
    if (me == MASTER) fprintf(screen,"\nCreating the universe of sub-communicators as defined by the user in the input file (or in the restart file)...\n");
    
    // Create cumulative number of processors vector
    vector<int> cumnp;
    cumnp.push_back(SCnp[0]);
    for (int i=1; i<nsc; i++) cumnp.push_back(cumnp[i-1]+SCnp[i]);
    
    // Find the subcomm to which the current processor pertains. Assign it a color (subcomm number) and key (rank in the subcomm)
    for (int i=0; i<nsc; i++){
        if (me >= cumnp[i]) color++;
    }
    if (me<cumnp[0]) key = me;
    else key = me - cumnp[color-1];
    
    //splitting the communicator
    MPI_Comm_split(MPI_COMM_WORLD, color, key, &subcomm);
    
    // Create vector with IDs of all submasters
    subMS.push_back(0);
    for (int i=1; i<nsc; i++) {
        subMS.push_back(subMS[i-1]+SCnp[i-1]);  //Issue: off by one error, check MASKE src
    }
           
    // All processors tell every other processor their color (subcomm id) for later use (when sending subcomm-specific stuff, e.g. random number seed)
    color_each = new int[nprocs];
    color_each[me] = color;
    MPI_Allgather(MPI_IN_PLACE,1,MPI_INT,color_each,1,MPI_INT,MPI_COMM_WORLD);    
}

    
// ---------------------------------------------------------------
// Printing info about the inputcprs class (possibly useful for debugging)
void Universe::printall()
{
	fprintf(screen,"\n---------ALL ABOUT UNIVERSE----------\n");
	//fprintf(screen,"inputcprs filename =  %s\n",fname.c_str());
	fprintf(screen,"---------------------------------------\n\n");
    //MPI_Comm_free(&subcomm);
	
}
