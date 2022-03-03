#include "simulations.h"
#include <sstream>
#include "error.h"
//#include "universe.h"

//#include "chemistry.h"
//#include "store.h"
//#include "lammpsIO.h"
//#include "error.h"
/*#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
*/

#include <string.h>

using namespace DETO_NS;

Simulations::Simulations(DETO *deto) : Pointers(deto)
{

    // The inputcprs class is run by all processors in COMM_WORLD. This reads the id of the processor
    MPI_Comm_rank(MPI_COMM_WORLD, &me);


}



// ---------------------------------------------------------------
// Class destructor
Simulations::~Simulations()
{
    
}






// ---------------------------------------------------------------
// Printing info about the inputcprs class (possibly useful for debugging)
void Simulations::printall()
{
	fprintf(screen,"\n---------ALL ABOUT Simulations----------\n");
	//fprintf(screen,"inputcprs filename =  %s\n",fname.c_str());
	fprintf(screen,"---------------------------------------\n\n");
	
}
