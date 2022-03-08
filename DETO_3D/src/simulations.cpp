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
void Simulations::add(std::string read_string)
{
// record name in a vector
    // record type in a vector and, depending on type, record associated parameters (none if "run", some if "cstgs"
    // record repeat yes/no in a boolean vector
    // record repeat attributes in a vector (empty entry if repeat = no, otherwsie whatever is needed for the repeat (I think a file to load something from)
    
}


// ---------------------------------------------------------------
void Simulations::add_attribute(std::string read_string)
{
// find simulation id matching specified name
    // add string to a vector of vectors of string (one list of attributes per simulation) at the mathcing simulation ID
    
}



// ---------------------------------------------------------------
// Printing info about the inputcprs class (possibly useful for debugging)
void Simulations::printall()
{
	fprintf(screen,"\n---------ALL ABOUT Simulations----------\n");
	//fprintf(screen,"inputcprs filename =  %s\n",fname.c_str());
	fprintf(screen,"---------------------------------------\n\n");
	
}
