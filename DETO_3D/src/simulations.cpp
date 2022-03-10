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
    std::string sim_name, sim_type, repeatfname;
    bool is_repeat;
    std::istringstream lss(read_string);
    lss >> sim_name >> sim_type;
    sim_names.push_back(sim_name);
    sim_types.push_back(sim_type);
    sim_attributes.push_back(std::vector<std::string>());
    // record type in a vector and, depending on type, record associated parameters (none if "run", some if "cstgs"
    // record repeat yes/no in a boolean vector
    // record repeat attributes in a vector (empty entry if repeat = no, otherwsie whatever is needed for the repeat (I think a file to load something from)
    
}


// ---------------------------------------------------------------
void Simulations::add_attribute(std::string read_string)
{
// find simulation id matching specified name
    // add string to a vector of vectors of string (one list of attributes per simulation) at the mathcing simulation ID
    std::string sim_ID;
    std::istringstream lss(read_string);
    lss >> sim_ID;
    for(int i = 0; i < sim_names.size(); i++){
        if(strcmp(sim_ID.c_str(),sim_names[i].c_str()) == 0){
            sim_attributes[i].push_back(lss.str());
            fprintf(screen,"simulations -- in %s adding attribute %s \n",sim_ID.c_str(),lss.str().c_str());
        }
    }
}



// ---------------------------------------------------------------
// Printing info about the inputcprs class (possibly useful for debugging)
void Simulations::printall()
{
	fprintf(screen,"\n---------ALL ABOUT Simulations----------\n");
	//fprintf(screen,"inputcprs filename =  %s\n",fname.c_str());
	fprintf(screen,"---------------------------------------\n\n");
	
}
