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
    std::string sim_name, sim_type, repstr, repYN;
    
    std::istringstream lss(read_string);
    lss >> sim_name >> sim_type >> repstr >> repYN;
    sim_names.push_back(sim_name);
    sim_types.push_back(sim_type);
    if (strcmp(repstr.c_str(), "repeat") != 0) {
        std::string msg = "Error: keyword \"repeat\" must follow simulation type in input script\n";
        error->errsimple(msg);
    }
    if (strcmp(repYN.c_str(), "yes") == 0){
        sim_is_repeat.push_back(1);
        std::string repfname;
        lss >> repfname;
        sim_repeat_file.push_back(repfname);
    }
    else if (strcmp(repYN.c_str(), "no") == 0) {
        sim_is_repeat.push_back(0);
        sim_repeat_file.push_back("NULL");
    }
    else {
        std::string msg = "Error: repeat must be \"yes\" or \"no\", case sensitive. Instead I found "+repYN+"\n";
        error->errsimple(msg);
    }
    
    
    // default values for additional inputs related to the simulations - will be overwritten if sim is "cstgs"
    cstgs_varname.push_back("NULL");
    cstgs_type.push_back("NULL");
    cstgs_par1.push_back(-1.);
    cstgs_par2.push_back(-1.);
    cstgs_tol.push_back(100000000);
    
    if (strcmp(sim_type.c_str(), "cstgs") == 0) {
        int pos = cstgs_varname.size()-1;
        std::string varname, tolstr, ttype;
        double par1, par2, tol;
        lss >> varname >> ttype >> par1 >> par2;
        cstgs_varname[pos] = varname;
        cstgs_type[pos] = ttype;
        cstgs_par1[pos] = par1;
        cstgs_par2[pos] = par2;
        if (strcmp(ttype.c_str(), "binary_chop") == 0) {
            lss >> tolstr >> tol;
            if (strcmp(tolstr.c_str(), "tol") != 0) {
                std::string msg = "Error: tolerance must be specified through \"tol\" keyword when using binary_chop";
                error->errsimple(msg);
            }
            cstgs_tol[pos] = tol;
        }
    }
    
    cstgs_crit.push_back("NULL");   // this will be provided separaetly by the user through dedicated "Criterion sim_ID string" command in the input scrupt
    
    sim_attributes.push_back(std::vector<std::string>());
    
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
            std::string read_string2;
            std::getline(lss, read_string2);
            sim_attributes[i].push_back(read_string2);
            if (me==MASTER){
            fprintf(screen,"simulations -- in %s adding attribute %s \n",sim_names[i].c_str(),sim_attributes[i][sim_attributes[i].size()-1].c_str());
            }
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
