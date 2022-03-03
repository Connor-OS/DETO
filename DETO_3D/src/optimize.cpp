#include "optimize.h"
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

Optimize::Optimize(DETO *deto) : Pointers(deto)
{

    // The inputcprs class is run by all processors in COMM_WORLD. This reads the id of the processor
    MPI_Comm_rank(MPI_COMM_WORLD, &me);


}



// ---------------------------------------------------------------
// Class destructor
Optimize::~Optimize()
{
    
}




// ---------------------------------------------------------------
// Printing info about the inputcprs class (possibly useful for debugging)
void Optimize::read_chimap(std::string mapfname)
{
    std::ifstream mapFile(mapfname.c_str());
    bool found_nchi = false;
    
    if (!mapFile.is_open()) {
        err_msg = "ERROR: cannot read file \""+mapfname+"\"";
        error->errsimple(err_msg);
    }
    else{
        MPI_Barrier(MPI_COMM_WORLD);
       // READ FILE (all processors need to know this)
        while (!mapFile.eof()) {
            MPI_Barrier(MPI_COMM_WORLD);
            std::getline (mapFile, read_string);
            if (!read_string.empty()){
                std::istringstream lss(read_string);
                while (lss >> word) {
                    if (strncmp(word.c_str(), "#", 1) == 0) break;
                    else if (strcmp(word.c_str(), "num_chi") == 0) {
                        lss >> nchi;
                        found_nchi = true;
                    }
                    else if (strcmp(word.c_str(), "chi") == 0) {
                        if (found_nchi==true){
                            // Insert map keys
                            chi_map.insert(std::pair<std::string, std::vector<double> > (word, std::vector<double>()));
                            while (lss >> word) {
                                chi_map.insert(std::pair<std::string, std::vector<double> > (word, std::vector<double>()));
                            }
                            // Populate chi_map 
                            for(int i=0;i<nchi;i++){
                                std::getline (mapFile, read_string);
                                if (!read_string.empty()){
                                    std::istringstream lss(read_string);
                                    for(it = chi_map.begin(); it != chi_map.end(); it++){
                                        // double value;
                                        double tempval;
                                        lss >> tempval;
                                        chi_map[it->first].push_back(tempval);
                                    }
                                }
                                else{
                                    err_msg = "ERROR: Fewer than specified chi's given";
                                    error->errsimple(err_msg);
                                }
                            }
                        }
                        else {
                            err_msg = "ERROR: num_chi not found in map_chi file before specifying table of chi's and related quantities to set";
                            error->errsimple(err_msg);
                        }
                    }
                    else {
                        err_msg = "Unrecognized command in map_chi file: "+read_string;
                        error->errsimple(err_msg);
                    }
                }
 
            }
        }
        
        mapFile.close();

        if (me == MASTER) {
            fprintf(screen,"DONE Reading map chi file\n");
            // Print Chi map to screen
            fprintf(screen,"\nChi Map is:\n");
            for(it = chi_map.begin(); it != chi_map.end(); it++){
                fprintf(screen,"%s  ",it->first.c_str());
            }
            fprintf(screen,"\n-------------------------\n");
            for(int i=0;i<nchi;i++){
                for(it = chi_map.begin(); it != chi_map.end(); it++){
                    fprintf(screen,"%e\t",chi_map[it->first][i]);
                }
                fprintf(screen,"\n");
            }
            fprintf(screen,"\n");
        }
        
        
    }
}




// ---------------------------------------------------------------
// Printing info about the inputcprs class (possibly useful for debugging)
void Optimize::printall()
{
	fprintf(screen,"\n---------ALL ABOUT OPTIMIZE----------\n");
	//fprintf(screen,"inputcprs filename =  %s\n",fname.c_str());
	fprintf(screen,"---------------------------------------\n\n");
	
}
