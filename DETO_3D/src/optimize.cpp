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
    // struct Chi_map chi_map;

    if (!mapFile.is_open())
    {
        err_msg = "ERROR: cannot read file \"" + mapfname + "\"";
        error->errsimple(err_msg);
    }
    else
    {
        MPI_Barrier(MPI_COMM_WORLD);
        // READ FILE (all processors need to know this)
        while (!mapFile.eof())
        {
            MPI_Barrier(MPI_COMM_WORLD);
            std::getline(mapFile, read_string);
            if (!read_string.empty())
            {
                std::istringstream lss(read_string);
                while (lss >> word)
                {
                    if (strncmp(word.c_str(), "#", 1) == 0)
                        break;
                    else if (strcmp(word.c_str(), "num_chi") == 0)
                    {
                        lss >> nchi;
                        found_nchi = true;
                    }
                    else if (strcmp(word.c_str(), "PROPERTIES:") == 0)
                    {
                        if (found_nchi == true)
                        {
                            // Insert map keys
                            bool material_set = false;
                            bool type_set = false;
                            while (lss >> word)
                            {
                                chi_map.properties.push_back(word);
                                chi_map.values.push_back(std::vector<double>());
                                for(int i=0; i<nchi; i++) {
                                    chi_map.values[chi_map.values.size()-1].push_back(0);
                                }

                                if (word == "material")
                                {
                                    material_set = true;
                                }
                                if (word == "type")
                                {
                                    type_set = true;
                                }
                            }
                            if (material_set == false || type_set == false)
                            {
                                err_msg = "ERROR: Please specify material and type in chi map";
                                error->errsimple(err_msg);
                            }
                            // Populate chi_map
                            for (int i = 0; i < nchi; i++)
                            {
                                std::getline(mapFile, read_string);
                                if (!read_string.empty())
                                {
                                    std::istringstream lss(read_string);
                                    for (int j=0; j<chi_map.properties.size(); j++)
                                    {
                                        // double value;
                                        double tempval;
                                        lss >> tempval;
                                        chi_map.values[j][i] = tempval;
                                    }
                                }
                                else
                                {
                                    err_msg = "ERROR: Fewer than specified chi's given";
                                    error->errsimple(err_msg);
                                }
                            }
                        }
                        else
                        {
                            err_msg = "ERROR: num_chi not found in map_chi file before specifying table of chi's and related quantities to set";
                            error->errsimple(err_msg);
                        }
                    }
                    else
                    {
                        err_msg = "Unrecognized command in map_chi file: " + read_string;
                        error->errsimple(err_msg);
                    }
                }
            }
        }
        mapFile.close();
    }
}

void Optimize::initalize_chi()
{
    for (int i=0; i<nchi; i++) {
        for (int j=0; j<chi_map.properties.size(); j++) {
            // lammpsIO->lammpsdo("set  type %d %s %g",type,chi_map.properties[j],chi_map.values[j][i])
            if (me == MASTER) {
                fprintf(screen,"set  type %d %s %g",type,chi_map.properties[j],chi_map.values[j][i]);
            }
        }
    }
}
    


// ---------------------------------------------------------------
// Printing info about the inputcprs class (possibly useful for debugging)
void Optimize::printall()
{
    fprintf(screen, "\n---------ALL ABOUT OPTIMIZE----------\n");
    // fprintf(screen,"inputcprs filename =  %s\n",fname.c_str());
    fprintf(screen,"Chi Map\n");
    for(int i=0; i<chi_map.properties.size(); i++) {
        fprintf(screen,"%s ",chi_map.properties[i].c_str());
    }
    fprintf(screen, "\n-------------------------\n");
    for (int i = 0; i < nchi; i++) {
        for(int j=0; j<chi_map.properties.size(); j++) {
            fprintf(screen,"%.2f ",chi_map.values[j][i]);
        }
        fprintf(screen,"\n");
    }
    
    fprintf(screen, "---------------------------------------\n\n");
}
