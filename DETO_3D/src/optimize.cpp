#include "optimize.h"
#include <sstream>
#include "error.h"
#include "lammpsIO.h"
//#include "universe.h"

//#include "chemistry.h"
//#include "store.h"
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
    bool found_nmat = false;
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
                    else if (strcmp(word.c_str(), "num_mater") == 0)
                    {
                        lss >> nmat;
                        found_nmat = true;
                        for(int i=0; i<nmat; i++){
                            vol_constraint.push_back(-1.0);
                            local_vol_constraint.push_back(-1.0);
                            local_vol_radius.push_back(-1.0);
                        }
                    }
                    else if (strcmp(word.c_str(), "PROPERTIES:") == 0)
                    {
                        if (found_nchi == true && found_nmat == true)
                        {
                            // Insert map keys
                            bool material_set = false;
                            bool type_set = false;
                            bool chi_set = false;
                            while (lss >> word)
                            {
                                chi_map.properties.push_back(word);
                                chi_map.values.push_back(std::vector<double>());
                                for(int i=0; i<nchi; i++) {
                                    chi_map.values[chi_map.values.size()-1].push_back(0);
                                }
                                if (word == "chi")
                                {
                                    chi_set = true;
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
                            if (material_set == false || type_set == false || chi_set == false)
                            {
                                err_msg = "ERROR: Please specify chi, material, and type in chi map";
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


void Optimize::add_constraint(std::string read_string)
{   
    std::string constraint_type;
    int mat_ID;
    double constraint, radius;
    std::istringstream lss(read_string);
    lss >> mat_ID >> constraint_type >> constraint;
    if(constraint_type == "volume") {
        fprintf(screen,"adding vol %f to material %d\n",constraint,mat_ID);
        vol_constraint[mat_ID-1] = constraint;
    }
    else if(constraint_type == "local_volume") {
        lss >> radius;
        fprintf(screen,"adding vol %f and %f to material %d\n",constraint,radius,mat_ID);
        local_vol_constraint[mat_ID-1] = constraint;
        local_vol_radius[mat_ID-1] = radius;
    }
    else {
        err_msg = "Unrecognized constraint type specified in: " + read_string;
        error->errsimple(err_msg);
    }
}

// ---------------------------------------------------------------
// function initializing values of chi from chi_map
void Optimize::initalize_chi(int natoms)
{   
    for(int i=0; i<natoms; i++){
        auto type_index = find(chi_map.properties.begin(), chi_map.properties.end(), "type");

        auto chi_index = find(chi_map.properties.begin(), chi_map.properties.end(), "chi");
        
        // chi.push_back(chi value at this point)
    }
    //if vol constraint set. i.e target_vol:
        //apply constraint with bisecting alg
}
    


// ---------------------------------------------------------------
// running the optimization
void Optimize::optrun()
{   
    int* atomID;
    int* atomtype;
    natoms = lammpsIO->extract_natoms();
    atomID = (int*)lammpsIO->extract_atom_varaiable("id");
    atomtype = (int*)lammpsIO->extract_atom_varaiable("type");
    for(int i=0; i<natoms; i++) {
        aID.push_back(*(atomID+i));
        atype.push_back(*(atomtype+i));
    }
    // if (me == MASTER) {
    //     fprintf(screen,"\n\nAtoms in simulation: %d\n\n",natoms);
    //     for(int i=0; i<natoms; i++) {
    //         fprintf(screen,"ID: %d type: %d\n",aID[i],atype[i]);
    //         }
    // }
    // initialize chi values from types and set other type-related quantities in the chi_map file too
    // initialize_chi();
    // In this class we need a per-particle vector of chi values, so we need to extract from LAMMPS all the particles with type included in the chi_map list. Actually, we need a chi vector for each material tpye (also, listed in chi_map file)
    // TODO: when reading chi_map also read number of materials. If any of the specified material ID in chi_map does not fit the number (e.g., you specify material 0 and 1 but you had given only num_mater = 1) then produce error
    // When defining the chi vector, you can make as long as all atoms in lammps, but assign 2 x max_chi (from chi_map) to the atoms whose type is not included din the chi_map (so they wil visualize as chi = max, i.e. solid, in OVITO)
    
    
    // for certain optimization types, e.g. GA, we may have to create a first set of chi_vectors to then initite the while loop
    // initialize_chipop();
    
    // start the while loop (exit conditions may be a tolerance on chi changes or a number of steps)
    // while (optimization not converged) {
    
        // create a population of M chi_vectors (e.g. take the current chi_vector and randomly move particle to chi category above or below, or use sensitivies, or cross-breed and mutate for GA, etc...).  In doing this, we must respect constraints on chi_vectors (e.g. total chi for each material or else)
    
        // for each individual chi_vector in the population, run all the simulations and evaluate its objectives (this is to be done in parallel, so each subcomm takes care of a fracion of the M individuals
    
        // compute the overall objective associated to each individual chi vector (i.e. a combination of the single objectives coming from all the simulations for that chi vector)
    
        // pick one or more chi_vectors to generate the next population
    
    // end of while loop
    // }
    
    // some criterion to pick the winner chi_vector (e.g. the one with lowest objective value)
    
}


// ---------------------------------------------------------------
// Printing info about the inputcprs class (possibly useful for debugging)
void Optimize::printall()
{
    fprintf(screen, "\n---------ALL ABOUT OPTIMIZE----------\n");
    //Print chi map
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
    //Print constraints
    fprintf(screen,"\nGloal volume constraints: \n");
    for(int i=0; i<nmat; i++) {
        fprintf(screen,"f%d : %f ",i,vol_constraint[i]);
    }
    fprintf(screen,"\n");
    fprintf(screen,"Local volume constraints: \n");
    for(int i=0; i<nmat; i++) {
        fprintf(screen,"f%d : %f radius: %f",i,local_vol_constraint[i],local_vol_radius[i]);
    }
    
    fprintf(screen, "\n---------------------------------------\n\n");
}
