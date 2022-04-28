#include "optimize.h"
#include <sstream>
#include "error.h"
#include "lammpsIO.h"
#include "universe.h"

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
    bool found_nmat = false;
    // bool found_nmat = false; Depreciated
    bool material_set = false; // boolian to determine if the user has defined material in chi map or set homogenious material
    bool type_set = false;
    bool chi_set = false;

    if (!mapFile.is_open()) {
        err_msg = "ERROR: cannot read file \"" + mapfname + "\"";
        error->errsimple(err_msg);
    }
    else {
        MPI_Barrier(MPI_COMM_WORLD);
        // READ chi_map file
        while (!mapFile.eof()) {
            MPI_Barrier(MPI_COMM_WORLD);
            std::getline(mapFile, read_string);
            if (!read_string.empty()) {
                std::istringstream lss(read_string);
                while (lss >> word) {
                    //First it is necessary to find both the num_chi and num_mater commands
                    if (strncmp(word.c_str(), "#", 1) == 0) break;
                    else if (strcmp(word.c_str(), "num_mat") == 0) {
                        lss >> nmat;
                        found_nmat = true;
                    }
                    // Now we fill in the properties for the chi_map
                    else if (strcmp(word.c_str(), "PROPERTIES:") == 0) {
                        if (found_nmat == true) {
                            // Check and insert properties
                            std::vector<std::string> check{"chi","material","type"};
                            for(int i=0; i<3; i++) {
                                lss >> word;
                                if(word != check[i]) {
                                    err_msg = "ERROR: chi map must contain PROPERTIES: chi material type in this order";
                                    error->errsimple(err_msg);
                                }
                            }
                            while (lss >> word) {
                                chi_map.properties.push_back(word);
                            }
                            // Populate chi_map
                            while (!mapFile.eof()) {
                                std::getline(mapFile, read_string);
                                read_string = read_string.substr(0, read_string.find("#"));
                                if (!read_string.empty()) {
                                    std::istringstream lss(read_string);
                                    double tempchi,temp;
                                    std::string tempmat;
                                    int temptype;
                                    tot_nchi = 0;
                                    lss >> tempchi >> tempmat >> temptype;
                                    bool mat_found =false;
                                    for(int mat=0; mat<chi_map.material.size(); mat++) {
                                        if(std::strcmp(chi_map.material[mat].c_str(),tempmat.c_str()) == 0) {
                                            mat_found = true;
                                            tot_nchi += 1;
                                            chi_map.nchi[mat] += 1;
                                            chi_map.chis[mat].push_back(tempchi);
                                            chi_map.types[mat].push_back(temptype);
                                            for(int prop=0; prop < chi_map.properties.size(); prop++) {
                                                lss >> temp;
                                                chi_map.values[mat][prop].push_back(temp);
                                            }
                                        }
                                    }
                                    if(mat_found == false) {
                                        tot_nchi += 1;
                                        chi_map.nchi.push_back(1);
                                        chi_map.material.push_back(tempmat);
                                        chi_map.chis.push_back(std::vector<double>{tempchi});
                                        chi_map.types.push_back(std::vector<int>{temptype});
                                        chi_map.values.push_back(std::vector<std::vector<double>>());
                                        for(int prop=0; prop < chi_map.properties.size(); prop++) {
                                            lss >> temp;
                                            chi_map.values[chi_map.values.size()-1].push_back(std::vector<double>{temp});
                                        }    
                                    }
                                }
                            }
                            if(nmat != chi_map.material.size()) {
                                err_msg = "ERROR: Specified number of materials was " + std::to_string(nmat) + " but found " + std::to_string(chi_map.material.size());
                                error->errsimple(err_msg);
                            }
                            // Sort chi_map
                            struct Chi_map chi_map_sorted = chi_map;
                            for(int i=0; i<nmat; i++) {
                                std::sort(chi_map_sorted.chis[i].begin(), chi_map_sorted.chis[i].end());
                                for(int srt_j=0; srt_j<chi_map.nchi[i]; srt_j++) {
                                    for(int unsrt_j=0; unsrt_j<chi_map.nchi[i]; unsrt_j++) {
                                        if(chi_map_sorted.chis[i][srt_j] == chi_map.chis[i][unsrt_j]) {
                                            chi_map_sorted.types[i][srt_j] = chi_map.types[i][unsrt_j];
                                            for(int k=0; k< chi_map.properties.size(); k++) {
                                                chi_map_sorted.values[i][k][srt_j] = chi_map.values[i][k][unsrt_j];
                                            }
                                        }
                                    }
                                }
                            }
                            chi_map = chi_map_sorted;
                            // Compute chi_min chi_max and chi_avg using per material sorted chi
                            for(int i=0; i<nmat; i++) {
                                chi_map.chi_max.push_back(chi_map.chis[i][chi_map.nchi[i]-1]);
                                chi_map.chi_min.push_back(chi_map.chis[i][0]);
                                double chi_sum = 0;
                                for(int j=0; j<chi_map.nchi[i]; j++) {
                                    chi_sum += chi_map.chis[i][j];
                                }
                                chi_map.chi_avg.push_back(chi_sum/chi_map.nchi[i]);
                            }
                            // Initalise constraints
                            for(int i=0; i<nmat; i++){
                                vol_constraintYN.push_back(false);
                                vol_constraint.push_back(0);
                                constraint_method.push_back("NULL");
                                local_vol_constraintYN.push_back(false);
                                local_vol_constraint.push_back(0);
                                local_vol_radius.push_back(0);
                                local_constraint_method.push_back("NULL");
                            }
                        }
                        else {
                            err_msg = "ERROR: num_mat not found in map_chi file before specifying table of chi's and related quantities to set";
                            error->errsimple(err_msg);
                        }
                    }
                    else {
                        err_msg = "ERROR: Unrecognized command in map_chi file: " + read_string;
                        error->errsimple(err_msg);
                    }
                }
            }
        }
        mapFile.close();
    }
}


// ---------------------------------------------------------------
void Optimize::add_constraint(std::string read_string)
{   
    std::string material,constraint_type,method;
    double constraint, radius;
    std::istringstream lss(read_string);
    lss >> material >> constraint_type >> method >> constraint;
    if(constraint_type == "totchi") {
        for(int i=0; i<nmat; i++) {
            if(chi_map.material[i] == material) {
                vol_constraintYN[i] = true;
                vol_constraint[i] = constraint;
                constraint_method[i] = method;
            }
        }
    }
    else if(constraint_type == "local_totchi") {
        lss >> radius;
        for(int i=0; i<nmat; i++) {
            if(chi_map.material[i] == material) {
                local_vol_constraintYN[i] = true;
                local_vol_constraint[i] = constraint;
                local_vol_radius[i] = radius;
                local_constraint_method[i] = method;
            }
        }
    }
    else {
        err_msg = "Unrecognized constraint type specified in: " + read_string;
        error->errsimple(err_msg);
    }
}


// ---------------------------------------------------------------
// function initializing values of chi from chi_map
void Optimize::initialize_chi()
{
    std::string tolmp;
    tolmp = "compute tempID all property/atom id"; // temp compute to get id of all atoms
    lammpsIO->lammpsdo(tolmp);
    
    tolmp = "compute tempType all property/atom type";  // temp compute reading type per atom
    lammpsIO->lammpsdo(tolmp);
    
    tolmp = "dump tdID all custom 1 dump.temp_"+universe->SCnames[universe->color]+" c_tempID c_tempType x y z";    //temp dump to update variables and computes
    lammpsIO->lammpsdo(tolmp);
    lammpsIO->lammpsdo("run 1");     // a run1 in lammps to dump the temp and so prepare variables and computes

    MPI_Barrier(MPI_COMM_WORLD);
    
    // extracting unsorted atom ids
    natoms = (int)lammpsIO->lmp->atom->natoms; // total number of atoms in lammps, all groups all fixes)
    nlocal = (int)lammpsIO->lmp->atom->nlocal; // number of atoms in current processore, all types all fixes)
    
    aID = (double *)lammps_extract_compute(lammpsIO->lmp,(char *) "tempID",1,1);
      // NB: each processor in subcom pulls out the ID of their atoms. We will put them all into a single vector, IDuns, to be managed by the submaster. The way that seems to work is to scan aID of each processor looking for the first nlocal atoms with non-zero id. Ids after those are random. The first nonzero nlocal ids are passed to the submaster, which eventually sorts them.
    atype = (double *)lammps_extract_compute(lammpsIO->lmp,(char *) "tempType",1,1);
    
    
    int tID[nlocal], ttype[nlocal];
    for(int i=0; i<nlocal; i++) {
        if (aID[i]<1){
            err_msg = "ERROR: zero or negative atom ID recorded from LAMMPS. If this happens, LAMMPS is probably creating non-contiguous ID and type vectors per processor. Change DETO code to cope with that...";
            error->errsimple(err_msg);
        }
        tID[i] = (int) aID[i];
        ttype[i] = (int) atype[i];
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    //MPI send size of each tID to submaster. Submaster creates array with enough space and assigns positions to accept IDs.  All procs then send tIDarr to sub master
    nploc = universe->SCnp[universe->color];
    key = universe->key;
    nID_each = new int[nploc];
    nID_each[key] = nlocal;   //each processor in subcomm records its number of atoms (nlocal) at location = key of its nID_each
    
    if (key>0) {
        int dest = 0;
        MPI_Send(&nID_each[key], 1, MPI_INT, dest, 1, (universe->subcomm));
    }
    if (key==0 && nploc>1) {
        for (int source=1; source<nploc; source++) {
            MPI_Recv(&nID_each[source], 1, MPI_INT, source, 1, (universe->subcomm), &status);
        }
    }
    
    // pass local ID arrays to submaster
    IDpos = new int[nploc];  // position of local tID array in submaster's unsorted list of IDs
    if (key==0) {
        IDpos[0]=0;
        for (int i=1; i<nploc; i++) {
            IDpos[i] = IDpos[i-1] + nID_each[i-1];
        }
    }
    // MPI send ID's and types from each processor to the submaster. 
    if (key==0){
        IDuns = new int[natoms];
        typeuns = new int[natoms];
    }
    if (key>0) {
        int dest = 0;
        MPI_Send(&tID[0], nlocal, MPI_INT, dest, 1, (universe->subcomm));
        MPI_Send(&ttype[0], nlocal, MPI_INT, dest, 2, (universe->subcomm));
    }
    if (key==0) {
        for (int i=0; i<nID_each[0]; i++) {
            IDuns[i] = tID[i];
            typeuns[i] = ttype[i];
        }
        for (int source=1; source<nploc; source++) {
            MPI_Recv(&IDuns[IDpos[source]], nID_each[source], MPI_INT, source, 1, (universe->subcomm), &status);
            MPI_Recv(&typeuns[IDpos[source]], nID_each[source], MPI_INT, source, 2, (universe->subcomm), &status);
        }
    }

    // DEBUG: delete this before running final version of the code
    // sleep(me);
    // if (key==0) {
    //     fprintf(screen,"\nID and type constructed at submaster\n-------------------\n");
    //     for(int i=0; i<natoms; i++) {
    //         fprintf(screen,"PROC %d ID: %d type %d\n",me,IDuns[i],typeuns[i]);
    //     }
    // }

    
    //  initialise a chi vector on the submaster
    if(key == 0) {
        // todo: compute chi_map.chi_avg and chi_map.chi_max
        for(int i=0; i<natoms; i++) {
            bool type_found = false;
            for(int j=0; j<nmat; j++) {
                for(int k=0; k<chi_map.types[j].size(); k++) {
                    if(typeuns[i] == chi_map.types[j][k]) {
                        chi.push_back(chi_map.chis[j][k]);
                        type_found = true;
                    }
                }
            }
            if(type_found == false) {
                chi.push_back(2);    // todo: push_back(fabs(chi_map.chi_avg)+chi_map.chi_max+fabs(chi_map.chi_max)).... do this per-material use greatest chi_map.chi_avg and chi_map.chi_max between all materials
            }
        }
        if(me == MASTER) {
            fprintf(screen,"\n------------------\nUnnormalised chi\n-------------------\n\n");
            for (int i=0; i<natoms; i++) {
                fprintf(screen,"%f\n",chi[i]);
            }
        }
        //Enforce Volume constraint, or local volume constraint if any are defined here
        constrain_vol();
        constrain_local_vol();    //todo: to be implemented... not in a rush though
        
    }
    if(me == MASTER) {
        fprintf(screen,"\n------------------\nNormalised chi\n-------------------\n\n");
        for (int i=0; i<natoms; i++) {
            fprintf(screen,"%f\n",chi[i]);
        }
    }

    // TODO: processors must communicate their IDS to the other and local MASTER (viz. proc with key == 0 in current bcomm) must assemble the local IDs and types into its global vectors aID, atype).  Then key == 0 must associate chi to each type.
    // for each lammps type assign corresponding chi from chi_map
    
    
    
    // In this class we need a per-particle vector of chi values, so we need to extract from LAMMPS all the particles with type included in the chi_map list. Actually, we need a chi vector for each material tpye (also, listed in chi_map file)
    // TODO: when reading chi_map also read number of materials. If any of the specified material ID in chi_map does not fit the number (e.g., you specify material 0 and 1 but you had given only num_mater = 1) then produce error
    // When defining the chi vector, you can make as long as all atoms in lammps, but assign 2 x max_chi (from chi_map) to the atoms whose type is not included din the chi_map (so they wil visualize as chi = max, i.e. solid, in OVITO)
    
    // delete IDuns;
    // delete IDpos;
    // delete nID_each;

}
    

// ---------------------------------------------------------------
// Function to constrain the total chi in the simulation, does not contain material
void Optimize::constrain_vol() // This function needs to have different modes for different volume constraints
//It is best not to use this function with negative chi as it is untested
{
    for(int i=0; i<nmat; i++) {
        double chi_sum = 0.;
        if(vol_constraintYN[i] == true){
            if(std::strcmp(constraint_method[i].c_str(),"scale") == 0) {
                double l1 = 0.;  // chi_min - chi_map.chi_avg/10
                double l2 = 1000;//chi_map.chi_max[i]*(1/chi_map.chi_min[i]);   // chi_map.chi_max //maybe this constraint only works for positive chi
                double lmid; //different constraint types can allow us to accomidate more
                double chi_constrained[natoms];
                while (l2-l1 > 1e-8){
                    chi_sum = 0.;
                    lmid = 0.5*(l2+l1);
                    for (int j=0; j<natoms; j++) {
                        chi_constrained[j] = chi[j]*lmid;
                        if(chi_constrained[j] > chi_map.chi_max[i]) chi_constrained[j] = chi_map.chi_max[i];
                        if(chi_constrained[j] < chi_map.chi_min[i]) chi_constrained[j] = chi_map.chi_min[i];
                        chi_sum += chi_constrained[j];
                    }
                    if ((chi_sum - vol_constraint[i]*(double)natoms) > 0) l2 = lmid;
                    else l1 = lmid;
                }
                for (int j=0; j<natoms; j++) {
                    chi[j] = chi_constrained[j];
                }
            }
            else if(std::strcmp(constraint_method[i].c_str(),"shift") == 0) {
                double l1 = -chi_map.chi_max[i];
                double l2 = chi_map.chi_max[i];
                double lmid;
                double chi_constrained[natoms];
                while (l2-l1 > 1e-8){
                    chi_sum = 0.;
                    lmid = 0.5*(l2+l1);
                    for (int j=0; j<natoms; j++) {
                        chi_constrained[j] = chi[j]+lmid;
                        if(chi_constrained[j] > chi_map.chi_max[i]) chi_constrained[j] = chi_map.chi_max[i];
                        if(chi_constrained[j] < chi_map.chi_min[i]) chi_constrained[j] = chi_map.chi_min[i];
                        chi_sum += chi_constrained[j];
                    }
                    if ((chi_sum - vol_constraint[i]*(double)natoms) > 0) l2 = lmid;
                    else l1 = lmid;
                }
                for (int j=0; j<natoms; j++) {
                    chi[j] = chi_constrained[j];
                }
            }
            else {
                err_msg = "ERROR: Unrecognised constraint method";
                error->errsimple(err_msg);
            }
            if ((chi_sum - vol_constraint[i]*(double)natoms) > 0.1 || (chi_sum - vol_constraint[i]*(double)natoms) < -0.1) {
                err_msg = "ERROR: was not able to apply volume constraint sucessfully";
                error->errsimple(err_msg);
            }
        }
    }
}


// ---------------------------------------------------------------
// running the optimization
void Optimize::constrain_local_vol()
{

}


// ---------------------------------------------------------------
// running the optimization
void Optimize::optrun()
{
    // initialize chi values from types and set other type-related quantities in the chi_map file too

    initialize_chi();
    
   
    
    
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
    for(int i=0; i<nmat; i++) {
        fprintf(screen,"\nMaterial: %s\n", chi_map.material[i].c_str());
        fprintf(screen,"Chi type ");
        for(int j=0; j<chi_map.properties.size(); j++) {
            fprintf(screen,"%s ",chi_map.properties[j].c_str());
        }
        fprintf(screen, "\n-------------------------\n");
        for (int j = 0; j < chi_map.nchi[i]; j++) {
            fprintf(screen,"%f %d ",chi_map.chis[i][j],chi_map.types[i][j]);
            for(int k=0; k<chi_map.properties.size(); k++) {
                fprintf(screen,"%.2f ",chi_map.values[i][k][j]);
            }
            fprintf(screen,"\n");
        }
    }
    //Print constraints
    fprintf(screen,"\nGloal volume constraints: \n");
    for(int i=0; i<nmat; i++) {
        fprintf(screen,"%s  f=%f ",chi_map.material[i].c_str(),vol_constraint[i]);
    }
    fprintf(screen,"\n");
    fprintf(screen,"Local volume constraints: \n");
    for(int i=0; i<nmat; i++) {
        fprintf(screen,"%s  f=%f radius=%f ",chi_map.material[i].c_str(),local_vol_constraint[i],local_vol_radius[i]);
    }
    
    fprintf(screen, "\n---------------------------------------\n\n");
}
