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
    fprintf(screen,"OKAY till HERE\n");
    std::ifstream mapFile(mapfname.c_str());
    bool found_nchi = false;
    bool found_nmat = false;

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
                    if (strncmp(word.c_str(), "#", 1) == 0) break;
                    else if (strcmp(word.c_str(), "num_chi") == 0) {
                        lss >> nchi;
                        found_nchi = true;
                    }
                    else if (strcmp(word.c_str(), "num_mater") == 0) {
                        lss >> nmat;
                        found_nmat = true;
                        for(int i=0; i<nmat; i++){
                            vol_constraintYN.push_back(false);
                            vol_constraint.push_back(-1.0);
                            local_vol_constraintYN.push_back(false);
                            local_vol_constraint.push_back(-1.0);
                            local_vol_radius.push_back(-1.0);
                        }
                    }
                    else if (strcmp(word.c_str(), "PROPERTIES:") == 0) {
                        if (found_nchi == true && found_nmat == true) {
                            // Insert map keys
                            bool material_set = false;
                            bool type_set = false;
                            bool chi_set = false;
                            std::vector<int> delete_pos; 
                            int pos = 0;
                            fprintf(screen,"OKAY till HERE!\n");
                            while (lss >> word) {
                                if (word == "chi") {
                                    chi_set = true;
                                    delete_pos.push_back(pos);
                                }
                                else if (word == "material") {
                                    material_set = true;
                                    delete_pos.push_back(pos);
                                }
                                else if (word == "type") {
                                    type_set = true;
                                    delete_pos.push_back(pos);
                                }
                                chi_map.properties.push_back(word);
                                chi_map.values.push_back(std::vector<double>());
                                for(int i=0; i<nchi; i++) {
                                    chi_map.values[chi_map.values.size()-1].push_back(0);
                                }
                                pos++;
                            }
                            if (material_set == false || type_set == false || chi_set == false) {
                                err_msg = "ERROR: Please specify chi, material, and type in chi map";
                                error->errsimple(err_msg);
                            }

                            fprintf(screen,"OKAY till HERE!!\n");
                            // Populate chi_map
                            for (int i = 0; i < nchi; i++) {
                                std::getline(mapFile, read_string);
                                if (!read_string.empty()) {
                                    std::istringstream lss(read_string);
                                    double tempval;
                                    int tempint;
                                    // lss >> word;
                                    for (int j=0; j<chi_map.properties.size(); j++) {
                                        if (strcmp(chi_map.properties[j].c_str(),"chi") == 0) {
                                            fprintf(screen,"found chi!");
                                            lss >> tempval;
                                            chi_map.chis.push_back(tempval);
                                        }
                                        else if (strcmp(chi_map.properties[j].c_str(),"material") == 0) {
                                            fprintf(screen,"found material!");
                                            lss >> word;
                                            chi_map.material.push_back(word);
                                        }
                                        else if (strcmp(chi_map.properties[j].c_str(),"type") == 0) {
                                            fprintf(screen,"found type!");
                                            lss >> tempint;
                                            chi_map.types.push_back(tempint);                                           
                                        }
                                        else {
                                            lss >> tempval;
                                            chi_map.values[j][i] = tempval; 
                                        }
                                    }
                                }
                                else {
                                    err_msg = "ERROR: Fewer than specified chi's given";
                                    error->errsimple(err_msg);
                                }
                                fprintf(screen,"OKAY till HERE!!!\n");
                            }
                            // // Move chi and type into separate vectors
                            sort(delete_pos.begin(), delete_pos.end(), std::greater<int>());
                            for(int i=0; i < delete_pos.size(); i++) {
                                chi_map.values.erase(chi_map.values.begin()+delete_pos[i]);
                                chi_map.properties.erase(chi_map.properties.begin()+delete_pos[i]);                                 
                            }
                        }
                        else {
                            err_msg = "ERROR: num_chi not found in map_chi file before specifying table of chi's and related quantities to set";
                            error->errsimple(err_msg);
                        }
                    }
                    else {
                        err_msg = "Unrecognized command in map_chi file: " + read_string;
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
    std::string constraint_type;
    int mat_ID;
    double constraint, radius;
    std::istringstream lss(read_string);
    lss >> mat_ID >> constraint_type >> constraint;
    if(constraint_type == "totchi") {
        vol_constraint[mat_ID-1] = constraint;
    }
    else if(constraint_type == "local_totchi") {
        lss >> radius;
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
        // todo: compute chi_avg and chi_max
        for(int i=0; i<natoms; i++) {
            bool type_found = false;
            for(int j=0; j<chi_map.types.size(); j++) {
                if(typeuns[i] == chi_map.types[j]) {
                    chi.push_back(chi_map.chis[j]);
                    type_found = true;
                }
            }
            if(type_found == false) {
                chi.push_back(2);    // todo: push_back(fabs(chi_avg)+chi_max).... do this per-material
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
    
    delete IDuns;
    delete IDpos;
    delete nID_each;

}
    

// ---------------------------------------------------------------
// running the optimization
void Optimize::constrain_vol()
{
    for(int i=0; i<nmat; i++) {
        if(vol_constraint[i] > 0){   // replace this condition with if(flag)...
            double l1 = -1.;
            double l2 = 1.;
            double lmid,chi_sum;
            double chi_constrained[natoms];
            while (l2-l1 > 1e-8){
                chi_sum = 0.;
                lmid = 0.5*(l2+l1);
                for (int j=0; j<natoms; j++) {
                    chi_constrained[j] = chi[j]+lmid;
                    if(chi_constrained[j] > 1) chi_constrained[j] = 1;
                    if(chi_constrained[j] < 0) chi_constrained[j] = 0;
                    chi_sum += chi_constrained[j];
                }
                if ( (chi_sum - vol_constraint[i]*(double)natoms) > 0 ) l2 = lmid;
                else l1 = lmid;
            }
            for (int j=0; j<natoms; j++) {
                chi[j] = chi_constrained[j];
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
    fprintf(screen,"Chi type material ");
    for(int i=0; i<chi_map.properties.size(); i++) {
        fprintf(screen,"%s ",chi_map.properties[i].c_str());
    }
    fprintf(screen, "\n-------------------------\n");
    for (int i = 0; i < nchi; i++) {
        fprintf(screen,"%f %d %s ",chi_map.chis[i],chi_map.types[i],chi_map.material[i].c_str());
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
