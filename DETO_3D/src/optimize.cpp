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
    bool found_nchi = false;
    bool found_nmat = false;

    if (!mapFile.is_open()) {
        err_msg = "ERROR: cannot read file \"" + mapfname + "\"";
        error->errsimple(err_msg);
    }
    else {
        MPI_Barrier(MPI_COMM_WORLD);
        // READ FILE (all processors need to know this)
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
                            vol_constraint.push_back(-1.0);
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
                            while (lss >> word) {
                                chi_map.properties.push_back(word);
                                chi_map.values.push_back(std::vector<double>());
                                for(int i=0; i<nchi; i++) {
                                    chi_map.values[chi_map.values.size()-1].push_back(0);
                                }
                                if (word == "chi") {
                                    chi_set = true;
                                }
                                if (word == "material") {
                                    material_set = true;
                                }
                                if (word == "type") {
                                    type_set = true;
                                }
                            }
                            if (material_set == false || type_set == false || chi_set == false) {
                                err_msg = "ERROR: Please specify chi, material, and type in chi map";
                                error->errsimple(err_msg);
                            }
                            // Populate chi_map
                            for (int i = 0; i < nchi; i++) {
                                std::getline(mapFile, read_string);
                                if (!read_string.empty()) {
                                    std::istringstream lss(read_string);
                                    for (int j=0; j<chi_map.properties.size(); j++) {
                                        // double value;
                                        double tempval;
                                        lss >> tempval;
                                        chi_map.values[j][i] = tempval;
                                    }
                                }
                                else {
                                    err_msg = "ERROR: Fewer than specified chi's given";
                                    error->errsimple(err_msg);
                                }
                            }
                            // Move chi and type into separate vectors
                            int c_i,t_i;
                            for(int i=0; i<chi_map.properties.size(); i++) {
                                if(chi_map.properties[i] == "chi") {
                                    c_i =  i;
                                    for(int j=0; j<nchi; j++) {
                                        chi_map.chis.push_back(chi_map.values[i][j]);
                                    }
                                }
                                if(chi_map.properties[i] == "type") {
                                    t_i = i;
                                    for(int j=0; j<nchi; j++) {
                                        chi_map.types.push_back((int)chi_map.values[i][j]);
                                    }
                                }
                            }
                            chi_map.values.erase(chi_map.values.begin()+c_i,chi_map.values.begin()+t_i);
                            chi_map.properties.erase(chi_map.properties.begin()+c_i,chi_map.properties.begin()+t_i); 
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
    if(constraint_type == "volume") {
        vol_constraint[mat_ID-1] = constraint;
    }
    else if(constraint_type == "local_volume") {
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
    nploc =universe->SCnp[universe->color];
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
            IDpos[i] =IDpos[i-1]+nID_each[i-1];
        }
    }
    IDuns = new int[natoms];
    if (key>0) {
        int dest = 0;
        MPI_Send(&tID[0], nlocal, MPI_INT, dest, 1, (universe->subcomm));
    }
    if (key==0) {
        for (int i=0; i<nID_each[0]; i++) {
            IDuns[i] = tID[i];
        }
        for (int source=1; source<nploc; source++) {
            MPI_Recv(&IDuns[IDpos[source]], nID_each[source], MPI_INT, source, 1, (universe->subcomm), &status);
        }
    }

    
    /*
    
    
    if (me != MASTER && universe-->color==0) {
        MPI_Send(&nlocal, 1, MPI_INT, MASTER, 0, universe->subcomm);
        MPI_Send(&localID, nlocal, MPI_INT, MASTER, 1, universe->subcomm);
        MPI_Send(&localtype, nlocal, MPI_INT, MASTER, 2, universe->subcomm);
    }
    else if (me == MASTER) {
        int atom_index = 0;
        for(int i=0; i<nlocal; i++) {
            atomID[atom_index] = localID[atom_index];
            atomtype[atom_index] = localtype[atom_index];
            atom_index++;
        }
        for(int i=0; i<2; i++) {
            int temp;
            MPI_Recv(&temp, 1, MPI_INT, i+1, 0, MPI_COMM_WORLD ,MPI_STATUS_IGNORE);
            int tempID[temp],temptype[temp];
            MPI_Recv(&tempID, temp, MPI_INT, i+1, 1, MPI_COMM_WORLD ,MPI_STATUS_IGNORE);
            MPI_Recv(&temptype, temp, MPI_INT, i+1, 2, MPI_COMM_WORLD ,MPI_STATUS_IGNORE);
            for(int j=0; j<temp; j++) {
                atomID[atom_index] = tempID[j];
                atomtype[atom_index] = temptype[j];
                atom_index++;
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    

    MPI_Bcast(&atomID, natoms, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&atomtype, natoms, MPI_INT, MASTER, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
*/
    
    // extract vector of lammps types
    // aID = (int*)lammpsIO->gather_atom_varaiable("id");
    // atype = (int*)lammpsIO->gather_atom_varaiable("type");

    sleep(me);
    for(int i=0; i<natoms; i++) {
        //fprintf(screen,"PROC %d ID: %d type %d\n",me,atomID[i],atomtype[i]);
        if (key==0){fprintf(screen,"PROC %d ID: %d\n",me,IDuns[i]);} //,atomtype[i]);
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
    fprintf(screen,"Chi type ");
    for(int i=0; i<chi_map.properties.size(); i++) {
        fprintf(screen,"%s ",chi_map.properties[i].c_str());
    }
    fprintf(screen, "\n-------------------------\n");
    for (int i = 0; i < nchi; i++) {
        fprintf(screen,"%f %d ",chi_map.chis[i],chi_map.types[i]);
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
