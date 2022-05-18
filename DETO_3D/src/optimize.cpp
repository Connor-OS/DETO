#include "optimize.h"
#include <sstream>
#include "error.h"
#include "lammpsIO.h"
#include "universe.h"
#include "output.h"
#include "simulations.h"
#include <iostream>

//#include "chemistry.h"
//#include "store.h"
//#include "error.h"
/*#include <stdlib.h>
#include <stdio.h>
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
// add a constraint per material that can later be enforced on the design
void Optimize::add_constraint(std::string read_string)
{   
    std::string material,constraint_type,method;
    double constraint, radius;
    std::istringstream lss(read_string);
    lss >> material >> constraint_type >> method >> constraint;
    if(constraint_type == "avgchi") {
        for(int i=0; i<nmat; i++) {
            if(chi_map.material[i] == material) {
                vol_constraintYN[i] = true;
                vol_constraint[i] = constraint;
                constraint_method[i] = method;
            }
        }
    }
    else if(constraint_type == "local_avgchi") {
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
// function setting the properties of the optimization
void Optimize::set_opt_type(std::string read_string)
{
    std::istringstream lss(read_string);
    lss >> opt_type;
    if(strcmp(opt_type.c_str(),"sensitivity") == 0) {
        pop_size = 1;
        lss >> opt_par1 >> opt_par2;
    }
    else if(strcmp(opt_type.c_str(),"genetic") == 0) {
        lss >> pop_size >> opt_par1 >> opt_par2;
    }
    if(me == MASTER) {
        fprintf(screen,"\noptimization paramaters\n");
        fprintf(screen,"opt_type: %s\npop_size: %d\n",opt_type.c_str(),pop_size);
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
    
    IDuns = new int[natoms];
    typeuns = new int[natoms];
    
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

    //  initialise a chi vector on the master
    if(me == MASTER) {
        for(int i=0; i<natoms; i++) {
            bool type_found = false;
            for(int j=0; j<nmat; j++) {
                for(int k=0; k<chi_map.types[j].size(); k++) {
                    if(typeuns[i] == chi_map.types[j][k]) {
                        chi.push_back(chi_map.chis[j][k]);
                        mat.push_back(chi_map.material[j]);
                        mat_index.push_back(j);
                        type_found = true;
                    }
                }
            }
            if(type_found == false) {
                chi.push_back(2);    // todo: push_back(fabs(chi_map.chi_avg)+chi_map.chi_max+fabs(chi_map.chi_max)).... do this per-material use greatest chi_map.chi_avg and chi_map.chi_max between all materials
                mat.push_back("non-opt");
                mat_index.push_back(-1);
            }
        }
        //Enforce Volume constraint, or local volume constraint if any are defined here
        // chi = constrain_avg_chi(chi);
        // chi = constrain_local_avg_chi(chi);    //todo: to be implemented... not in a rush though
    }
    
    //Broadcast chi vector to each processor
    chi_each = new double[natoms];
    mat_each = new int[natoms];
    if(me == MASTER) {
        chi_each = &chi[0];
        mat_each = &mat_index[0];
    }

    MPI_Bcast(&chi_each[0],natoms,MPI_DOUBLE,0,MPI_COMM_WORLD); // Most likely remove this infavour of only broadcasting what each subcom needs
    MPI_Bcast(&mat_each[0],natoms,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&IDuns[0],natoms,MPI_INT,0,MPI_COMM_WORLD);

    if(me > 0) {
        for(int i=0; i<natoms; i++) {
            chi.push_back(chi_each[i]);
            mat_index.push_back(mat_each[i]);
            if(mat_each[i] > -1) mat.push_back(chi_map.material[mat_each[i]]);
            else mat.push_back("non-opt");
        }
    }

    // output->toplog("\n------------------\nInitial chi\n-------------------\n\n");
    // for (int i=0; i<natoms; i++) {
    //     std::string logmsg = "";
    //     std::ostringstream ss;
    //     ss << IDuns[i] << " " << chi[i] << " " << mat[i]; 
    //     logmsg = logmsg+ss.str(); ss.str(""); ss.clear();
    //     output->toplog(logmsg);
    // }
}


// ---------------------------------------------------------------
// function initializing population of chi from chi_map
void Optimize::initialize_chipop() 
{
    if(me == MASTER) {
        // To generate inital population mutate a population at a special inital rate
        double inital_mutation_rate = 0.1;
        for(int i=0; i<pop_size; i++) {
            chi_pop.push_back(chi);
            mat_pop.push_back(mat);
            mat_index_pop.push_back(mat_index);
            for(int j=0; j<natoms; j++) {
                if(strcmp(mat[j].c_str(),"non-opt") != 0 && i > 0) {
                    double mutation_chance = ((double) rand() / (RAND_MAX));
                    //(mat_index[j] + rand()%(nmat-1))%nmat  // is it typical to select the same material
                    if(mutation_chance < inital_mutation_rate) {
                        int mat_mut = (int)(rand() %nmat);
                        int chi_mut = (int)(rand() %chi_map.nchi[mat_mut]);
                        mat_index_pop[i][j] = mat_mut;
                        mat_pop[i][j] = chi_map.material[mat_mut];
                        chi_pop[i][j] = chi_map.chis[mat_mut][chi_mut];
                    }
                }
            }
        }

        for(int i=0; i<pop_size; i++) {
            constrain_avg_chi(i);
            constrain_local_avg_chi(chi_pop[i]);
        }

        // output->toplog("\n------------------\nChi Population\n-------------------\n\n");
        // for (int i=0; i<natoms; i++) {
        //     std::string logmsg = "";
        //     for(int j=0; j<pop_size; j++) {
        //         std::ostringstream ss;
        //         ss << "\t" << chi_pop[j][i] << " " << mat_pop[j][i] << "\t|"; 
        //         logmsg = logmsg+ss.str(); ss.str(""); ss.clear();
        //     }
        //     output->toplog(logmsg);
        // }
    }
}


// ---------------------------------------------------------------
// function initializing population of chi from chi_map
void Optimize::split_pop() 
{
    //split chi population and send to submasters
    pop_sizeps = 0;
    for(int i=0; i<pop_size; i++) {
        if(me == MASTER) {
            chi_each = &chi_pop[i][0];
            mat_each = &mat_index_pop[i][0];
            if(universe->subMS[i%(universe->nsc)] != MASTER) {
                MPI_Send(&chi_each[0], natoms, MPI_DOUBLE, universe->subMS[i%(universe->nsc)], 0, MPI_COMM_WORLD);
                MPI_Send(&mat_each[0], natoms, MPI_INT, universe->subMS[i%(universe->nsc)], 1, MPI_COMM_WORLD);
            }
        }
        if(me == universe->subMS[i%(universe->nsc)]) {
            if(universe->subMS[i%(universe->nsc)] != MASTER) {
                MPI_Recv(&chi_each[0], natoms, MPI_DOUBLE, MASTER, 0, MPI_COMM_WORLD, &status);
                MPI_Recv(&mat_each[0], natoms, MPI_INT, MASTER, 1, MPI_COMM_WORLD, &status);
            }
            pop_sizeps++;
            chi_popps.push_back(std::vector<double>());
            mat_index_popps.push_back(std::vector<int>());
            mat_popps.push_back(std::vector<std::string>());
            for(int j=0; j<natoms; j++) {
                chi_popps[chi_popps.size()-1].push_back(chi_each[j]);
                mat_index_popps[mat_index_popps.size()-1].push_back(mat_each[j]);
                if(mat_each[j] > -1) mat_popps[mat_popps.size()-1].push_back(chi_map.material[mat_each[j]]);
                else mat_popps[mat_popps.size()-1].push_back("non-opt");
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    // communicate sub_pop size to slaves processors
    if(key == 0) {
        for(int j=1; j<universe->SCnp[universe->color]; j++) {
            MPI_Send(&pop_sizeps, 1, MPI_INT, universe->subMS[universe->color]+j, 1, MPI_COMM_WORLD);
        }
    }
    if(key > 0) MPI_Recv(&pop_sizeps, 1, MPI_INT, universe->subMS[universe->color], 1, MPI_COMM_WORLD, &status);
    //Submasters send to individual slaves
    for(int i=0; i<pop_sizeps; i++) {
        if(key == 0) {
            chi_each = &chi_popps[i][0];
            mat_each = &mat_index_popps[i][0];
            for(int j=1; j<universe->SCnp[universe->color]; j++) {
                MPI_Send(&chi_each[0], natoms, MPI_DOUBLE, universe->subMS[universe->color]+j, i, MPI_COMM_WORLD);
                MPI_Send(&mat_each[0], natoms, MPI_INT, universe->subMS[universe->color]+j, i+pop_sizeps, MPI_COMM_WORLD);
            }
        }
        if(key > 0) {
            MPI_Recv(&chi_each[0], natoms, MPI_DOUBLE, universe->subMS[universe->color], i, MPI_COMM_WORLD, &status);
            MPI_Recv(&mat_each[0], natoms, MPI_INT, universe->subMS[universe->color], i+pop_sizeps, MPI_COMM_WORLD, &status);
            chi_popps.push_back(std::vector<double>());
            mat_popps.push_back(std::vector<std::string>());
            mat_index_popps.push_back(std::vector<int>());
            for(int j=0; j<natoms; j++) {
                chi_popps[chi_popps.size()-1].push_back(chi_each[j]);
                mat_index_popps[mat_index_popps.size()-1].push_back(mat_each[j]);
                if(mat_each[j] > -1) mat_popps[mat_popps.size()-1].push_back(chi_map.material[mat_each[j]]);
                else mat_popps[mat_popps.size()-1].push_back("non-opt");
            }
        }
    }
    output->toplog("\n------------------\nChi sub-population\n-------------------\n\n");
    for (int i=0; i<natoms; i++) {
        std::string logmsg = "";
        std::ostringstream ss;
        ss << IDuns[i];
        for(int j=0; j<pop_sizeps; j++) {
            ss << "\t" << chi_popps[j][i] << " " << mat_popps[j][i] << "\t|"; 
            logmsg = logmsg+ss.str(); ss.str(""); ss.clear();
        }
        output->toplog(logmsg);
    }
}

// ---------------------------------------------------------------
// Function to constrain the total chi in the simulation, does not contain material
void Optimize::constrain_avg_chi(int id) // This function needs to have different modes for different volume constraints
//It is best not to use this function with negative chi as it is untested
// To do add better warning messages telling the user if this constraint is poorly used.
// Distinguish local chi
{
    double chi_constrained[natoms];
    for(int i=0; i<natoms; i++) {
        chi_constrained[i] = chi_pop[id][i];
    }
    for(int i=0; i<nmat; i++) {
        double chi_sum = 0.;
        int natoms_mat = 0;
        if(vol_constraintYN[i] == true){
            if(std::strcmp(constraint_method[i].c_str(),"scale") == 0) {
                double l1 = 0.;  // chi_min - chi_map.chi_avg/10
                double l2 = 1000; //chi_map.chi_max[i]*(1/chi_map.chi_min[i]);   // chi_map.chi_max //maybe this constraint only works for positive chi
                double lmid; //different constraint types can allow us to accomidate more
                while (l2-l1 > 1e-8){
                    natoms_mat = 0;
                    chi_sum = 0.;
                    lmid = 0.5*(l2+l1);
                    for (int j=0; j<natoms; j++) {
                        if(mat_pop[id][j] == chi_map.material[i]) {
                            chi_constrained[j] = chi_pop[id][j]*lmid;
                            if(chi_constrained[j] > chi_map.chi_max[i]) chi_constrained[j] = chi_map.chi_max[i];
                            if(chi_constrained[j] < chi_map.chi_min[i]) chi_constrained[j] = chi_map.chi_min[i];
                            chi_sum += chi_constrained[j];
                            natoms_mat++;
                        }
                    }
                    if ((chi_sum - vol_constraint[i]*(double)natoms_mat) > 0) l2 = lmid;
                    else l1 = lmid;
                }
            }
            else if(std::strcmp(constraint_method[i].c_str(),"shift") == 0) {
                double l1 = -chi_map.chi_max[i];
                double l2 = chi_map.chi_max[i];
                double lmid;
                while (l2-l1 > 1e-8){
                    natoms_mat = 0;
                    chi_sum = 0.;
                    lmid = 0.5*(l2+l1);
                    for (int j=0; j<natoms; j++) {
                        if(mat[j] == chi_map.material[i]) {
                            chi_constrained[j] = chi[j]+lmid;
                            if(chi_constrained[j] > chi_map.chi_max[i]) chi_constrained[j] = chi_map.chi_max[i];
                            if(chi_constrained[j] < chi_map.chi_min[i]) chi_constrained[j] = chi_map.chi_min[i];
                            chi_sum += chi_constrained[j];
                            natoms_mat++;
                        }
                    }
                    if ((chi_sum - vol_constraint[i]*(double)natoms_mat) > 0) l2 = lmid;
                    else l1 = lmid;
                }
            }
            else {
                err_msg = "ERROR: Unrecognised constraint method";
                error->errsimple(err_msg); // return an error flag because inside submaster
                error_flag = true;
            }
            if ((chi_sum - vol_constraint[i]*(double)natoms_mat) > 0.1 || (chi_sum - vol_constraint[i]*(double)natoms_mat) < -0.1) {
                err_msg = "ERROR: was not able to apply constraint sucessfully to %s",chi_map.material[i].c_str();
                error->errsimple(err_msg);
            }
            for (int j=0; j<natoms; j++) {
                chi_pop[id][j] = chi_constrained[j];
            }
        }
    }
}


// ---------------------------------------------------------------
// running the optimization
std::vector<double> Optimize::constrain_local_avg_chi(std::vector<double> chi)
{
    return chi;
}

// ---------------------------------------------------------------
// running the optimization
void Optimize::load_chi(int id)
// This function maps a continuous chi vector to it's closest discrete counterpart and then sets the lammps instance to be using this chi arrangment
// this should be called imediatly before any lammps run takes place (probably with the chi vector passeed as an argument)
{   
    //Set atom types
    for(int i=0; i<natoms; i++) {
        int k = 0;
        for(int j=0; j<chi_map.material.size(); j++) {
            if(std::strcmp(mat_popps[id][i].c_str(),chi_map.material[j].c_str()) == 0) {
                int l1 = 0;
                int l2 = chi_map.nchi[j]-1;
                while(l2-l1 > 1) {
                    int k = (l1+l2)/2;
                    if(chi_popps[id][i] < chi_map.chis[j][k]) l2 = k;
                    else l1 = k;
                }
                if((chi_map.chis[j][l1]+chi_map.chis[j][l2])/2 > chi_popps[id][i]) k = l1;
                else k = l2;
                std::string set_type = "set atom " + std::to_string(IDuns[i]) + " type " + std::to_string(chi_map.types[j][k]);
                lammpsIO->lammpsdo(set_type);
            }
        }
    }
    lammpsIO->lammpsdo("run 0");
    // Create pairs and bonds
    for(int i=0; i<chi_map.types.size(); i++) {
        for(int j=0; j<chi_map.types[i].size(); j++) lammpsIO->lammpsdo("group " + std::to_string(chi_map.types[i][j]) + " type " + std::to_string(chi_map.types[i][j]));
    }
    
    lammpsIO->lammpsdo("delete_bonds all multi remove");
    for(int i=0; i<potentials.size(); i++) {
        lammpsIO->lammpsdo(potentials[i].c_str());
    }
    for(int i=0; i<chi_map.types.size(); i++) {
        for(int j=0; j<chi_map.types[i].size(); j++) lammpsIO->lammpsdo("group " + std::to_string(chi_map.types[i][j]) + " delete");
    }
}


// ---------------------------------------------------------------
// evaluate the objective function for the current chi 
void Optimize::evaluate_objective()
{
    for(int i=0; i<sims->sim_attributes.size(); i++) {
        for(int j=0; j<sims->n_repeats[i]; j++) {
            for(int k=0; k<sims->sim_obj_names[i][j].size(); k++) {
                fprintf(screen,"Objective %s simulation: %d repeat: %d value: %f\n",sims->sim_obj_LMPnames[i][j][k].c_str(),i,j,sims->sim_obj_val[i][j][k]);
            }
        }
    }
    // Need to create a new vector opt_obj_eval to be filled in by submasters then comunicated back to MASTER for design update

    // use a run 0 to evaluate in LAMMPS through all subprocessors
}


// ---------------------------------------------------------------
// running the optimization
void Optimize::optrun()
{
    // initialize chi values from types and set other type-related quantities in the chi_map file too
    
    initialize_chi();
    
    // for certain optimization types, e.g. GA, we may have to create a first set of chi_vectors to then initite the while loop
    if(strcmp(opt_type.c_str(),"genetic") == 0) {
        initialize_chipop();
        split_pop();
        if(me == MASTER) fprintf(screen,"DONE Generating chi population\n\n");
    }

    MPI_Barrier(MPI_COMM_WORLD);
    
    for(int i=0; i<pop_sizeps; i++) {
        if(key == 0) fprintf(screen,"Running sim %d.%d\n",universe->color,i+1);
        load_chi(i);
        // lammpsIO->lammpsdo("run 10");
        lammpsIO->lammpsdo("write_data data.init."+std::to_string(universe->color)+std::to_string(i));
        sims->run();
        // evaluate_objective(); // each chi configuration should have a filled out obj_val vector by this point as all simulations will have been run
    }
    
    MPI_Barrier(MPI_COMM_WORLD);

    if(me == MASTER) fprintf(screen,"\nCOMMPLETED SIM RUN\n");
    
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
