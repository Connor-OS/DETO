#include "optimize.h"
#include <sstream>
#include "error.h"
#include "lammpsIO.h"
#include "universe.h"
#include "output.h"
#include "simulations.h"
#include "update.h"
#include <iostream>
#include <algorithm>
#include <chrono>
// #include <memory>

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
void Optimize::read_chimap(string mapfname)
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
            getline(mapFile, read_string);
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
                            vector<string> check{"chi","material","type"};
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
                                getline(mapFile, read_string);
                                read_string = read_string.substr(0, read_string.find("#"));
                                if (!read_string.empty()) {
                                    std::istringstream lss(read_string);
                                    double tempchi,temp;
                                    string tempmat;
                                    int temptype;
                                    tot_nchi = 0;
                                    lss >> tempchi >> tempmat >> temptype;
                                    bool mat_found =false;
                                    for(int i=0; i<chi_map.material.size(); i++) {
                                        if(strcmp(chi_map.material[i].c_str(),tempmat.c_str()) == 0) {
                                            mat_found = true;
                                            tot_nchi += 1;
                                            chi_map.nchi[i] += 1;
                                            chi_map.chis[i].push_back(tempchi);
                                            chi_map.types[i].push_back(temptype);
                                            for(int prop=0; prop < chi_map.properties.size(); prop++) {
                                                lss >> temp;
                                                chi_map.values[i][prop].push_back(temp);
                                            }
                                        }
                                    }
                                    if(mat_found == false) {
                                        tot_nchi += 1;
                                        chi_map.nchi.push_back(1);
                                        chi_map.material.push_back(tempmat);
                                        chi_map.chis.push_back(vector<double>{tempchi});
                                        chi_map.types.push_back(vector<int>{temptype});
                                        chi_map.values.push_back(vector<vector<double>>());
                                        for(int prop=0; prop < chi_map.properties.size(); prop++) {
                                            lss >> temp;
                                            chi_map.values[chi_map.values.size()-1].push_back(vector<double>{temp});
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
void Optimize::add_constraint(string read_string)
{   
    string material,constraint_type,method;
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


// // ---------------------------------------------------------------
// // function setting the properties of the optimization
// void Optimize::set_opt_type(string read_string)
// {
//     std::istringstream lss(read_string);
//     lss >> opt_type;
//     if(strcmp(opt_type.c_str(),"sensitivity") == 0) {
//         pop_size = 1;
//         lss >> opt_par1 >> opt_par2;
//     }
//     else if(strcmp(opt_type.c_str(),"genetic") == 0) {
//         lss >> pop_size >> opt_par1 >> opt_par2;
//     }
//     if(me == MASTER) {
//         fprintf(screen,"\noptimization paramaters\n");
//         fprintf(screen,"opt_type: %s\npop_size: %d\n",opt_type.c_str(),pop_size);
//     }
// }

// ---------------------------------------------------------------
// function initializing values of chi from chi_map and types currently loaded into lammps instance
void Optimize::initialize_chi(vector<double>& chi,vector<int>& mat)
{
    if(IDuns) delete[] IDuns;
    lammpsIO->lammpsdo("compute tempID all property/atom id");
    lammpsIO->lammpsdo("compute tempType all property/atom type");
    lammpsIO->lammpsdo("dump tdID all custom 1 dump.temp_"+universe->SCnames[universe->color]+" c_tempID c_tempType x y z");
    lammpsIO->lammpsdo("run 1");     // a run1 in lammps to dump the temp and so prepare variables and computes
    lammpsIO->lammpsdo("undump tdID");

    MPI_Barrier(MPI_COMM_WORLD);
    
    // extracting unsorted atom ids
    natoms = (int)lammpsIO->lmp->atom->natoms; // total number of atoms in lammps, all groups all fixes)
    nlocal = (int)lammpsIO->lmp->atom->nlocal; // number of atoms in current processore, all types all fixes)
    
    aID = (double *)lammps_extract_compute(lammpsIO->lmp,(char *) "tempID",1,1);
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
    MPI_Bcast(&IDuns[0],natoms,MPI_INT,0,MPI_COMM_WORLD);

    //  initialise a chi vector on the master
    for(int i=0; i<natoms; i++) {
        bool type_found = false;
        for(int j=0; j<nmat; j++) {
            for(int k=0; k<chi_map.types[j].size(); k++) {
                if(typeuns[i] == chi_map.types[j][k]) {
                    chi.push_back(chi_map.chis[j][k]);
                    mat.push_back(j);
                    type_found = true;
                }
            }
        }
        if(type_found == false) {
            chi.push_back(2);    // todo: push_back(fabs(chi_map.chi_avg)+chi_map.chi_max+fabs(chi_map.chi_max)).... do this per-material use greatest chi_map.chi_avg and chi_map.chi_max between all materials
            mat.push_back(-1);
        }
    }
    //Enforce Volume constraint, or local volume constraint if any are defined here
    // chi = constrain_avg_chi(chi);
    // chi = constrain_local_avg_chi(chi);    //todo: to be implemented... not in a rush though

    MPI_Barrier(MPI_COMM_WORLD);
    delete[] IDpos;
    delete[] nID_each;
    delete[] typeuns;
    lammpsIO->lammpsdo("uncompute tempID");
    lammpsIO->lammpsdo("uncompute tempType");
    if(me != MASTER || dto->wplog == false) return;
    output->toplog("\n------------------\nInitial chi\n-------------------\n\n");
    for (int i=0; i<natoms; i++) {
        string logmsg = "";
        std::ostringstream ss;
        ss << IDuns[i] << " " << chi[i] << " " << mat[i]; 
        logmsg = logmsg+ss.str(); ss.str(""); ss.clear();
        output->toplog(logmsg);
    }
}


// ---------------------------------------------------------------
// function initializing population of chi from chi_map
void Optimize::initialize_chipop() 
{
    if(me == MASTER) {
        // To generate inital population mutate a population at a special inital rate
        // for(int i=0; i<pop_size; i++) {
        //     chi_pop.push_back(chi);
        //     mat_pop.push_back(mat);
        //     for(int j=0; j<natoms; j++) {
        //         if(mat[j] != -1 && i > 0) {
        //             double mutation_chance = ((double) rand() / (RAND_MAX));
        //             //(mat_index[j] + rand()%(nmat-1))%nmat  // is it typical to select the same material
        //             int mat_mut = (int)(rand() %nmat);
        //             int chi_mut = (int)(rand() %chi_map.nchi[mat_mut]);
        //             mat_pop[i][j] = mat_mut;
        //             chi_pop[i][j] = chi_map.chis[mat_mut][chi_mut];
        //         }
        //     }
        // }
        // for(int i=0; i<pop_size; i++) {
        //     constrain_avg_chi(i);
        //     constrain_local_avg_chi(chi_pop[i]);
        // }
        if(pop_size > 1) {
            // generate a population using a Latin Hypercube
            // vector of possible chi values
            vector<double> values;
            for(int i=0; i<chi_map.material.size(); i++) {
                for(int j=0; j<chi_map.chis[i].size(); j++) {
                    values.push_back(chi_map.chis[i][j]);
                }
            }
            vector<vector<double>> sample_space;
            int c = 0;
            while(chi_pop.size() < pop_size) {
                sample_space.clear();
                for(int i=0; i<natoms; i++) {
                    sample_space.push_back(values);
                }
                for(int i=0; i<values.size(); i++) {
                    chi_pop.push_back(vector<double>());
                    mat_pop.push_back(vector<int>());
                    for(int j=0; j<natoms; j++) {
                        if(mat[j] == -1) {
                            chi_pop[c].push_back(2);
                            mat_pop[c].push_back(-1);
                        }
                        else {
                            int sample = (int) rand() % (values.size()-i);
                            chi_pop[c].push_back(sample_space[j][sample]);
                            sample_space[j].erase(sample_space[j].begin()+sample);
                            mat_pop[c].push_back(0);
                        }
                    }
                    c++;
                }
            }
            
            while(chi_pop.size() > pop_size) {
                chi_pop.erase(chi_pop.end());
                mat_pop.erase(mat_pop.end());
            }
        }
        if(pop_size == 1) {
            chi_pop.push_back(chi);
            mat_pop.push_back(mat);
        }

        if(me == MASTER) fprintf(screen,"DONE Generating chi population\n\n");
        if(me != MASTER || !dto->wplog) return;
        output->toplog("\n------------------\nChi Population\n-------------------\n\n");
        for (int i=0; i<natoms; i++) {
            string logmsg = "";
            for(int j=0; j<pop_size; j++) {
                std::ostringstream ss;
                ss << "\t" << chi_pop[j][i] << " " << mat_pop[j][i] << "\t|"; 
                logmsg = logmsg+ss.str(); ss.str(""); ss.clear();
            }
            output->toplog(logmsg);
        }
    }
}


// ---------------------------------------------------------------
// function spliting chi_pop into subgroups and sending to individual subcomms
void Optimize::split_pop() 
{
    //code for deleting 2d array
    //Free each sub-array
    if(chi_popps) {
        for(int i=0; i<pop_sizeps[universe->color]; i++) {
            delete[] chi_popps[i]; 
            delete[] mat_popps[i]; 
        }
        //Free the array of pointers
        delete[] chi_popps; 
        delete[] mat_popps;
    }

    // assign sections of population to different subcommunicators, supports changing population size
    pop_sizeps.clear();
    pop_sizeps_cum.clear();
    for(int i=0; i<universe->nsc; i++) pop_sizeps.push_back(int(pop_size/universe->nsc));
    for(int i=0; i<pop_size%universe->nsc; i++) pop_sizeps[i]++; 
    pop_sizeps_cum.push_back(0);
    for(int i=0; i<universe->nsc-1; i++) pop_sizeps_cum.push_back(pop_sizeps_cum[i]+pop_sizeps[i]);

    //code for creating 2d array TODO: make this a function call that I can call anywhere in the code 
    int n1 = std::max(pop_sizeps[universe->color],1), n2 = natoms;

    chi_popps = new double*[n1];
    for(int i = 0; i < n1; i++){
        chi_popps[i] = new double[n2];
    }
    mat_popps = new int*[n1];
    for(int i = 0; i < n1; i++){
        mat_popps[i] = new int[n2];
    }

    //split chi population and send to submasters
    for(int i=0; i<universe->nsc; i++) {
        for(int j=0; j<pop_sizeps[i]; j++) {
            if(me == MASTER) {
                for(int k=0; k<universe->SCnp[i]; k++) {
                    if(universe->subMS[i]+k != 0) {
                        // fprintf(screen,"sending inside subcomm %i to member %i chi %i\n",i,k,j);
                        MPI_Send(&chi_pop[j+pop_sizeps_cum[i]][0], natoms, MPI_DOUBLE, universe->subMS[i]+k, j, MPI_COMM_WORLD);
                        MPI_Send(&mat_pop[j+pop_sizeps_cum[i]][0], natoms, MPI_INT, universe->subMS[i]+k, j+1, MPI_COMM_WORLD);
                    }
                }
            }
            if(universe->color == i) {
                if(me != MASTER) {
                    // fprintf(screen,"recieving inside subcomm %i to member %i chi %i\n",i,me-universe->subMS[i],j);
                    MPI_Recv(&chi_popps[j][0], natoms, MPI_DOUBLE, MASTER, j, MPI_COMM_WORLD, &status);
                    MPI_Recv(&mat_popps[j][0], natoms, MPI_INT, MASTER, j+1, MPI_COMM_WORLD, &status);
                }
                else {
                    for(int k=0; k<natoms; k++) {
                        chi_popps[j][k] = chi_pop[j][k];
                        mat_popps[j][k] = mat_pop[j][k];
                    }
                }
            }
        }
        // in case of population less than subcomms copy first chi_pop
        if(pop_sizeps[i] == 0) {
            if(me == MASTER) {
                for(int k=0; k<universe->SCnp[i]; k++) {
                    if(universe->subMS[i]+k != 0) {
                        // fprintf(screen,"sending inside subcomm %i to member %i chi %i\n",i,k,j);
                        MPI_Send(&chi_pop[0][0], natoms, MPI_DOUBLE, universe->subMS[i]+k, 0, MPI_COMM_WORLD);
                        MPI_Send(&mat_pop[0][0], natoms, MPI_INT, universe->subMS[i]+k, 1, MPI_COMM_WORLD);
                    }
                }
            }
            if(universe->color == i) {
                // fprintf(screen,"recieving inside subcomm %i to member %i chi %i\n",i,me-universe->subMS[i],j);
                MPI_Recv(&chi_popps[0][0], natoms, MPI_DOUBLE, MASTER, 0, MPI_COMM_WORLD, &status);
                MPI_Recv(&mat_popps[0][0], natoms, MPI_INT, MASTER, 1, MPI_COMM_WORLD, &status);
            }
        }
    }
    if(!dto->wplog) return;
    output->toplog("\n------------------\nChi sub-population\n-------------------\n\n");
    for (int i=0; i<natoms; i++) {
        string logmsg = "";
        std::ostringstream ss;
        ss << IDuns[i];
        for(int j=0; j<std::max(pop_sizeps[universe->color],1); j++) {
            if(mat_popps[j][i] != -1) {
                ss << "\t" << chi_popps[j][i] << " " << chi_map.material[mat_popps[j][i]] << "\t|"; 
                logmsg = logmsg+ss.str(); ss.str(""); ss.clear();
            }
            else{
                ss << "\tnon-opt\t|"; 
                logmsg = logmsg+ss.str(); ss.str(""); ss.clear();
            }
        }
        output->toplog(logmsg);
    }
}

// ---------------------------------------------------------------
// Function to constrain the total chi in the simulation, does not contain material
void Optimize::constrain_avg_chi(int id)
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
            if(strcmp(constraint_method[i].c_str(),"scale") == 0) {
                double l1 = 0.;  // chi_min - chi_map.chi_avg/10
                double l2 = 1000; //chi_map.chi_max[i]*(1/chi_map.chi_min[i]);   // chi_map.chi_max //maybe this constraint only works for positive chi
                double lmid; //different constraint types can allow us to accomidate more
                while (l2-l1 > 1e-8){
                    natoms_mat = 0;
                    chi_sum = 0.;
                    lmid = 0.5*(l2+l1);
                    for (int j=0; j<natoms; j++) {
                        if(mat_pop[id][j] == i) {
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
            else if(strcmp(constraint_method[i].c_str(),"shift") == 0) {
                double l1 = -chi_map.chi_max[i];
                double l2 = chi_map.chi_max[i];
                double lmid;
                while (l2-l1 > 1e-8){
                    natoms_mat = 0;
                    chi_sum = 0.;
                    lmid = 0.5*(l2+l1);
                    for (int j=0; j<natoms; j++) {
                        if(mat_pop[id][j] == i) {
                            chi_constrained[j] = chi_pop[id][j]+lmid;
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
// Applying the local volume constraint from (Aage et all 2017)
vector<double> Optimize::constrain_local_avg_chi(vector<double> chi)
{
    return chi;
}

// ---------------------------------------------------------------
// running the optimization
void Optimize::load_chi(int id)
// This function maps a continuous chi vector to it's closest discrete counterpart and then sets the lammps instance to be using this chi arrangment
// this should be called imediatly before any lammps run takes place (probably with the chi vector passeed as an argument)
{   
    if(sims->ndim == 2) lammpsIO->lammpsdo("read_dump dump.init_config 1 x y vx vy add keep box yes trim yes");
    else lammpsIO->lammpsdo("read_dump dump.init_config 1 x y z vx vy vz add keep box yes trim yes");
    //Set atom types
    for(int i=0; i<natoms; i++) {
        int k;
        for(int j=0; j<chi_map.material.size(); j++) {
            if(mat_popps[id][i] == j) {
                int l1 = 0;
                int l2 = chi_map.nchi[j]-1;
                while(l2-l1 > 1) {
                    k = int((l1+l2)/2);
                    if(chi_popps[id][i] < chi_map.chis[j][k]) l2 = k;
                    else l1 = k;
                }
                if((chi_map.chis[j][l1]+chi_map.chis[j][l2])/2 > chi_popps[id][i]) k = l1;
                else k = l2;
                string set_type = "set atom " + std::to_string(IDuns[i]) + " type " + std::to_string(chi_map.types[j][k]);
                lammpsIO->lammpsdo(set_type);
            }
        }
    }
    // Create pairs and bonds    
    lammpsIO->lammpsdo("delete_bonds all multi remove");
    for(int i=0; i<potentials.size(); i++) {
        lammpsIO->lammpsdo(potentials[i].c_str());
    }
}


// ---------------------------------------------------------------
// evaluate the objective function for the current chi 
double Optimize::evaluate_objective()
{
    for(int i=0; i<sims->n_sims; i++) {
        for(int j=0; j<sims->n_repeats[i]; j++) {
            for(int k=0; k<sims->sim_obj_names[i][j].size(); k++) {
                lammpsIO->lammpsdo("variable " + sims->sim_obj_names[i][j][k] + " equal " + std::to_string(sims->sim_obj_val[i][j][k])); //need to ad repeats probably by removing a layer here
            }
        }
    }
    lammpsIO->lammpsdo("variable obj equal " + obj_function);
    return *(double *)lammpsIO->extract_varaiable("obj");
}

// ---------------------------------------------------------------
// communicate the objective function evaluations back to the master and rank fitness 
void Optimize::communicate_objective(int* fitness)
{
    //move evals on MASTER
    if(me == MASTER) {
        for(int i=0; i<pop_sizeps[0]; i++) {
            opt_objective_eval[i] = opt_objective_evalps[i];
        }
    }
    //communicate evals from submasters
    if(key == 0) {
        for(int i=1; i<universe->nsc; i++) {
            if(universe->color == i) {
                MPI_Send(&opt_objective_evalps[0], pop_sizeps[i], MPI_DOUBLE, MASTER, i, MPI_COMM_WORLD);
            }
            if(me == MASTER) {
                MPI_Recv(&opt_objective_eval[pop_sizeps_cum[i]], pop_sizeps[i], MPI_DOUBLE, universe->subMS[i], i, MPI_COMM_WORLD, &status);
            }
        }
    }
    //communicate thermo info
    if(key == 0) {
        for(int i=1; i<universe->nsc; i++) {
            if(universe->color == i) {
                for(int j=0; j<pop_sizeps[i]; j++) {
                    MPI_Send(&thermo_val[pop_sizeps_cum[i]+j][0], thermo_string.size(), MPI_DOUBLE, MASTER, j, MPI_COMM_WORLD);
                }
            }
            for(int j=0; j<pop_sizeps[i]; j++) {
                if(me == MASTER) {
                    MPI_Recv(&thermo_val[pop_sizeps_cum[i]+j][0], thermo_string.size(), MPI_DOUBLE, universe->subMS[i], j, MPI_COMM_WORLD, &status);
                }
            }
        }
    }
    if(me == MASTER) {
        vector<double> opt_objective_eval_sorted (opt_objective_eval,opt_objective_eval+pop_size);
        sort(opt_objective_eval_sorted.begin(), opt_objective_eval_sorted.end());
        for(int i=0; i<pop_size; i++) {
            int j=0;
            while(opt_objective_eval_sorted[i] != opt_objective_eval[j]) j++;
            fitness[i] = j;
        }
    }
    MPI_Bcast(&fitness[0],pop_size,MPI_INT,0,MPI_COMM_WORLD);
}


// ---------------------------------------------------------------
// running the optimization
void Optimize::optrun()
{
    chi = vector<double>();
    mat = vector<int>();
    // initialize chi values from types and set other type-related quantities in the chi_map, then initalize a population of chi's based on sim_type
    initialize_chi(chi,mat);
    initialize_chipop();

    // store inital configuration to be reset back to
    if(universe->color == 0) lammpsIO->lammpsdo("write_dump all custom dump.init_config id x y z diameter type vx vy vz");
    for(int i=0; i<potentials.size(); i++) {
        lammpsIO->lammpsdo(potentials[i].c_str());
    }

    // begin optimization loop
    if(me == MASTER) fprintf(screen,"Starting Optimization\n");

    int n1 = pop_size, n2 = thermo_string.size();
    thermo_val = new double*[n1];
    for(int i = 0; i < n1; i++){
        thermo_val[i] = new double[n2];
    }

    int step = 0;
    while(step < 500) {
        if(me == MASTER) fprintf(screen,"Step %d\n\n",step);
        auto step_start = high_resolution_clock::now();
        // split population among avilable subcomms
        split_pop();
        //initalise container to hold individual evaluations of the objective function
        //remove previous instance
        if(opt_objective_eval) delete[] opt_objective_eval;
        if(opt_objective_evalps) delete[] opt_objective_evalps;
        //assign space in memory to hold objective function evaluations
        if(key == 0) {
            opt_objective_evalps = new double[pop_sizeps[universe->color]];
            if(me == MASTER) {
                opt_objective_eval = new double[pop_size];
            }
        }
        //run simulations on all population members
        for(int id=0; id<pop_sizeps[universe->color]; id++) {
            //reset simulation and load configuration
            load_chi(id);
            //run simulations and extract objectives
            sims->run();
            //extract thermo
            for(int i=0; i<thermo_string.size(); i++) {
                thermo_val[pop_sizeps_cum[universe->color]+id][i] = *(double *)lammpsIO->extract_varaiable(thermo_string[i].c_str());                
            }
            //evaluate combined objective function
            if(key == 0) {
                opt_objective_evalps[id] = evaluate_objective();
                if(me==universe->nsc) {
                    fprintf(screen,"Completed: %d/%d population members\n",(id+1)*universe->nsc,pop_size);
                    fprintf(screen,"\x1b[A");
                }
            }
        }

        //communicate all objective function evaluations back to the master and rank fitness
        int fitness[pop_size];
        communicate_objective(fitness);
        //update to next chi_pop
        update->update_chipop(chi_pop,mat_pop,opt_objective_eval,fitness);
        //write dumps
        if(me == MASTER) output->writethermo(step,opt_objective_eval,fitness);
        output->writedump(step,pop_size,fitness);
        step++;
        auto step_end = high_resolution_clock::now();
        auto step_duration = duration_cast<std::chrono::milliseconds>(step_end-step_start);      
        if(me == MASTER) fprintf(screen,"duration: %.3f seconds\n\n", step_duration.count()/1000.0);
    }

    //optimisation finished
    if(me == MASTER) fprintf(screen,"\nCOMMPLETED OPTIMIZATION RUN\n");

    if(opt_objective_eval) delete[] opt_objective_eval;
    if(opt_objective_evalps) delete[] opt_objective_evalps;
    if(chi_popps) {
        for(int i=0; i<pop_sizeps[universe->color]; i++) {
            delete[] chi_popps[i]; 
            delete[] mat_popps[i]; 
        }
        //Free the array of pointers
        delete[] chi_popps; 
        delete[] mat_popps;
    }
    delete[] IDuns;
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
    //Print population 
    fprintf(screen,"\nChi Population: %i\n",pop_size);
    fprintf(screen,"Per subcom:\n");
    for(int i=0; i<universe->nsc; i++) {
        fprintf(screen,"subcomm: %i: %i\n",i,pop_sizeps[i]);
    }
    
    //Print constraints
    fprintf(screen,"\nGloal volume constraints: \n");
    for(int i=0; i<nmat; i++) {
        fprintf(screen,"%s  f=%f \n",chi_map.material[i].c_str(),vol_constraint[i]);
    }
    fprintf(screen,"\n");
    fprintf(screen,"Local volume constraints: \n");
    for(int i=0; i<nmat; i++) {
        fprintf(screen,"%s  f=%f radius=%f \n",chi_map.material[i].c_str(),local_vol_constraint[i],local_vol_radius[i]);
    }

    fprintf(screen,"Objective function is %s\n",obj_function.c_str());

    fprintf(screen, "\n---------------------------------------\n\n");
}
