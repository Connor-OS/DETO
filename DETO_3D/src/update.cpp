#include <string.h>
#include "update.h"
#include "simulations.h"
#include "optimize.h"
#include "lammpsIO.h"
#include "universe.h"
#include "error.h"
#include "output.h"
#include <memory>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>


using namespace DETO_NS;

Update::Update(DETO *deto) : Pointers(deto)
{

    // The inputcprs class is run by all processors in COMM_WORLD. This reads the id of the processor
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    key = universe->key;
}

// ---------------------------------------------------------------
// Class destructorgit
Update::~Update()
{
}

// ---------------------------------------------------------------
// function setting the properties of the optimization
void Update::set_opt_type(std::string read_string)
{
    /// Warning pop_size growth not supported
    //reduced functionality to deal with multi sims of differnet styles i.e opt_type vector
    // lost some flexability in the design here
    std::istringstream lss(read_string);
    lss >> opt_par1 >> opt_par2;
    string read_in;
    while(lss >> read_in)
        opt_type.push_back(read_in);
    pop_size = 1;

    // string keyword;
    // if(strcmp(opt_type.c_str(),"pertibation") == 0) {
    //     pop_size = 1;
    //     lss >> opt_par1 >> opt_par2;
    //     //opt_par1 = move_limit
    // }
    // else if(strcmp(opt_type.c_str(),"gradient_descent") == 0) {
    //     pop_size = 1;
    //     lss >> opt_par1 >> opt_par2;
    // }
    
    // else if(strcmp(opt_type.c_str(),"genetic") == 0) {
    //     gen_elitism = 0;
    //     lss >> pop_size >> opt_par1 >> opt_par2 >> opt_par3;
    //     // opt_par1 = tournement size cross, opt_par2 = crossover rate, opt_par3 = Mutation rate
    //     while(lss >> keyword) {
    //         if(std::strcmp(keyword.c_str(),"elitism") == 0) {
    //             lss >> gen_elitism;
    //         }
    //     }
    //     if(pop_size%2 == 1 || gen_elitism%2 == 1) {
    //         err_msg = "ERROR: please use even numbers only for population size and elitism";
    //         error->errsimple(err_msg);
    //     }
    // }
    // else if(strcmp(opt_type.c_str(),"monte-carlo") == 0) {
    //     lss >> pop_size >> opt_par1 >> opt_par2;
    // }
    // else {
    //     err_msg = "ERROR: opt_type\"" + opt_type + "\" not recognised";
    //     error->errsimple(err_msg);
    // }
    optimize->pop_size = pop_size;
    if(me == MASTER) {
        fprintf(screen,"\noptimization paramaters\n");
        for(int i=0; i<opt_type.size(); i++) fprintf(screen,"opt_type: %s\npop_size: %d\n",opt_type[i].c_str(),pop_size);
    }
}


// ---------------------------------------------------------------
// genetic algorythm for update
void Update::genetic(const std::vector<std::vector<double>>& chi_pop,const std::vector<std::vector<int>>& mat_pop,const double* opt_obj_eval, const int* fitness)
{
// Not Implemented
}


// ---------------------------------------------------------------
// monte-carlo algorythm for update
void Update::monte_carlo(const std::vector<std::vector<double>>& chi_pop,const std::vector<std::vector<int>>& mat_pop,const double* opt_obj_eval)
{
// Not Implemented
}


// ---------------------------------------------------------------
// monte-carlo algorythm for update
void Update::gradient_descent(const double* chi,const int* mat, vector<double>& dchi)
{

    LAMMPS_NS::LAMMPS *lmptr;
    lmptr = lammpsIO->lmp;
    int natoms = dchi.size();

    int tagintsize;
    int64_t nbonds;
    int* bonds;
    double* atom_pos;

    atom_pos = new double[3 * natoms];
    lammps_gather_atoms(lmptr, "x", 1, 3, atom_pos);

    tagintsize = lammps_extract_setting(lmptr, "tagint");
    
    if (tagintsize == 4)
        nbonds = *(int32_t *)lammps_extract_global(lmptr, "nbonds");
     else
        nbonds = *(int64_t *)lammps_extract_global(lmptr, "nbonds");
    bonds = new int[nbonds * 3];
    lammps_gather_bonds(lmptr, bonds);
    
    for (int i=0; i<nbonds*3; i++) std::cout << bonds[i] << std::endl;

    // Compute bond length
    std::vector<double> l(nbonds, 0.0);
    for (size_t k = 0; k < nbonds; ++k) {
        if (bonds[3 * k] == 0) {
            l[k] = 1.0;
        } else {
            size_t i = bonds[3 * k + 1]  - 1;
            size_t j = bonds[3 * k + 2]  - 1;
            l[k] = std::sqrt(std::pow(atom_pos[3 * i] - atom_pos[3 * j], 2) +
                             std::pow(atom_pos[3 * i + 1] - atom_pos[3 * j + 1], 2) +
                             std::pow(atom_pos[3 * i + 2] - atom_pos[3 * j + 2], 2));
        }
    }

    // Compute sensitivity
    double k_0 = 100.0;
    std::fill(dchi.begin(), dchi.end(), 0.0);

    for (size_t k = 0; k < nbonds; ++k) {
        size_t i = bonds[3 * k + 1] - 1;
        size_t j = bonds[3 * k + 2] - 1;
        dchi[i] += 0.5 * chi[i] * std::pow(chi[j], 2) * std::pow((1.0 - l[k]), 2);
        dchi[j] += 0.5 * chi[j] * std::pow(chi[i], 2) * std::pow((1.0 - l[k]), 2);
    }

    filter_sensitivities(chi,dchi,bonds,nbonds);

    delete[] bonds;
    delete[] atom_pos;

    //normalise sensitivity between 0-1
    double min = 100;
    double max = 0;
    for (int i=0; i<natoms; i++) {
        if(mat[i]!=-1) {
            if(dchi[i] < min) min = dchi[i];
            if(dchi[i] > max) max = dchi[i];
        }
    }

    double dchi_max = 0; double dchi_min = 1;
    for (int i=0; i<natoms; i++) {
        dchi[i] = (dchi[i]-min)/(max-min);
        if(mat[i]!=-1 && dchi[i] > dchi_max) dchi_max = dchi[i];
        if(mat[i]!=-1 && dchi[i] < dchi_min) dchi_min = dchi[i];
    }

    // for(int i=0; i< dchi.size(); i++) std::cout << dchi[i] << std::endl;
}


// ---------------------------------------------------------------
// update using the pertibation method
void Update::perturbation(const double* chi,const int* mat,const double* opt_obj_eval,vector<double>& dchi, int sim)
{   
    key = universe->key;
    double obj_eval;
    if(key == 0) {
        for(int i=1; i<universe->nsc; i++) {
            if(me == MASTER) {
                MPI_Send(&opt_obj_eval[0], 1, MPI_DOUBLE, universe->subMS[i], i, MPI_COMM_WORLD);
                obj_eval = opt_obj_eval[0];
            }
            if(universe->color == i) MPI_Recv(&obj_eval, 1, MPI_DOUBLE, MASTER, i, MPI_COMM_WORLD, &status);
        }
    }
    vector<int> pert_sizeps;
    vector<int> pert_sizeps_cum;

    pert_sizeps.clear();
    pert_sizeps_cum.clear();
    for(int i=0; i<universe->nsc; i++) pert_sizeps.push_back(int(natoms/universe->nsc));
    for(int i=0; i<natoms%universe->nsc; i++) pert_sizeps[i]++;
    pert_sizeps_cum.push_back(0);
    for(int i=0; i<universe->nsc-1; i++) pert_sizeps_cum.push_back(pert_sizeps_cum[i]+pert_sizeps[i]);

    if(key==0) {
        dchips = new double[pert_sizeps[universe->color]];
    }


    // perturbation loop
    for(int i=0; i<pert_sizeps[universe->color]; i++) {
        int pert_ID = pert_sizeps_cum[universe->color]+i;
        //  load inital config
        optimize->load_chi(0); // we may be able to improve the efficency of this. But is it worth it for the savings we will gain?
        // upate single chi and revaluate objective
        int index = optimize->chi_map.lookup(chi[pert_ID],mat[pert_ID]);
        if(index == -1) {
            if(key == 0)  dchips[i] = 0.0;
        }
        else if(index < optimize->chi_map.nchi[mat[pert_ID]]-1) {
            std::string set_type = "set atom " + std::to_string(optimize->IDuns[pert_ID]) + " type " + std::to_string(optimize->chi_map.types[mat[pert_ID]][index+1]);
            lammpsIO->lammpsdo(set_type);
            lammpsIO->lammpsdo("delete_bonds all multi remove");
            for(int i=0; i<optimize->potentials.size(); i++) {
                lammpsIO->lammpsdo(optimize->potentials[i].c_str());
            }
            sims->run_one(sim);
            // if(key == 0) dchips[i] = std::max(((obj_eval - optimize->evaluate_objective())/0.05),0.0); //brDchi needs to be repalced with dist to nearest chi
            if(key == 0) dchips[i] = (obj_eval - optimize->evaluate_objective())/0.05; //brDchi needs to be repalced with dist to nearest chi
        } 
        else {
            std::string set_type = "set atom " + std::to_string(optimize->IDuns[pert_ID]) + " type " + std::to_string(optimize->chi_map.types[mat[pert_ID]][index-1]);
            lammpsIO->lammpsdo(set_type);
            lammpsIO->lammpsdo("delete_bonds all multi remove");
            for(int i=0; i<optimize->potentials.size(); i++) {
                lammpsIO->lammpsdo(optimize->potentials[i].c_str());
            }
            sims->run_one(sim);
            if(key == 0) dchips[i] = (obj_eval - optimize->evaluate_objective())/0.05; //brDchi needs to be repalced with dist to nearest chi        }
        }
        if(me==universe->nsc) {
            fprintf(screen,"Completed: %d/%d perturbations\n",(i+1)*universe->nsc,natoms);
            // fprintf(screen,"\x1b[A");
        }
    }
    optimize->load_chi(0); // load inital config after perterbations complete

    if(me == MASTER) {
        for(int i=0; i<pert_sizeps[0]; i++) {
            dchi[i] = dchips[i];
        }
    }
    //communicate evals from submasters
    if(key == 0) {
        for(int i=1; i<universe->nsc; i++) {
            if(universe->color == i) MPI_Send(&dchips[0], pert_sizeps[i], MPI_DOUBLE, MASTER, i, MPI_COMM_WORLD);
            if(me == MASTER) MPI_Recv(&dchi[pert_sizeps_cum[i]], pert_sizeps[i], MPI_DOUBLE, universe->subMS[i], i, MPI_COMM_WORLD, &status);
        }
    }

    
    if(me == MASTER) {
        double min = 100;
        double max = -100;
        for (int i=0; i<natoms; i++) {
            if(mat[i]!=-1) {
                if(dchi[i] < min) min = dchi[i];
                if(dchi[i] > max) max = dchi[i];
            }
        }
        // std::cout << "max " << max << " min " << min << "\n";
        double dchi_max = 0; double dchi_min = 1;
        for (int i=0; i<natoms; i++) {
            dchi[i] = (dchi[i]-min)/(max-min);
            if(mat[i]!=-1 && dchi[i] > dchi_max) dchi_max = dchi[i];
            if(mat[i]!=-1 && dchi[i] < dchi_min) dchi_min = dchi[i];
        }
        // std::cout << "max_sens " << dchi_max << " min_sens " << dchi_min << "\n";
    }

    if(key==0) {
        delete [] dchips;
    }
}


void Update::filter_sensitivities(const double* chi, vector<double>& dchi, const int* bonds, const int nbonds)
{
    // Simple filtering using bondlist as the neighbor list, assuming rmin = 1.5 and l = 1
    //
    double rmin = 1.5;
    double l = 1.0;

    std::vector<double> dchi_filt(dchi.size(), 0.0);
    std::vector<double> tot(dchi.size(), rmin);

    for (size_t i = 0; i < dchi.size(); ++i) {
        dchi_filt[i] = dchi[i] * chi[i] * rmin;
    }


    for (size_t k = 0; k < nbonds; ++k) {
        size_t i = bonds[3 * k + 1] - 1;
        size_t j = bonds[3 * k + 2] - 1;
        double fac = rmin - l;
        tot[i] += fac;
        tot[j] += fac;
        dchi_filt[i] += dchi[j] * chi[j] * fac;
        dchi_filt[j] += dchi[i] * chi[i] * fac;
    }

    // Calculate and store the final result
    for (size_t i = 0; i < dchi.size(); ++i) {
        std::cout << chi[i] << " " << tot[i] << " " << dchi[i] << " " << dchi_filt[i] << std::endl;
        dchi[i] = dchi_filt[i] / (chi[i] * tot[i] + 0.001);
    }
}


// ---------------------------------------------------------------
// return an updated chi vector based on current vector and a sensitivity vector
void Update::sensitivity_update(const vector<double> dchi, const double* chi, const int* mat)
{
    // Updating chi while respecting constraints on volume fraction and range
    double l1 = 0.;
    double l2 = 1.;
    double lmid;
    double move = opt_par1;
    double vol_frac = opt_par2;
    double chi_sum;
    int n_part_opt;
    while (l2-l1 > 1e-15){
        n_part_opt = 0;
        chi_sum = 0.;
        lmid = 0.5*(l2+l1);
        for (int i=0; i<natoms; i++) {
            if(mat[i] == -1) {
                // chi_sum += 1;
                chi_next[i] = 2;
            }
            else {
                chi_next[i] = chi[i] * sqrt(dchi[i]/lmid);
                if (chi_next[i] > chi[i] + move) chi_next[i] = chi[i] + move;
                if (chi_next[i] > 1.) chi_next[i] = 1.;
                if (chi_next[i] < chi[i] - move) chi_next[i] = chi[i] - move;
                if (chi_next[i] < 0) chi_next[i] = 0;
                chi_sum += chi_next[i];
                n_part_opt++;
            }
        }        
        if ((chi_sum - vol_frac*(double)n_part_opt) > 0) l1 = lmid;
        else l2 = lmid;
    }
}


// ---------------------------------------------------------------
// update chi population
void Update::update_chipop(vector<vector<double>>& chi_pop, vector<vector<int>>& mat_pop,const double* opt_obj_eval, const int* fitness)
{
    natoms = (int)lammpsIO->lmp->atom->natoms;
    chi_next = new double[natoms];

    double* chi = optimize->chi_popps[0];
    int* mat = optimize->mat_popps[0];

    //initlaise an sensitivity vector of zeros for each sim
    vector<vector<double>> dchi = vector<vector<double>>();
    for(int i=0; i<sims->n_sims; i++) dchi.push_back(vector<double>(natoms, 0.0));

    for(int sim=0; sim<opt_type.size(); sim++) {
        // select method of computing sensitivity
        //Direct towards user specified update method (some methods e.g genetic require only master, some require all subcomms)
        if(strcmp(opt_type[sim].c_str(),"pertibation") == 0) {
            perturbation(chi,mat,opt_obj_eval,dchi[sim],sim);
        }
        if(me == MASTER){
            if(strcmp(opt_type[sim].c_str(),"gradient_descent") == 0) {
                gradient_descent(chi,mat,dchi[sim]);
            }
            // if(strcmp(opt_type.c_str(),"genetic") == 0) {
            //     genetic(chi_pop,mat_pop,opt_obj_eval,fitness);
            // }
            // else if(strcmp(opt_type.c_str(),"monte-carlo") == 0) {
            //     monte_carlo(chi_pop,mat_pop,opt_obj_eval);
            // }
        }
    }

    if(me == MASTER) {
        // average all sensitivity vectors
        vector<double> dchi_total = vector<double>(natoms, 0);
        for(int i=0; i<natoms; i++) {
            for(int j=0; j<dchi.size(); j++) {
                dchi_total[i] += dchi[j][i]/dchi.size();
                // std::cout << dchi[j][i] << " ";
            }
            // std::cout << dchi_total[i] << std::endl;
        } 
        sensitivity_update(dchi_total,chi,mat);

        double chi_sum = 0;
        int n_part_opt = 0;
        for (int i=0; i<natoms; i++) {
            chi_pop[0][i] = chi_next[i];
            if(mat_pop[0][i] != -1) {
                chi_sum += chi_pop[0][i];
                n_part_opt++;
            }
        }
        optimize->vol_frac = chi_sum/n_part_opt;
    }

    // delete [] dchi;
    delete [] chi_next;
    MPI_Barrier(MPI_COMM_WORLD);
}


void Update::printall()
{
    fprintf(screen, "\n---------ALL ABOUT UPDATE----------\n");
}