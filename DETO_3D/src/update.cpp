#include <string.h>
#include "update.h"
#include "simulations.h"
#include "optimize.h"
#include "lammpsIO.h"
#include "universe.h"
#include "error.h"
#include "output.h"
#include <memory>


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
    std::istringstream lss(read_string);
    lss >> opt_type;
    string keyword;
    if(strcmp(opt_type.c_str(),"pertibation") == 0) {
        pop_size = 1;
        lss >> opt_par1 >> opt_par2;
        //opt_par1 = move_limit
    }
    else if(strcmp(opt_type.c_str(),"gradient_descent") == 0) {
        pop_size = 1;
        lss >> opt_par1 >> opt_par2;
    }
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
    else {
        err_msg = "ERROR: opt_type\"" + opt_type + "\" not recognised";
        error->errsimple(err_msg);
    }
    optimize->pop_size = pop_size;
    if(me == MASTER) {
        fprintf(screen,"\noptimization paramaters\n");
        fprintf(screen,"opt_type: %s\npop_size: %d\n",opt_type.c_str(),pop_size);
    }
}


// // ---------------------------------------------------------------
// // genetic algorythm for update
// void Update::genetic(const std::vector<std::vector<double>>& chi_pop,const std::vector<std::vector<int>>& mat_pop,const double* opt_obj_eval, const int* fitness)
// {
//     fprintf(screen,"\n");
//     chi_next.clear();
//     mat_next.clear();
//     //Selection
//     int selection_size = pop_size-gen_elitism;
//     std::vector<int> selection;
//     selection.clear();
//     int best;   
//     for(int i=0; i<selection_size; i++) {
//         best = rand() % pop_size;
//         for(int j=0; j<opt_par1-1; j++) {
//             int eval = rand() % pop_size;
//             if(opt_obj_eval[eval] < opt_obj_eval[best]) {
//                 best = eval;
//             }
//         }
//         chi_next.push_back(chi_pop[best]);
//         mat_next.push_back(mat_pop[best]);
//     }
//     //Crossover
//     for(int i=0; i<selection_size; i+=2) {
//         int j = i+1;
//         double temp = 0;
//         int tempint = 0;
//         // attempt to crossover at each gene individually
//         double cross_chance = ((double) rand() / (RAND_MAX));
//         if(cross_chance < opt_par2) {
//             for(int k=0; k<natoms; k++) {
//                 if(rand()%2 == 1) {
//                     temp = chi_next[i][k];
//                     chi_next[i][k] = chi_next[j][k];
//                     chi_next[j][k] = temp;
//                     tempint = mat_next[i][k];
//                     mat_next[i][k] = mat_next[j][k];
//                     mat_next[j][k] = tempint;
//                 }
//             }
//         }
//     }
//     //Mutation
//     for(int i=0; i<selection_size; i++) {
//         for(int j=0; j<natoms; j++) {
//             if(mat_pop[i][j] != -1) {
//                 double mutation_chance = ((double) rand() / (RAND_MAX));
//                 if(mutation_chance < opt_par3) {
//                     int mat_mut = (int)(rand() %optimize->nmat);
//                     int chi_mut = (int)(rand() %optimize->chi_map.nchi[mat_mut]);
//                     mat_next[i][j] = mat_mut;
//                     chi_next[i][j] = optimize->chi_map.chis[mat_mut][chi_mut];
//                 }
//             }
//         }
//     }
//     //add the elite solutions back in
//     for(int i=0; i<gen_elitism; i++) {
//         chi_next.push_back(chi_pop[fitness[i]]);
//         mat_next.push_back(mat_pop[fitness[i]]);
//     }
// }


// ---------------------------------------------------------------
// monte-carlo algorythm for update
void Update::monte_carlo(const std::vector<std::vector<double>>& chi_pop,const std::vector<std::vector<int>>& mat_pop,const double* opt_obj_eval)
{
    // double best_eval = opt_objective_eval[0];
    //     int best = 0;
    //     for(int i=0; i<optimize->pop_sizeps[0]; i++) {
    //         if(opt_objective_eval[i] < best_eval) {
    //             best = i;
    //             best_eval = opt_objective_eval[best];
    //         }
    //     }
    //     fprintf(screen,"\nBest choice is: %d value: %f\n",best, opt_objective_eval[best]);
    //     optimize->chi = optimize->chi_pop[best];
}

// ---------------------------------------------------------------
// monte-carlo algorythm for update
void Update::gradient_descent(const std::vector<std::vector<double>>& chi_pop,const std::vector<std::vector<int>>& mat_pop)
{
    //then should update chi_pop based on this sensitivity directly
    //create next chi vector;

    double* chi = optimize->chi_popps[0];
    int* mat = optimize->mat_popps[0];

    //should extract the sensitivities directly from a variable assigned in lammps
    //naively combine all sensitivities linearly we should warn the user about this or ofer them options to use more complex formula
    dchi = new double[natoms];
    for(int i=0; i<sims->n_sims; i++) {
        for(int j=0; j<sims->n_repeats[i]; j++) {
            for(int k=0; k<sims->sim_sens_names[i][j].size(); k++) {
                for (int l=0; l<natoms; l++) {
                    dchi[l] += sims->sim_sens_val[i][j][k][l];
                }
            }
        }
    }
    //normalise sensitivity between 0-1
    double min = 100;
    double max = 0;
    for (int i=0; i<natoms; i++) {
        if(mat[i]!=-1) {
            if(dchi[i] < min) min = dchi[i];
            if(dchi[i] > max) max = dchi[i];
        }
    }
    std::cout << "max " << max << " min " << min << "\n";
    double dchi_max = 0; double dchi_min = 1;
    for (int i=0; i<natoms; i++) {
        dchi[i] = (dchi[i]-min)/(max-min);
        if(mat[i]!=-1 && dchi[i] > dchi_max) dchi_max = dchi[i];
        if(mat[i]!=-1 && dchi[i] < dchi_min) dchi_min = dchi[i];
    }
    // for (int i=0; i<natoms; i++) std::cout << dchi[i] << "\n";
    std::cout << "max_sens " << dchi_max << " min_sens " << dchi_min << "\n";
    
    sensitivity_update(dchi,chi,mat);

    delete[] dchi;

    


    // this should be easy okay :)


}


// ---------------------------------------------------------------
// update using the pertibation method
void Update::perturbation(const vector<vector<double>>& chi_pop,const vector<vector<int>>& mat_pop,const double* opt_obj_eval)
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
    double* chi = optimize->chi_popps[0];
    int* mat = optimize->mat_popps[0];

    pert_sizeps.clear();
    pert_sizeps_cum.clear();
    for(int i=0; i<universe->nsc; i++) pert_sizeps.push_back(int(natoms/universe->nsc));
    for(int i=0; i<natoms%universe->nsc; i++) pert_sizeps[i]++;
    pert_sizeps_cum.push_back(0);
    for(int i=0; i<universe->nsc-1; i++) pert_sizeps_cum.push_back(pert_sizeps_cum[i]+pert_sizeps[i]);

    if(key==0) {
        dchips = new double[pert_sizeps[universe->color]];
        if(me==MASTER) dchi = new double[natoms];
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
            sims->run();
            if(key == 0) dchips[i] = -std::max(((obj_eval - optimize->evaluate_objective())/0.05),0.0); //brDchi needs to be repalced with dist to nearest chi
        } 
        else {
            std::string set_type = "set atom " + std::to_string(optimize->IDuns[pert_ID]) + " type " + std::to_string(optimize->chi_map.types[mat[pert_ID]][index-1]);
            lammpsIO->lammpsdo(set_type);
            lammpsIO->lammpsdo("delete_bonds all multi remove");
            for(int i=0; i<optimize->potentials.size(); i++) {
                lammpsIO->lammpsdo(optimize->potentials[i].c_str());
            }
            sims->run();
            if(key == 0) dchips[i] = -std::max(((optimize->evaluate_objective() - obj_eval)/0.05),0.0); //brDchi needs to be repalced with dist to nearest chi
        }

        if(me==universe->nsc) {
            fprintf(screen,"Completed: %d/%d perturbations\n",(i+1)*universe->nsc,natoms);
            fprintf(screen,"\x1b[A");
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

    sensitivity_update(dchi,chi,mat);

    if(key==0) {
        delete [] dchips;
        if(me==MASTER) delete[] dchi;
    }
}


// ---------------------------------------------------------------
// return an updated chi vector based on current vector and a sensitivity vector
void Update::sensitivity_update(const double* dchi, const double* chi, const int* mat)
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
    output->toplog("\n------------------\nnew dchi , chi, mat\n-------------------\n\n");
    for (int i=0; i<natoms; i++) {
        string logmsg = "";
        std::ostringstream ss;
        ss << dchi[i] << " " << chi_next[i] << " " << mat[i]; 
        logmsg = logmsg+ss.str(); ss.str(""); ss.clear();
        output->toplog(logmsg);
    }

}


// ---------------------------------------------------------------
// update chi population
void Update::update_chipop(vector<vector<double>>& chi_pop, vector<vector<int>>& mat_pop,const double* opt_obj_eval, const int* fitness)
{
    natoms = (int)lammpsIO->lmp->atom->natoms;
    chi_next = new double[natoms];
    //Direct towards user specified update method (some methods e.g genetic require only master, some require all subcomms)
    if(strcmp(opt_type.c_str(),"pertibation") == 0) {
        perturbation(chi_pop,mat_pop,opt_obj_eval);
    }
    if(me == MASTER){
        if(strcmp(opt_type.c_str(),"gradient_descent") == 0) {
            gradient_descent(chi_pop,mat_pop);
        }
        // if(strcmp(opt_type.c_str(),"genetic") == 0) {
        //     genetic(chi_pop,mat_pop,opt_obj_eval,fitness);
        // }
        // else if(strcmp(opt_type.c_str(),"monte-carlo") == 0) {
        //     monte_carlo(chi_pop,mat_pop,opt_obj_eval);
        // }
        
        double chi_sum = 0;
        int n_part_opt = 0;
        for (int i=0; i<natoms; i++) {
            // std::cout << chi_next[i] << "\n";
            chi_pop[0][i] = chi_next[i];
            if(mat_pop[0][i] != -1) {
                chi_sum += chi_pop[0][i];
                n_part_opt++;
            }
        }
        optimize->vol_frac = chi_sum/n_part_opt;
    }
    delete [] chi_next;
    MPI_Barrier(MPI_COMM_WORLD);
}


void Update::printall()
{
    fprintf(screen, "\n---------ALL ABOUT UPDATE----------\n");
}