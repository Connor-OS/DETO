#include <string.h>
#include "update.h"
#include "simulations.h"
#include "optimize.h"
#include "universe.h"
#include "error.h"


using namespace DETO_NS;

Update::Update(DETO *deto) : Pointers(deto)
{

    // The inputcprs class is run by all processors in COMM_WORLD. This reads the id of the processor
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    key = universe->key;

}

// ---------------------------------------------------------------
// Class destructor
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
    if(strcmp(opt_type.c_str(),"sensitivity") == 0) {
        pop_size = 1;
        lss >> opt_par1 >> opt_par2;
    }
    else if(strcmp(opt_type.c_str(),"genetic") == 0) {
        lss >> pop_size >> opt_style >> opt_par1 >> opt_par2;
    }
    else if(strcmp(opt_type.c_str(),"monte-carlo") == 0) {
        lss >> pop_size >> opt_par1 >> opt_par2;
    }
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


// ---------------------------------------------------------------
// evaluate the objective function for the current chi 
void Update::evaluate_objective(int id)
{
    //not needed here
}


// ---------------------------------------------------------------
// genetic algorythm for update
void Update::genetic(const std::vector<std::vector<double>>& chi_pop,const std::vector<std::vector<int>>& mat_pop,const double* opt_obj_eval)
{
    // fprintf(screen,"obj_eval\n");
    // for(int j=0; j<pop_size; j++) {
    //     fprintf(screen,"%.3f ",opt_obj_eval[j]);
    // }
    fprintf(screen,"\n");
    // fprintf(screen,"Starting tournement selection\n" );
    chi_next.clear();
    //tournement selection
    std::vector<int> selection;
    selection.clear();
    int best;
    int tour_size = 3;    
    for(int i=0; i<2*pop_size; i++) {
        best = rand() % pop_size;
        for(int j=0; j<tour_size-1; j++) {
            int eval = rand() % pop_size;
            if(opt_obj_eval[eval] < opt_obj_eval[best]) {
                best = eval;
            }
        }
        selection.push_back(best);
    }
    // fprintf(screen,"Winning selection: " );
    // for(int i=0; i<2*pop_size; i++) fprintf(screen,"%d ",selection[i]);
    // fprintf(screen,"\n" );

    //crossover
    std::vector<double> chi;
    std::vector<int> mat;
    for(int i=0; i<pop_size; i++) {
        for(int j=0; j<natoms; j++) {
            int inherit = 2*i+rand()%2;
            chi.push_back(chi_pop[selection[inherit]][j]);
            mat.push_back(mat_pop[selection[inherit]][j]);
        }
        chi_next.push_back(chi);
        mat_next.push_back(mat);
        chi.clear();
        mat.clear();
        // for(int k=0; k<100; k++) {
        //     fprintf(screen,"p1:%.0f p2:%.0f ch:%.0f\n",chi_pop[selection[2*i]][k],chi_pop[selection[2*i+1]][k],chi_next[i][k]);
        // }
        // fprintf(screen,"next crossover:\n");
    }

    // fprintf(screen,"Chi_next_pop\n");
    // for(int i=0; i<natoms; i++) {
    //     for(int j=0; j<pop_size; j++) {
    //         fprintf(screen,"%.2f ",chi_next[j][i]);
    //     }
    //     fprintf(screen,"\n");
    // }

    //mutation
    double mutation_rate = 0.01;
    for(int i=0; i<pop_size; i++) {
        for(int j=0; j<natoms; j++) {
            if(mat_pop[i][j] != -1) {
                double mutation_chance = ((double) rand() / (RAND_MAX));
                //(mat_index[j] + rand()%(nmat-1))%nmat  // is it typical to select the same material
                if(mutation_chance < mutation_rate) {
                    int mat_mut = (int)(rand() %optimize->nmat);
                    int chi_mut = (int)(rand() %optimize->chi_map.nchi[mat_mut]);
                    mat_next[i][j] = mat_mut;
                    chi_next[i][j] = optimize->chi_map.chis[mat_mut][chi_mut];
                }
            }
        }
    }

}


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
// genetic algorythm for update
void Update::sensitivity()
{
    
}


// ---------------------------------------------------------------
// update chi population
std::vector<std::vector<double>> Update::update_chipop(const std::vector<std::vector<double>>& chi_pop,const std::vector<std::vector<int>>& mat_pop,const double* opt_obj_eval)
{
    // fprintf(screen,"Inside update Object\n");
    // pop_size = chi_pop.size();
    natoms = chi_pop[0].size();
    // fprintf(screen,"pop_size is: %d natoms is: %d\n",pop_size,natoms);
    // fprintf(screen,"Chi_pop\n");
    // for(int i=0; i<natoms; i++) {
    //     for(int j=0; j<pop_size; j++) {
    //         fprintf(screen,"%.2f ",chi_pop[j][i]);
    //     }
    //     fprintf(screen,"\n");
    // }
    // fprintf(screen,"obj_eval\n");
    // for(int j=0; j<pop_size; j++) {
    //     fprintf(screen,"%.3f ",opt_obj_eval[j]);
    // }
    // fprintf(screen,"\n");
    
    if(strcmp(opt_type.c_str(),"genetic") == 0) {
        genetic(chi_pop,mat_pop,opt_obj_eval);
    }
    else if(strcmp(opt_type.c_str(),"monte-carlo") == 0) {
        monte_carlo(chi_pop,mat_pop,opt_obj_eval);
    }
    else if(strcmp(opt_type.c_str(),"sensitivity") == 0) {
        sensitivity();
    }

    // for(int i=0; i<natoms; i++) {
    //     for(int j=0; j<pop_size; j++) {
    //         fprintf(screen,"%.2f ",chi_next[j][i]);
    //     }
    //     fprintf(screen,"\n");
    // }

    return chi_next;
}


void Update::printall()
{
    fprintf(screen, "\n---------ALL ABOUT UPDATE----------\n");
}