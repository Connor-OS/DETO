#include "deto.h"
#include "error.h"
#include "inputdeto.h"
#include "universe.h"
#include "output.h"
#include "lammpsIO.h"
#include "optimize.h"
#include "simulations.h"
#include "update.h"


using namespace DETO_NS;

// ---------------------------------------------------------------
// Initialize class
DETO::DETO(int narg, char **arg)
{
    screen = NULL;
    screen=stdout;
    wplog = false;
    
    plog = NULL;
    /*
    nulog_flag = false;
    */
    
    MPI_Comm_rank(MPI_COMM_WORLD, &me);

    
    if (me==MASTER) {
        fprintf(screen,"\n---------STARTING SIMULATION----------\n");
        fprintf(screen,"argc = %d\n",narg);
        for (int i=0; i<narg; i++) {
            fprintf(screen,"arg[%d] = %s\n",i,arg[i]);
        }
        fprintf(screen,"---------------------------------------\n\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    
    error = new Error(this);
    inputdeto = new Inputdeto(this,narg,arg);
    universe = new Universe(this);
    output = new Output(this);
    lammpsIO = new LammpsIO(this);
    optimize = new Optimize(this);
    sims = new Simulations(this);
    update = new Update(this);
}

// ---------------------------------------------------------------
// Class destructor
DETO::~DETO()
{
    if (me==MASTER) fprintf(screen,"Deleting deto class\n");
    MPI_Barrier(MPI_COMM_WORLD);
    
    
     if (me==MASTER) {
        fprintf(screen,"Deleting inputdeto class\n");
        inputdeto->printall();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    delete inputdeto;
    
    
    if (me==MASTER) fprintf(screen,"Deleting error class\n");
    MPI_Barrier(MPI_COMM_WORLD);
    delete error;
    
    if (me==MASTER) {
        fprintf(screen,"Deleting output class\n");
        output->printall();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    delete output;
    
    
    if (me==MASTER) {
        fprintf(screen,"Deleting lammpsIO class\n");
        lammpsIO->printall();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    delete lammpsIO;
    
    
    if (me==MASTER) {
        fprintf(screen,"Deleting optimize class\n");
        optimize->printall();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    delete optimize;
    
    if (me==MASTER) {
        fprintf(screen,"Deleting simulations class\n");
        sims->printall();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    delete sims;
    
    if (me==MASTER) {
        fprintf(screen,"Deleting update class\n");
        update->printall();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    delete update;

    if (me==MASTER) {
        fprintf(screen,"Deleting universe class\n");
        universe->printall();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Comm_free(&(universe->subcomm));
    delete universe;
    
}



// ---------------------------------------------------------------
// Print stuff
void DETO::printall()
{
    if (me==MASTER) {
        fprintf(screen,"\n---------END OF SIMULATION----------\n");
        //fprintf(screen,"water to cement ratio =  %f\n",wc);
        fprintf(screen,"---------------------------------------\n\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);
}
