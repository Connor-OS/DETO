#include "deto.h"

/*
#include "memory.h"
 */
#include "error.h"
#include "inputdeto.h"
#include "universe.h"
#include "output.h"
#include "lammpsIO.h"
 /*

#include "chemistry.h"
#include "solution.h"
#include "fix.h"
#include "fix_delete.h"
#include "fix_Cfoo.h"
#include "relax.h"
#include "setconc.h"
#include "krun.h"
#include "randm.h"
#include "fix_nucleate.h"
#include "store.h"

#ifdef MASKE_WITH_NUFEB
#include "fix_nufeb.h"
#endif
#ifdef MASKE_WITH_SPECIATION
#include "spec.h"
#endif
 */

using namespace DETO_NS;

// ---------------------------------------------------------------
// Initialize class
DETO::DETO(int narg, char **arg)
{
    screen = NULL;
    screen=stdout;
    wplog = false;
    /*
    plog = NULL;
    
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
    

    
    /*
    step = 0;
    tempo = 0;
    doublestep = (double)step;
    kB = 1.;
    hpl = 1.;
    
    Rtypes.clear();
    Ttypes.clear();
	
	memory = new Memory(this);
    */
    error = new Error(this);
    inputdeto = new Inputdeto(this,narg,arg);
    universe = new Universe(this);
    output = new Output(this);
    lammpsIO = new LammpsIO(this);
     /*
    
    chem = new Chemistry(this);
    solution = new Solution(this);
    fix = new Fix(this);
    fix_del = new Fix_delete(this);
#ifdef MASKE_WITH_NUFEB
    fix_nufeb = new Fix_nufeb(this);
#endif
    krun = new Krun(this);
    randm = new Randm(this);
    fix_cfoo = new Fix_Cfoo(this);
    relax = new Relax(this);
    setconc = new Setconc(this);
#ifdef MASKE_WITH_SPECIATION
    spec = new Spec(this);
#endif
    fix_nucl = new Fix_nucleate(this);
    store = new Store(this);
     */
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
    
    if (me==MASTER) fprintf(screen,"Deleting output class\n");
    MPI_Barrier(MPI_COMM_WORLD);
    delete output;
    
    
    if (me==MASTER) {
        fprintf(screen,"Deleting lammpsIO class\n");
        lammpsIO->printall();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    delete lammpsIO;
    
    /*
    if (me==MASTER) {
        fprintf(screen,"Deleting chemistry class\n");
        chem->printall();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    delete chem;
    if (me==MASTER) {
        fprintf(screen,"Deleting solution class\n");
        solution->printall();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    delete solution;
    if (me==MASTER) {
        fprintf(screen,"Deleting fix class\n");
        fix->printall();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    delete fix;
    if (me==MASTER) {
        fprintf(screen,"Deleting krun class\n");
        krun->printall();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    delete krun;
    if (me==MASTER) {
        fprintf(screen,"Deleting fix_delete class\n");
        fix_del->printall();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    delete fix_del;
    if (me==MASTER) {
        fprintf(screen,"Deleting randm class\n");
        randm->printall();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    delete randm;
    if (me==MASTER) {
        fprintf(screen,"Deleting Cfoo class\n");
        fix_cfoo->printall();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    delete fix_cfoo;
    if (me==MASTER) {
        fprintf(screen,"Deleting relax class\n");
        relax->printall();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    delete relax;
#ifdef MASKE_WITH_NUFEB
    delete fix_nufeb;
    if (me==MASTER) {
        fprintf(screen,"Deleting fix_nufeb class\n");
        //relax->printall();
    }
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    if (me==MASTER) {
        fprintf(screen,"Deleting fix_nucleate class\n");
        fix_nucl->printall();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    delete fix_nucl;
    if (me==MASTER) {
        fprintf(screen,"Deleting store class\n");
        store->printall();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    delete store;
    
     */
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
