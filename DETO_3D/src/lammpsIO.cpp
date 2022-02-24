#include "lammpsIO.h"
#include "universe.h"
//#include "error.h"
/*#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>

#include "clinker.h"
#include "pyramids.h"
#include "hydration.h"
#include "memory.h"
#include "csh.h"
#include "c2s.h"
#include "c3s.h"
#include "water.h"
#include <string.h>
*/

using namespace DETO_NS;


LammpsIO::LammpsIO(DETO *deto) : Pointers(deto)
{
    // The inputcprs class is run by all processors in COMM_WORLD. This reads the id of the processor
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    
    if (me == MASTER) fprintf(screen,"Generating lammpsIO class\n");

    lammps_active=false;
    lmpThSt = "thermo_style custom step atoms";
    /*
    units = "real";
    atomstyle = "ellipsoid";
    timestep = 1.;
     */
}



// ---------------------------------------------------------------
// Class destructor
LammpsIO::~LammpsIO()
{
    delete lmp;
}



// ----------------------------------------------------------------
// Create a new Lammps instance
void LammpsIO::create()
{
    // if lammps has not been already activated by this proc (just for safety) and if the processor pertains to a subcomm for which lammps was turend on in the input file..
    if (!lammps_active) {
        
        // option to print or not print lammps output to screen (useful for debugging but littering MASKE's output otherwsie)
        bool noscreen = true;
        int argn = 0;
        char **words;
        if (noscreen) {
            words = (char**)calloc (3,sizeof(char*));
            words[0] = "placeholder";
            words[1] = "-screen";
            words[2] = "none";
            argn = 3;
        }
        else words = NULL;

        lmp = new LAMMPS_NS::LAMMPS(argn,words,universe->subcomm);
        //lmp = new LAMMPS_NS::LAMMPS(0,NULL,universe->subcomm);
        lammps_active=true;
        
        delete [] words;
       
        
        fprintf(screen,"Proc[%d]: Lammps instance created, as part of subcomm %s\n",me,universe->SCnames[universe->color].c_str());
        
        // creates subcomm-specific lammps logfile (and temporary dumpfile name, currently turned off)
        std::string todo;
        todo = "log log."+universe->SCnames[universe->color]+".lammps";
        lammpsIO->lammpsdo(todo);

        
        
    }
}




// ---------------------------------------------------------------
// Printing info about the inputcprs class (possibly useful for debugging)
void LammpsIO::lammpsdo(std::string todo)
{
   lmp->input->one(todo.c_str());
    
}






// ---------------------------------------------------------------
// Printing info about the inputcprs class (possibly useful for debugging)
void LammpsIO::printall()
{
	fprintf(screen,"\n---------ALL ABOUT LAMMPSIO----------\n");
	//fprintf(screen,"inputcprs filename =  %s\n",fname.c_str());
	fprintf(screen,"---------------------------------------\n\n");
	
}
