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
// Extract number of atoms in simulation
double LammpsIO::extract_natoms()
{   
    double natoms = lammps_get_natoms(lmp);
    return natoms;
}


// ---------------------------------------------------------------
// Extract per atom variable of atoms in simulation
void* LammpsIO::extract_atom_varaiable(std::string toextract)
{   
    void* atom_properties = lammps_extract_atom(lmp,toextract.c_str());
    return atom_properties;
}


// ---------------------------------------------------------------
// Extract per atom variable of atoms in simulation
void* LammpsIO::extract_varaiable(std::string toextract)
{   
    void* variable = lammps_extract_variable(lmp,toextract.c_str(),NULL);
    return variable;
}


// ---------------------------------------------------------------
// Extract per atom variable of atoms in simulation
void* LammpsIO::extract_global(std::string toextract)
{   
    void* variable = lammps_extract_global(lmp,toextract.c_str());
    return variable;
}


// ---------------------------------------------------------------
// read and print bond information directly from lammps
void LammpsIO::print_bonds()
{   
    LAMMPS_NS::bigint num_bonds;
    int** bonds;
    num_bonds = lmp->atom->nbonds;
    bonds = lmp->atom->bond_type;
    fprintf(screen,"%ld\n",num_bonds);
    for(int i=0; i<num_bonds; i++) {
        fprintf(screen,"%d\n",*bonds[i]);
        *lmp->atom->bond_type[i] = 66; //we can modify the bond data directly from here
    }
}


// ---------------------------------------------------------------
// read and print bond information directly from lammps
void LammpsIO::set_type(int id,int type)
{
    lmp->atom->type[id] = type; //don't do this unless you are an expert
}

// ---------------------------------------------------------------
// Extract per atom variable of atoms in simulation
void* LammpsIO::gather_atom_varaiable(char * toextract)
{   
    int64_t natoms;
    int tagintsize = lammps_extract_setting(lmp, "tagint");
    if (tagintsize == 4)
        natoms = *(int32_t *)lammps_extract_global(lmp, "natoms");
    else
        natoms = *(int64_t *)lammps_extract_global(lmp, "natoms");
    void *atom_properties;
    atom_properties = malloc(natoms*2*tagintsize);
    lammps_gather_concat(lmp,toextract,1,1, atom_properties);

    void * properties = atom_properties;
    // free(atom_properties);
    
    return properties;
}


// ---------------------------------------------------------------
// Printing info about the inputcprs class (possibly useful for debugging)
void LammpsIO::printall()
{
	fprintf(screen,"\n---------ALL ABOUT LAMMPSIO----------\n");
	//fprintf(screen,"inputcprs filename =  %s\n",fname.c_str());
	fprintf(screen,"---------------------------------------\n\n");
	
}
