#ifndef LAMMPSIO_H
#define LAMMPSIO_H

#include "pointers.h"
#include "mpi.h"
#include <string>
#include <sstream>
#include "lammps.h"         // these are LAMMPS include files
#include "input.h"
#include "atom.h"
#include "library.h"



#define MASTER 0

namespace DETO_NS {
    using std::vector;
    using std::string;
    
	class LammpsIO : protected Pointers {
	public:
		
		 LammpsIO(class DETO *);
		~LammpsIO();
		
        LAMMPS_NS::LAMMPS *lmp;
        
        void create();  //Create a new Lammps instance
        void printall();
        void lammpsdo(string); //passes a command to lammps
        double extract_natoms();
        void* extract_global(string);
        void* extract_atom_varaiable(string);  // spelling error
        void* extract_varaiable(string);
        int extract_setting(string);
        void* gather_atom_varaiable(string);
        void print_bonds();
        void set_type(int,int);
        // std::unique_ptr<int[]> extract_bonds();
        void extract_bonds(int* & bonds);
        
        string lmpThSt;    // string recording the themo style from input: used to add and evaluate tem computes when needed in source code of lammps
        
        /*
        string units;
        string atomstyle;
        
        double timestep;
        stringstream ss;
         */
        bool lammps_active; // makes the current processor if lammps has been already activated for its subcomm
        /*
        //string tdump_fname;  // temporary dump file name, used by maske to evaluate computes etc for output
        */
        bool wllog;
	private:
		//string fname;              // open inputcprs file
        int me;     // id of the current processor (rank)
		
	};
	
}

#endif
