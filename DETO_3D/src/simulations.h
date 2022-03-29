#ifndef SIMULATIONS_H
#define SIMULATIONS_H

#include "pointers.h"
#include <string>
#include <vector>
//#include "mpi.h"
//#include <string>
//#include <vector>
/*#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
*/

#define MASTER 0

namespace DETO_NS {
	
	class Simulations : protected Pointers {
	public:
        Simulations(class DETO *);
		~Simulations();
    
        int me;     // id of the current processor (rank)
        std::vector<std::string> sim_names; //vector holding sim names
        
        std::vector<std::string> sim_types; //vector holding sim types (run or cstgs)
        std::vector<std::string> cstgs_varname;    // name of variables to be changed during a cstgs simulation
        std::vector<std::string> cstgs_type;    // type of change (binary_chop or lin_increment) for cstgs sims
        std::vector<double> cstgs_par1,cstgs_par2;    // parameters guiding update of var in cstgs sims
        std::vector<double> cstgs_tol;    // tolerance for cgtgs of binary_chop type
        std::vector<std::string> cstgs_crit; //criterion to escape a cstgs sim
        
        std::vector<int> n_repeats; //number of repeats
        std::vector<std::string> sim_repeat_file; //vector holding path to repeat file data

        
        std::vector<std::vector<std::string>> sim_attributes; //vector sim attributes
        
        std::vector<std::vector<std::vector<std::string>>> sim_obj_names; // Unique ID of each attribute
        std::vector<std::vector<std::vector<std::string>>> sim_obj_LMPnames; // LAMMPS varaible names linked to objective
        std::vector<std::vector<std::vector<double>>> sim_obj_val; //one list of objectives names and values for each simulation in a repat for each user-defined simulation in the input file  (NB: if repeat = no, then the second-level vectors will simply have only one entry, which is a vector of objectives names and value for that simulation)

        std::vector<std::vector<std::string>> sim_rep_vars; // name of LAMMPS variables for repeat
        std::vector<std::vector<std::vector<double>>> sim_rep_val; // values to be repeated at

        void printall();
        void add(std::string);
        
        void add_attribute(std::string);
        void add_objective(std::string);

        void read_repeat(std::string);
        
	private:
        std::string err_msg, read_string, word;
        
        
        //void check_name_unique(std::string);

		//std::string fname;              // open inputcprs file
        //int me;     // id of the current processor (rank)
		
	};
	
}

#endif
