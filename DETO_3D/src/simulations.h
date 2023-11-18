#ifndef SIMULATIONS_H
#define SIMULATIONS_H

#include "pointers.h"
#include "universe.h"
#include <string>
#include <vector>

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
        std::vector<std::string> cstgs_type;    // type of change (binary or increment) for cstgs sims
        std::vector<std::string> cstgs_incty;    // type of "increment" search: can be linear, power, or log
        std::vector<double> cstgs_par1,cstgs_par2, cstgs_par3, cstgs_par4;    // parameters guiding update of var in cstgs sims.
        // If  "increment linear" then par1 = start value, par2 = step, par3 = max number increments
        // if "increment power", then par1 = start value, par2 = step, par3 = exponent, par4 = max number increments
        // if "increment log", then par1 = start value, par2 = step, par3 = max number increments
        // if "binary" , then par1 = left value, par2 = right value, par3 = max number increments, par 4 = tolerance
       
        std::vector<std::string> cstgs_crit; //criterion to escape a cstgs sim
        std::vector<std::vector<std::string>> cstgs_crit_vnms; //name of lammps variables used in the criterion string
        
        int ndim;   //dimensionality of the simulation e.g 2D/3D
        int n_sims; //number of simulations defined in the input script
        std::vector<int> n_repeats; //number of repeats
        std::vector<std::string> sim_repeat_file; //vector holding path to repeat file data
        
        std::vector<std::vector<std::string>> sim_attributes; //vector sim attributes
        
        std::vector<std::vector<std::vector<std::string>>> sim_obj_names; // Unique ID of each objective
        std::vector<std::vector<std::vector<std::string>>> sim_obj_LMPnames; // LAMMPS varaible names linked to objective
        std::vector<std::vector<std::vector<double>>> sim_obj_val; //one list of objectives values for each simulation in a repeat for each user-defined simulation in the input file  (NB: if repeat = no, then the second-level vectors will simply have only one entry, which is a vector of objectives names and value for that simulation)

        std::vector<std::vector<std::vector<std::string>>> sim_sens_names; // Unique ID of each sensitivity
        std::vector<std::vector<std::vector<std::string>>> sim_sens_LMPnames; // LAMMPS varaible names linked to sensitivity
        std::vector<std::vector<std::vector<double*>>> sim_sens_val; //one list of sensitivity values for each simulation in a repeat for each user-defined simulation in the input file  (NB: if repeat = no, then the second-level vectors will simply have only one entry, which is a vector of sensitivity names and value for that simulation)

        std::vector<std::vector<std::string>> sim_rep_vars;      // vector of lammps-like variables in repeat file: one vector for each simulation
        std::vector<std::vector<std::vector<double>>> sim_rep_val;  // vector of values for each variable in repeat file, one for each simulation


        void printall();
        void add(std::string);
        
        void add_attribute(std::string);
        void add_objective(std::string);
        void add_sensitivity(std::string);
        void add_constraint(std::string);

        void read_repeat(std::string);

        void run();
        void run_one(int sim);
        
	private:
        std::string err_msg, read_string, word, read_string2;
        	
	};
	
}

#endif
