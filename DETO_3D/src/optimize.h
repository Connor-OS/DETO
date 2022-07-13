#ifndef OPTIMIZE_H
#define OPTIMIZE_H

#include "pointers.h"
#include "universe.h"
// #include "simulations.h"
#include <string>
#include <vector>
#include <map>
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
	
	class Optimize : protected Pointers {
	public:
        Optimize(class DETO *);
		~Optimize();
            
        int me;     // id of the current processor (rank)

        void read_chimap(std::string);
        void add_constraint(std::string);
        void set_opt_type(std::string);
        void optrun();
        void printall();

        std::vector<std::string> potentials; // vector containing all information regarding the potentials available for the optimization
        int pop_size; // size of population of solutions

        struct Chi_map // structure containing the information mapping chi to type and other properties
        {
            std::vector<std::string> material;         // vector containing names of "material" encountered in chi map

            std::vector<std::vector<double>> chis;       // vector containing "chi" values in chi map file column per material
            std::vector<std::vector<int>> types;         // vector containing "types" values in chi map file column per material
            
            std::vector<std::string> properties;    // vector containing names of "other properties" in chi map file
            std::vector<std::vector<std::vector<double>>> values;   // matrix containing values of "other properties" in chi map file per material
        
            std::vector<double> chi_max; // max value of chi specified in chi_map for each material
            std::vector<double> chi_min; // min value of chi specified in chi_map for each material
            std::vector<double> chi_avg; // average value of chi specified in chi_map for each material
            std::vector<int> nchi;   //vector containing the number of chi values per material in the optimization
        };
        struct Chi_map chi_map; // Instance of chi_map to be used throughout optimization run to convert chi to type
        int nmat;   //number of materials in the optimization

	private:

        std::vector<bool> flag_avgchi_cstr;   // vector (one per material) of flags, true if user has defined a volume constraint on average chi through the system for that material

        //Constraint and design variables
        std::vector<bool> vol_constraintYN; // vector specifying if a volume constraint is set for each material    
        std::vector<double> vol_constraint; // volume constraint number between 0 and 1 per material
        std::vector<bool> local_vol_constraintYN; // vector specifying if a local volume constraint is set for each material      
        std::vector<double> local_vol_constraint; // local volume constraint number between 0 and 1 per material
        std::vector<double> local_vol_radius;  // Radius over which local material volume is constrained per material
        std::vector<std::string> constraint_method; // type of application of constraint per material, can be either scale or shift
        std::vector<std::string> local_constraint_method; // type of application of local constraint per material, can be either scale or shift
        std::string err_msg, read_string, word;
        int tot_nchi; // total number of chi in chi map all materials
        int natoms; // number of atoms in LAMMPS
        int nlocal; // number of atoms in current proc from LAMMPS

        std::vector<double> chi;   // chi values per atom (those with type not in chi_map will be assigned chi > chi_max
        std::vector<int> mat;   // numerical index asociated to the material

        std::vector<std::vector<double>> chi_pop;   // vector of vector of chi values for the whole population
        std::vector<std::vector<int>> mat_pop;  // vector of vector of material values for the whole population

        std::vector<int> pop_sizeps; // pop_size per subcommunicator
        std::vector<int> pop_sizeps_cum; // pop_size per subcommunicator

        double ** chi_popps;   // holds chi values per subcomm
        int ** mat_popps;  //holds material index per subcomm

        //Optimisation variables
        double* opt_objective_eval;
        double* opt_objective_evalps;
        
        int init_ts; //timestep after initalisation
        double tol; // optimisation tolerance
        //
    
        double* aID;  // ID of all atoms in LAMMPS in current proc
        double* atype;  // type of all atoms in LAMMPS in current proc
        int* IDuns; // unsortd IDs of all atoms in LAMMPS
        int* typeuns; //unsorted types of all atoms in LAMMPS
        int* nID_each;  //array with number of IDs in each processor in current subcomm
        // double* chi_each; // array containing chi in each processor
        // int* mat_each; // array containing the material in each processor
        int* ID;

        int nploc;  // number of processors in local subcomm
        int* IDpos;     // position of local tID array in submaster's unsorted list of IDs        
        int key;   // rank of current proc in current subcomm

        bool error_flag = false;

		//std::string fname;              // open inputcprs file
        //int me;     // id of the current processor (rank)
		
        MPI_Status status;

        void initialize_chi();
        void initialize_chipop();
        void update_chipop();
        void split_pop();
        void constrain_avg_chi(int id); 
        std::vector<double> constrain_local_avg_chi(std::vector<double>); //todo: change to void
        void load_chi(int);
        void evaluate_objective(int id);
	};
	
}

#endif
