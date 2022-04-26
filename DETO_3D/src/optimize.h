#ifndef OPTIMIZE_H
#define OPTIMIZE_H

#include "pointers.h"
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
        void initialize_chi();
        void constrain_vol();
        void constrain_local_vol();
        void optrun();
        void printall();
        
	private:
        struct Chi_map // structure containing the information mapping chi to type and other properties
        {
            std::vector<std::string> material;         // vector containing names of "material" encountered in chi map

            std::vector<std::vector<double>> chis;       // vector containing "chi" values in chi map file column per material
            std::vector<std::vector<int>> types;         // vector containing "types" values in chi map file column per material
            
            std::vector<std::string> properties;    // vector containing names of "other properties" in chi map file
            std::vector<std::vector<std::vector<double>>> values;   // matrix containing values of "other properties" in chi map file per material
        };
        struct Chi_map chi_map; // Instance of chi_map
        struct Chi_map chi_map_sorted; // Instance of sorted chi_map per material

        std::vector<bool> flag_avgchi_cstr;   // vector (one per material) of flags, true if user has defined a volume constraint on average chi through the system for that material

        std::vector<bool> vol_constraintYN; // vector specifying if a volume constraint is set for each material    
        std::vector<double> vol_constraint; // volume constraint number between 0 and 1 per material
        std::vector<bool> local_vol_constraintYN; // vector specifying if a local volume constraint is set for each material      
        std::vector<double> local_vol_constraint; // local volume constraint number between 0 and 1 per material
        std::vector<double> local_vol_radius;  // Radius over which local material volume is constrained per material
        std::string err_msg, read_string, word;
        std::vector<double> chi_max; // max value of chi specified in chi_map for each material
        std::vector<double> chi_avg; // average value of chi specified in chi_map for each material
        std::vector<int> nchi;   //vector containing the number of chi values per material in the optimization
        int tot_nchi; // total number of chi in chi map all materials
        int nmat;   //number of materials in the optimization

        int natoms; // number of atoms in LAMMPS
        int nlocal; // number of atoms in current proc from LAMMPS
        std::vector<double> chi;   // chi values per atom (those with type not in chi_map will be assigne chi > chi_max

        double* aID;  // ID of all atoms in LAMMPS in current proc
        double* atype;  // type of all atoms in LAMMPS in current proc
        int* IDuns; // unsortd IDs of all atoms in LAMMPS
        int* typeuns; //unsorted types of all atoms in LAMMPS
        int* nID_each;  //array with number of IDs in each processor in current subcomm
        
        int nploc;  // number of processors in local subcomm
        
        int* IDpos;     // position of local tID array in submaster's unsorted list of IDs
        
        int key;   // rank of current proc in current subcomm

		//std::string fname;              // open inputcprs file
        //int me;     // id of the current processor (rank)
		
        MPI_Status status;
        
	};
	
}

#endif
