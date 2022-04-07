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
        void optrun();
        void printall();
        
	private:
        struct Chi_map // structure containing the information mapping chi to type and other properties
        {
            std::vector<double> chis;
            std::vector<int> types;
            std::vector<std::string> properties;
            std::vector<std::vector<double>> values;
        };
        struct Chi_map chi_map; // Instance of chi_map

        std::vector<double> vol_constraint; // volume constraint number between 0 and 1 per material
        std::vector<double> local_vol_constraint; // local volume constraint number between 0 and 1 and radius per material
        std::vector<double> local_vol_radius;
        
        std::string err_msg, read_string, word;
        double chi_max; // max value of chi specified in chi_map
        int nchi;   //number of chi values in the optimization
        int nmat;   //number of material values in the optimization

        int natoms; // number of atoms in LAMMPS
        int nlocal; // number of atoms in current proc from LAMMPS
        std::vector<double> chi;   // chi values per atom (those with type not in chi_map will be assigne chi > chi_max

        double* aID;  // ID of all atoms in LAMMPS in current proc
        double* atype;  // type of all atoms in LAMMPS in current proc
        int* IDuns; // unsortd IDs of all atoms in LAMMPS
        int* nID_each;  //array with number of IDs in each processor in current subcomm
        
        int nploc;  // number of processors in local subcomm
        //void check_name_unique(std::string);
        
        int* IDpos;     // position of local tID array in submaster's unsorted list of IDs
        
        int key;   // rank of current proc in current subcomm

		//std::string fname;              // open inputcprs file
        //int me;     // id of the current processor (rank)
		
        MPI_Status status;
        
	};
	
}

#endif
