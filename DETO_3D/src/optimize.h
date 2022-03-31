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

        struct Chi_map // structure containing the information mapping chi's to type and other properties
        {
            std::vector<std::string> properties;
            std::vector<std::vector<double>> values;
            // std::vector<int> types; //probably this vector will be included to streamline seting the atom types
        };
        struct Chi_map chi_map; // Instance of chi_map

        
        void read_chimap(std::string);
        void initalize_chi();
        void optrun();
        void printall();
        
	private:
        std::string err_msg, read_string, word;
        double chi_max; // max value of chi specified in chi_map
        int nchi;   //number of chi values in the optimization
        
        int natoms; // number of atoms in LAMMPS
        std::vector<dobule> chi;   // chi values per atom (those with type not in chi_map will be assigne chi > chi_max
        std::vector<int> aID;  // ID of all atoms in LAMMPS
        std::vector<int> atype;  // type of all atoms in LAMMPS
        
        
        //void check_name_unique(std::string);

		//std::string fname;              // open inputcprs file
        //int me;     // id of the current processor (rank)
		
	};
	
}

#endif
