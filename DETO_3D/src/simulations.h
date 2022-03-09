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
        std::vector<std::string> sim_types; //vector holding sim types
        std::vector<bool> sim_is_repeat; //vector of bool specifying if repeat
        std::vector<std::string> sim_repeat_file; //vector holding path to repeat file data

        std::vector<std::vector<std::string>> sim_attributes; //vector sim attributes

        void printall();
        void add(std::string);
        
        void add_attribute(std::string);
        
	private:
        //std::string err_msg, read_string, word;
        
        
        //void check_name_unique(std::string);

		//std::string fname;              // open inputcprs file
        //int me;     // id of the current processor (rank)
		
	};
	
}

#endif
