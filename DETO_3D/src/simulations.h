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
        
        
        void printall();
        
	private:
        //std::string err_msg, read_string, word;
        
        
        //void check_name_unique(std::string);

		//std::string fname;              // open inputcprs file
        //int me;     // id of the current processor (rank)
		
	};
	
}

#endif
