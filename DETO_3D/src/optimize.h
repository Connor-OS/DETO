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
        
        std::map<std::string, std::vector<double> > chi_map; // a map recording the content of the chi_map file provided by the user
        std::map<std::string, std::vector<double> >::iterator it; // map iterator
        
        void read_chimap(std::string);
        void printall();
        
	private:
        std::string err_msg, read_string, word;
        int nchi;   //number of chi values in the optimization
        
        //void check_name_unique(std::string);

		//std::string fname;              // open inputcprs file
        //int me;     // id of the current processor (rank)
		
	};
	
}

#endif
