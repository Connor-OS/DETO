#ifndef DETO_H
#define DETO_H

#include "stdio.h"
#include "stdlib.h"
#include <string>
#include "mpi.h"
#include <unistd.h>   //just for the sleep() function
#include <vector>
#include <fstream>

#define MASTER 0


namespace DETO_NS {
    
    //using namespace LAMMPS_NS;
    
    
	class DETO {
		
        
        
	public:
        class Error *error;
        class Inputdeto *inputdeto;
        class Universe *universe;
        class Output *output;
        class LammpsIO *lammpsIO;
        class Optimize *optimize;
        class Simulations *sims;
        class Update *update;
        
        FILE *screen;                  // screen output
        
        bool wplog;      // if true, each processor writes a processor specific log for debug
        FILE *plog;                  // processor-specific log file output
        std::string plogfname;

	DETO(int,char**);  // constructor
	~DETO();           // destructor
		
	void printall();
        
    private:
        int me;
	};
}

#endif
