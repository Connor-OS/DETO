#ifndef SET_inputmsk_H
#define SET_inputmsk_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "pointers.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include "mpi.h"
#include <vector>

#define MASTER 0

namespace MASKE_NS {
	
	class Inputmsk : protected Pointers {
	public:
		int narg;                    // # of command args
		char **arg;                  // parsed args for command
		
		 Inputmsk(class MASKE *, int, char **);
		~Inputmsk();
		
		void file();       // process the input script
		void printall();
        void execline(std::string);

	
        bool isrestart; // a flag to let the code know whether it is running a restart or not
        
	private:
		std::string fname;              // open inputcprs file
        int me;     // id of the current processor (rank)
		std::string inconfig,totrash,err_msg,read_string2,firstWord,word;              // input config file (data or xyz type)
        int tstep;  // lammps timestep read from xyz file
        bool time_strict, time_loose, time_first, time_last;    //flags to know how to read time from XYZ file
        bool foundSubcomm, foundstep,foundtempo;
        double(ru);
        
	};
	
}

#endif
