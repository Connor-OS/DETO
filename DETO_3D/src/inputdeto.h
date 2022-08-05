#ifndef SET_inputdeto_H
#define SET_inputdeto_H

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

namespace DETO_NS {
	using std::string;
	using std::getline;
	using std::to_string;
	
	class Inputdeto : protected Pointers {
	public:
		int narg;                    // # of command args
		char **arg;                  // parsed args for command
		
		 Inputdeto(class DETO *, int, char **);
		~Inputdeto();
		
		void file();       // process the input script
		void printall();
        void execline(string);

	
        bool isrestart; // a flag to let the code know whether it is running a restart or not
        
	private:
		string fname;              // open input file
        int me;     // id of the current processor (rank)
		string inconfig,totrash,err_msg,read_in,firstWord,word;              // input config file (data or xyz type)
        int tstep;  // lammps timestep read from xyz file
        bool time_strict, time_loose, time_first, time_last;    //flags to know how to read time from XYZ file
        bool foundSubcomm, foundstep,foundtempo;
        double(ru);
        
	};
	
}

#endif
