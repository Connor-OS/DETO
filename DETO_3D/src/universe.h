#ifndef UNIVERSE_H
#define UNIVERSE_H

#include "pointers.h"
#include "mpi.h"
#include <string>
#include <stdlib.h>
#include <vector>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <unistd.h>   //just for the sleep() function


#define MASTER 0

namespace DETO_NS {
    using std::vector;
    using std::string;
	
	class Universe : protected Pointers {
	public:
		//int narg;                    // # of command args
		//char **arg;                  // parsed args for command
		
		 Universe(class DETO *);
		~Universe();
		
		void create();       // creates the universe of sub-communicators
		void printall();
        
        int nsc;            // number of subcommunicators defined
        int color;          // processor-specific color, which identifies what subcomm the processor pertains to
        int key;            // processor-specific key, which indicates its rank within the subcomm
       
        MPI_Status status;
        
        vector<string> SCnames;    //vector containing the names of the subcommunicators
        vector<int> SCnp;    //vector containing the number of processors in each subcommunicator

        vector<int> subMS;    //vector containing the IDs of the submasters of each subcommunicator
        
        int *color_each; //vector containing the color, viz the subcomm id, of each processor
        //int *subcID;    // array containing the proc id of the masters of each subcomm
        
        MPI_Comm subcomm;
        int me;     // id of the current processor (rank)
        int nprocs;     // number of processor in the universe (Comm_size)

		
	private:
		//string fname;              // open inputcprs file

	};
	
}

#endif
