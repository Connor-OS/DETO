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
		
	/*class Memory *memory;          // memory allocation functions
     */
        class Error *error;
        class Inputdeto *inputdeto;
        class Universe *universe;
        class Output *output;
        class LammpsIO *lammpsIO;
        class Optimize *optimize;
        class Simulations *sims;
        /*
        class Chemistry *chem;
        class Solution *solution;
        class Fix *fix;
        class Krun *krun;
        class Randm *randm;
        class Fix_delete *fix_del;
        
        class Fix_Cfoo *fix_cfoo;
        class Relax *relax;
        class Setconc *setconc;
        class Fix_nucleate *fix_nucl;
        class Store *store;
#ifdef MASKE_WITH_NUFEB
	class Fix_nufeb *fix_nufeb;
#endif
#ifdef MASKE_WITH_SPECIATION
	class Spec *spec;
#endif
     */

        /*
        int step;       // the current simulation step
        double doublestep;  // same as step, but in double version, as needed to calculate the frequency of invoking transition (see isasctive functions in transition files)
        double tempo;   // current simulation time
        double kB;      // Boltzmann constant (user-defined in input, otherwise initialised to 1)
        double hpl;      // Planck constant (user-defined in input, otherwise initialised to 1)
        */
        
        FILE *screen;                  // screen output
        
        
        bool wplog;      // if true, each processor writes a processor specific log for debug
        FILE *plog;                  // processor-specific log file output
        std::string plogfname;
        /*
        bool wthermo;      // if true, the MASTER writes the thermo file
        FILE *thermo;                  // thermo outup file: written by master only
        std::string th_fname;

        std::vector<int> Rtypes;    // vector of types associated to real particles
        std::vector<int> Ttypes;    // vector of types associated to trial particles

	bool nulog_flag;     // true if nufeb log is to be generated
	std::ofstream nulog; // nufeb log file

	
         */
        
	DETO(int,char**);  // constructor
	~DETO();           // destructor
		
	void printall();
        
    private:
        int me;
	};
}

#endif
