#ifndef SET_output_H
#define SET_output_H

#include "pointers.h"
#include "stdio.h"
#include "mpi.h"
#include "stdlib.h"
#include <string>
#include <vector>
#include <iostream>
#include <sstream>

#define MASTER 0

namespace DETO_NS {
    class Output : protected Pointers {
        
    public:
        
        Output(class DETO *);
        ~Output();
        
        /*
         int th_every;  // thermo is written every so many steps
        
        std::vector<std::string> dumpID; // IDs of dumps in current subcomm
        std::vector<int> dump_every; //frequency of each dump in current subcomm
        std::vector<std::string> dump_string; // string called by dump  "wrtie_dump group custom filename id type..."
        std::vector<bool> dump_first;  //recognises if dump was called already or this is first call


        
        std::vector<std::string> th_qtts; // quantities to print in thermo file

#ifdef MASKE_WITH_NUFEB
	int nufeb_submaster_rank;
#endif
        */
        void createplog(std::string);
        void toplog(std::string);
        /*
        
        void createthermo(std::string);
        void writethermo(void);
        void add_thqtt(std::string);
        void writedump(int);
        void add_dump(std::string,int,std::string);
         */

        
    private:
        int me;
    };
}

#endif
