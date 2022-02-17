#ifndef ERROR_H
#define ERROR_H

#include "pointers.h"
#include "stdio.h"
#include "mpi.h"
#include "stdlib.h"
#include <string>

#define MASTER 0

namespace DETO_NS {
    class Error : protected Pointers {
		
	public:
        
		Error(class DETO *);
		~Error();
		
        void errsimple(std::string);
        
    private:
        int me;
	};
}

#endif
