#ifndef SET_POINTERS_H
#define SET_POINTERS_H

#include "deto.h"

namespace DETO_NS {

// universal defines inside namespace

#define FLERR __FILE__,__LINE__

#define MIN(A,B) ((A) < (B) ? (A) : (B))
#define MAX(A,B) ((A) > (B) ? (A) : (B))

#define PI 3.1415926535897932384626433832795

class Pointers {

    
public:

    Pointers(DETO *ptr) :
    dto(ptr),
    error(ptr->error),
    inputdeto(ptr->inputdeto),
    universe(ptr->universe),
    output(ptr->output),
    lammpsIO(ptr->lammpsIO),
    optimize(ptr->optimize),
    sims(ptr->sims),
    update(ptr->update),
    plog(ptr->plog),
    screen(ptr->screen){}
    virtual ~Pointers() {}
        
protected:
    DETO *dto;
    Error *&error;
    Inputdeto *&inputdeto;
    Universe *&universe;
    Output *&output;
    LammpsIO *&lammpsIO;
    Optimize *&optimize;
    Simulations *&sims;
    Update *&update;
    FILE *&screen;
    FILE *&plog;
};
     
}

#endif
