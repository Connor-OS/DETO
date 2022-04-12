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
    /*
     memory(ptr->memory),
     */
    inputdeto(ptr->inputdeto),
    universe(ptr->universe),
    output(ptr->output),
    lammpsIO(ptr->lammpsIO),
    optimize(ptr->optimize),
    sims(ptr->sims),
    /*
     
    chem(ptr->chem),
    solution(ptr->solution),
    fix(ptr->fix),
    fix_del(ptr->fix_del),
#ifdef MASKE_WITH_NUFEB
    fix_nufeb(ptr->fix_nufeb),
#endif
    krun(ptr->krun),
    randm(ptr->randm),
     */
    plog(ptr->plog),
    screen(ptr->screen){}
    /*
    thermo(ptr->thermo),
    fix_cfoo(ptr->fix_cfoo),
    relax(ptr->relax),
    setconc(ptr->setconc),
#ifdef MASKE_WITH_SPECIATION
    spec(ptr->spec),
#endif
    store(ptr->store),
    fix_nucl(ptr->fix_nucl) {}
        
     */
    virtual ~Pointers() {}
        
protected:
    DETO *dto;
    /*
     Memory *&memory;
     */
    Error *&error;
    Inputdeto *&inputdeto;
    Universe *&universe;
    Output *&output;
    LammpsIO *&lammpsIO;
    Optimize *&optimize;
    Simulations *&sims;
    /*
     
    Chemistry *&chem;
    Solution *&solution;
    Fix *&fix;
    Fix_delete *&fix_del;
    Fix_Cfoo *&fix_cfoo;
#ifdef MASKE_WITH_NUFEB
    Fix_nufeb *&fix_nufeb;
#endif
    Relax *&relax;
#ifdef MASKE_WITH_SPECIATION
    Spec *&spec;
#endif
    Krun *&krun;
    Randm *&randm;
    Setconc *&setconc;
    
    Fix_nucleate *&fix_nucl;
    Store *&store;
     */
    FILE *&screen;
    FILE *&plog;
    /*
    FILE *&thermo;
     */
};
     
}

#endif