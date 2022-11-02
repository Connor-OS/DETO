#include "output.h"
#include "universe.h"
#include "error.h"
#include "lammpsIO.h"
#include "optimize.h"
#include <string.h>
#include <iostream>
#include <fstream>

using namespace DETO_NS;

// ---------------------------------------------------------------
// Initialize class
Output::Output(DETO *deto) : Pointers(deto)
{
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    
    if (me == MASTER) fprintf(screen,"Generating output class\n");
    
    // th_every = 0;
    wrestart = false;
     


}



// ---------------------------------------------------------------
// Class destructor
Output::~Output()
{
    
}


// ---------------------------------------------------------------
// creates processor specific log file
void Output::createplog(string fname)
{
    plog = fopen (fname.c_str(),"w");
    fclose(plog);
    toplog("BEGINNING OF PROCESSOR SPECIFIC LOG");
    
}


// ---------------------------------------------------------------
// Write to processor-specific log file
void Output::toplog(string msg)
{
    plog = fopen (dto->plogfname.c_str(),"a");
    fprintf(plog,"%s\n",msg.c_str());
    fclose(plog);
}


// ---------------------------------------------------------------
// record new dump
void Output::add_dump(int devery, string dfile, string dstring, int n_fitest)
{
    dump_every.push_back(devery);
    dstring = dstring;
    dump_string.push_back(dstring);
    dump_file.push_back(dfile);
    dump_fitest.push_back(n_fitest);
    dump_first.push_back(true);
}

// ---------------------------------------------------------------
// write new entry in dump file
void Output::writedump(int step)
{
    for(int i=0; i<dump_every.size(); i++) {
        if(step == 0 && me==MASTER) {
            std::ofstream dump;
            dump.open(dump_file[i]);
            dump.close();
        }
         // write what to do when no n_fitest provided
        if(step%dump_every[i] == 0) {
            lammpsIO->lammpsdo("dump " + std::to_string(i) + " " +  dump_string[i]);
            std::string tolmp;
            tolmp = "dump_modify " + std::to_string(i) + " every 1 first yes append yes";
            lammpsIO->lammpsdo(tolmp);
            // write entry
            lammpsIO->lammpsdo("run 0");
            // close dump
            tolmp = "undump " + std::to_string(i);
            lammpsIO->lammpsdo(tolmp);
        }       
    }
}

// ---------------------------------------------------------------
// write new entry in dump file
void Output::writedump(int step, int pop_size, int* fitness)
{
    for(int i=0; i<dump_every.size(); i++) {
        if(step == 0 && me==MASTER) {
            std::ofstream dump;
            dump.open(dump_file[i]);
            dump.close();
        }
        MPI_Barrier(MPI_COMM_WORLD); // write what to do when no n_fitest provided
        if(step%dump_every[i] == 0) {
            if(dump_fitest[i] > 0) { //dump a proportion of the fitest solutions
                MPI_Barrier(MPI_COMM_WORLD);
                for(int j=0; j<dump_fitest[i]; j++) {
                    int k=0;
                    while(fitness[j] >= optimize->pop_sizeps_cum[k+1] && k < universe->nsc-1) k++;
                    if(universe->color == k) {
                        optimize->load_chi(fitness[j]-optimize->pop_sizeps_cum[k]);
                        lammpsIO->lammpsdo(dump_string[i] + " modify append yes");
                    }
                    MPI_Barrier(MPI_COMM_WORLD);
                }
            }
        }
    }
}


// ---------------------------------------------------------------
// write new restart file
void Output::writerestart()
{
    if(universe->color == 0) {
        lammpsIO->lammpsdo("write_data " + restart_file + " nocoeff");
    }
}


// ---------------------------------------------------------------
// write new entry in thermo file, currently a place holder function that allows to get a basic function
//TODO: add user defined thermo properties
void Output::writethermo(int step, double* objective_eval, int* fitness)
{
    // fprintf(screen,"writing %d %f to thermo file",step,objective_eval);
    if(step == 0) {
        thermo.open("thermo.objective");
        thermo << "step,best,mean,worst";
        for(int i=0; i<optimize->thermo_string.size(); i++) {
            thermo << "," << optimize->thermo_string[i];
        }
        thermo << std::endl;
    }
    else thermo.open("thermo.objective",std::ios_base::app);
    int pop_size = optimize->pop_size;
    double mean = 0;
    for(int i=0; i<pop_size; i++) {
        mean += objective_eval[i];
    }
    mean = mean/pop_size;

    thermo << step << "," << objective_eval[fitness[0]] << "," << mean << "," << objective_eval[fitness[optimize->pop_size-1]];
    for(int i=0; i<optimize->thermo_string.size(); i++) {
        thermo << "," << optimize->thermo_val[fitness[0]][i];
    }
    thermo << std::endl;


    fprintf(screen,"best\tmean\tworst");
    for(int i=0; i<optimize->thermo_string.size(); i++) {
        fprintf(screen,"\t%s",optimize->thermo_string[i].c_str());
    }
    fprintf(screen,"\n");
    fprintf(screen,"%.3f\t%.3f\t%.3f",objective_eval[fitness[0]],mean,objective_eval[fitness[optimize->pop_size-1]]);
    for(int i=0; i<optimize->thermo_string.size(); i++) {
        fprintf(screen,"\t%.3f",optimize->thermo_val[fitness[0]][i]);
    }
    fprintf(screen,"\n");

    thermo.close();
}


// ---------------------------------------------------------------
// Printing info about the inputcprs class (possibly useful for debugging)
void Output::printall()
{
    fprintf(screen, "\n---------ALL ABOUT OUTPUT----------\n");
    fprintf(screen,"Dumps %li\n",dump_every.size());
    for(int i=0; i<dump_every.size(); i++) {
        fprintf(screen,"Every %i dump: %s for fitest: %i\n",dump_every[i],dump_string[i].c_str(),dump_fitest[i]);
    }
    fprintf(screen, "\n---------------------------------------\n\n");
}