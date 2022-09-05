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

/*

// ---------------------------------------------------------------
// create thermo file (master only)
void Output::createthermo(string fname)
{
    thermo = fopen (fname.c_str(),"w");
    fprintf(thermo,"Step\tTime");
    for (int i=0; i<th_qtts.size(); i++) fprintf(thermo,"\t%s",th_qtts[i].c_str());
    fprintf(thermo,"\n");
    fclose(thermo);

#ifdef MASKE_WITH_NUFEB
    // Find out which process is the submaster for nufeb's subcomm
    if (me == MASTER) {
      for (int i = 0; i < fix->aCtype.size(); i++) {
	if (fix->aCtype[i] == "nufeb") {
	  int rank = 0;
	  for (int j = 0; j < universe->SCnames.size(); j++) {
	    if (universe->SCnames[j] == fix->aCscom[i]) break;
	    rank += universe->SCnp[j];
	  }
	  nufeb_submaster_rank = rank;
	}
      }
    }
#endif
}


// ---------------------------------------------------------------
// add quantity to output in thermo file
void Output::add_thqtt(string name)
{
    th_qtts.push_back(name);
}


// ---------------------------------------------------------------
// write new entry in thermo file (master only)
void Output::writethermo(void)
{
    if (me == MASTER){
        thermo = fopen (msk->th_fname.c_str(),"a");
        fprintf(thermo,"%d",msk->step);
        fprintf(thermo,"\t%e",msk->tempo);
    }
    
    for (int i=0; i<th_qtts.size(); i++) {
        std::istringstream iss(th_qtts[i]);
        string token;
        std::getline(iss, token, '_');   // reads whatever before the first underscore
        if (strcmp(token.c_str(),"conc")==0) {
            iss>>token; // now this contains molecule name
            // locate molecule position in chem
            for (int j=0; j<chem->Nmol; j++) {
                if (strcmp((chem->molnames[j]).c_str(),token.c_str())==0) {
                    if (me == MASTER) fprintf(thermo,"\t%e",chem->mol_cins[j]);
                }
            }
        } else if (strcmp(token.c_str(),"lmp")==0) {
            string vc; //either v or c, to identify lammps variable or compute
            std::getline(iss, vc, '_');
            iss>>token; // now this contains a variable or compute name from lammps
            double varc = 0.;
            
            // the first subcomm with lammps active records the lmp variable or compute value
            if (universe->flampSC == universe->color){
                string tolmp = "run 0";
                lammpsIO->lammpsdo(tolmp);
                 if (strcmp(vc.c_str(),"v")==0) {
                     fprintf(screen,"PROC %d varname %s",me,token.c_str());
                     varc = *((double *) lammps_extract_variable(lammpsIO->lmp,(char *)token.c_str(),0));
                 }
                 else if (strcmp(vc.c_str(),"c")==0) {
                     varc = *((double *) lammps_extract_compute(lammpsIO->lmp,(char *)token.c_str(),0,0));
                 }
            }
        
            // If first subcomm with lammps active is not the MASTER's one, communicate value to master
            if (universe->flampSC != 0){
            }
            //      if universe->flampSC == universe->color   && key==0
            //          send variable value out
            //      else
            //          recv from universe->flampID
            
            if (me == MASTER) fprintf(thermo,"\t%e",varc);
        }
#ifdef MASKE_WITH_NUFEB
	else if (strcmp(token.c_str(),"nufeb")==0) {
	  string var;
	  iss >> var;
	  if (var == "pH") {
	    if (fix_nufeb->init_flag && universe->key == 0) {
	      double ph = 0;
	      if (fix_nufeb->setup_flag)
		ph = fix_nufeb->kinetics->sh[0];
	      MPI_Send(&ph,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
	    }
	    if (me == MASTER && nufeb_submaster_rank > 0) {
	      double ph = 0;
	      MPI_Recv(&ph,1,MPI_DOUBLE,nufeb_submaster_rank,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	      fprintf(thermo,"\t%e",-std::log10(ph));
	    }
	  } else if (var.substr(0, 5) == "conc[") {
	    std::size_t closesb = var.find("]");
	    if (closesb != string::npos) {
	      int i = std::stoi(var.substr(5,closesb));
	      if (fix_nufeb->init_flag && universe->key == 0) {
		double conc = 0;
		if (fix_nufeb->setup_flag)
		  conc = fix_nufeb->kinetics->nus[i][0];
		MPI_Send(&conc,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
	      }
	      if (me == MASTER && nufeb_submaster_rank > 0) {
		double conc;
		MPI_Recv(&conc,1,MPI_DOUBLE,nufeb_submaster_rank,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		fprintf(thermo,"\t%e",conc);
	      }
	    }
	  }
	}
#endif
    }
    if (me == MASTER) fprintf(thermo,"\n");
    if (me == MASTER) fclose(thermo);
    
}

*/

// ---------------------------------------------------------------
// record new dump
void Output::add_dump(int devery, string dfile, string dstring, int n_fitest)
{
    dump_every.push_back(devery);
    dstring = "write_dump " + dstring;
    dump_string.push_back(dstring);
    dump_file.push_back(dfile);
    dump_fitest.push_back(n_fitest);
    dump_first.push_back(true);
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
        MPI_Barrier(MPI_COMM_WORLD);
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
// write new entry in thermo file, currently a place holder function that allows to get a basic function
//TODO: add user defined thermo properties
void Output::writethermo(int step, double* objective_eval, int* fitness)
{
    // fprintf(screen,"writing %d %f to thermo file",step,objective_eval);
    if(step == 0) {
        thermo.open("thermo.objective");
        thermo << "step,best,mean,worst" << std::endl;
    }
    else thermo.open("thermo.objective",std::ios_base::app);
    int pop_size = optimize->pop_size;
    double mean = 0;
    for(int i=0; i<pop_size; i++) {
        mean += objective_eval[i];
    }
    mean = mean/pop_size;

    thermo << step << "," << objective_eval[fitness[0]] << "," << mean << "," << objective_eval[fitness[optimize->pop_size-1]] << std::endl;


    fprintf(screen,"best\tmean\tworst\n");
    fprintf(screen,"%.3f\t%.3f\t%.3f\n",objective_eval[fitness[0]],mean,objective_eval[fitness[optimize->pop_size-1]]);

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