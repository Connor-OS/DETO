#include "inputdeto.h"
#include "error.h"
#include "deto.h"
#include "universe.h"
#include "lammpsIO.h"
#include "optimize.h"
#include "update.h"
#include "simulations.h"
#include "output.h"
#include <unistd.h>   //just for the sleep() function
//needed to convert strings to lower case
#include <algorithm>
//for cluster Rocket...
#include <string.h>

using namespace DETO_NS;

Inputdeto::Inputdeto(DETO *deto, int argc, char **argv) : Pointers(deto)
{
    // The inputcprs class is run by all processors in COMM_WORLD. This reads the id of the processor
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    
    if (me == MASTER) fprintf(screen,"Generating inputdeto object\n");
    
    if (argc <2) {
        string msg = "ERROR: no inputdeto file specified\n";
        error->errsimple(msg);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    isrestart = false;
    fname = argv[1];
    time_strict=false;
    time_loose=false;
    time_first=false;
    time_last=false;
    ru = 1.;
    foundSubcomm=false;
    foundstep=false;
    foundtempo=false;

}



// ---------------------------------------------------------------
// Class destructor
Inputdeto::~Inputdeto()
{
    
}







// ---------------------------------------------------------------
// Reading inputdeto file (not from restart)
void Inputdeto::file()
{
	
    if (me == MASTER) fprintf(screen, "\nLooking for \"DETOrestart.dat\" file in current folder...\n");
    
    string read_string;
    
    std::ifstream inFile("DETOrestart.dat");
    
    
    
    // Is this a restart or not?
    if (!inFile.is_open()) {
        inFile.open (fname.c_str(), std::ifstream::in);
        if (me == MASTER) fprintf(screen,"\"DETOrestart.dat\" not found. Proceeding to input file \"%s\"...\n",fname.c_str());
    }
    
	
	//crash if the input file cannot be read
	if (!inFile.is_open()) {
        err_msg = "ERROR: cannot read file \""+fname+"\"";
        error->errsimple(err_msg);
	}
    MPI_Barrier(MPI_COMM_WORLD);
	
    
    
   // READ FILE (all processors need to know this)
    while (!inFile.eof()) {
        MPI_Barrier(MPI_COMM_WORLD);
        getline (inFile, read_string);
        if (!read_string.empty()){
            //std::istringstream ss(read_string);
            //ss >> firstWord;
            execline(read_string);
        }
    }
    
    inFile.close();

    if (me == MASTER) fprintf(screen,"DONE Reading input file\n\n");
    MPI_Barrier(MPI_COMM_WORLD);
}




// ---------------------------------------------------------------
// reading initial configuration from xyz lammps-dump-like file
void Inputdeto::execline(string read_string)
{
    
    bool getout = false;
    
    // each line is read word by word. Some words are keywords after which a number of elements is expected. As soon as the symbol # is found the rest of line is trashed being a comment
    std::istringstream lss(read_string);
    while (lss >> word && !getout) {
        if (strncmp(word.c_str(), "#", 1) == 0) break;
        else if (strcmp(word.c_str(), "subcomm") == 0) {
            int nsubs;
            lss >> nsubs;
            universe->nsc = nsubs;
            //record names of subcommunicators
            for (int i=0; i<universe->nsc; i++) {
                std::stringstream ss;
                ss << i;
                string str = ss.str();
                universe->SCnames.push_back(str);
                universe->SCnp.push_back(0);
            }
            
            int nprocs;
            MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
            for (int i=0; i<nprocs; i++){
                int j;    // subcomm id
                j = i % nsubs;
                universe->SCnp[j] = universe->SCnp[j] + 1;
            }
            
            if (me == MASTER){
                for (int i=0; i<nsubs; i++){
                    string temp_string;
                    temp_string = universe->SCnames[i];
                    fprintf(screen, "\nSubcommunicator \"%s\" added (with %d processors)",temp_string.c_str(),universe->SCnp[i]);
                }
                fprintf(screen, "\n");
            }
        
            universe->create();
            
            fprintf(screen,"\nProcessor \"%d\"assigned with colour\"%d\"and key \"%d\"\n",me,universe->color,universe->key);
        
            lammpsIO->create();
            /*foundSubcomm = true;
             */
            
            MPI_Barrier(MPI_COMM_WORLD);
        }
        else if (strcmp(word.c_str(), "lammps") ==0)  {
            if (lammpsIO->lammps_active) {
                string subs, newstring;
                lss >> firstWord;
                // add subcomm name to dump name ----
                if (strcmp(firstWord.c_str(),"dump")==0) {
                    newstring = "dump ";
                    for (int i=0; i<4; i++) {
                        lss >> subs;
                        newstring = newstring+subs+" ";
                    }
                    lss >> subs;
                    subs = subs+"."+universe->SCnames[universe->color];
                    newstring = newstring+subs+" ";
                    while (lss >> subs) {
                        newstring = newstring+subs+" ";
                    }
                    //fprintf(screen, "\nString is \"%s\"\n",newstring.c_str());
                }
                else{
                    newstring = firstWord+" ";
                    while (lss >> subs) {
                        newstring = newstring+subs+" ";
                    }
                    if (strcmp(firstWord.c_str(),"thermo_style")==0) {
                        lammpsIO->lmpThSt = newstring;
                    }
                }
                //---
                //fprintf(screen, "In %s LAMMMPS doing: %s\n",(universe->SCnames[universe->color]).c_str(),newstring.c_str());
                lammpsIO->lammpsdo(newstring);
            }
            else {
                string err_msg;
                err_msg = "ERROR: lammps not active on a processor";
                error->errsimple(err_msg);
            }
        }
        else if (strcmp(word.c_str(), "opt_map_chi") == 0) {
            string mapfname;
            lss >> mapfname;
            optimize->read_chimap(mapfname);
        }
        else if (strcmp(word.c_str(), "read_potentials") == 0) {
            string potfname;
            lss >> potfname;
            std::ifstream potFile(potfname.c_str());
            if (!potFile.is_open()) {
                err_msg = "ERROR: cannot read file \""+potfname+"\"";
                error->errsimple(err_msg);
            }
            else{
                MPI_Barrier(MPI_COMM_WORLD);
               // READ FILE (all processors need to know this)
                while (!potFile.eof()) {
                    MPI_Barrier(MPI_COMM_WORLD);
                    getline (potFile, read_in);
                    // lammpsIO->lammpsdo(read_in);
                    optimize->potentials.push_back(read_in);
                }
            }
        }
        else if (strcmp(word.c_str(), "dump") == 0){
            if (lammpsIO->lammps_active) {
                bool comment_found = false;
                string dump_string,dump_file;
                int dump_every,n_fitest;
                n_fitest = 0;
                lss >> read_in;
                for(int i=0; i<2; i++) {
                    lss >> read_in;
                    dump_string = dump_string+" "+read_in;
                }
                lss >> dump_every;
                lss >> dump_file;
                dump_string = dump_string+" 1 "+dump_file;
                while (lss >> read_in && !comment_found){
                    if(strcmp(read_in.c_str(), "n_fitest") == 0){
                        lss >> n_fitest;
                    }
                    else if (strncmp(read_in.c_str(), "#", 1) != 0){
                        dump_string = dump_string+" "+read_in;
                    }
                    else  {
                        comment_found = true;
                        getout=true;
                    }
                }
                output->add_dump(dump_every,dump_file,dump_string,n_fitest);
            }
        }
        else if (strcmp(word.c_str(), "objective_function") == 0) {
            getline(lss, read_in);
            optimize->obj_function = read_in;
        }
        else if (strcmp(word.c_str(), "simulation") == 0) {
            getline(lss, read_in);
            sims->add(read_in);
        }
        else if (strcmp(word.c_str(), "add_attribute") == 0) {
            getline(lss, read_in);
            sims->add_attribute(read_in);
        }
        else if (strcmp(word.c_str(), "add_objective") == 0) {
            getline(lss, read_in);
            sims->add_objective(read_in);
        }
        else if (strcmp(word.c_str(), "add_sensitivity") == 0) {
            getline(lss, read_in);
            sims->add_sensitivity(read_in);
        }
        else if (strcmp(word.c_str(), "add_constraint") == 0) {
            getline(lss, read_in);
            optimize->add_constraint(read_in);
        }
        else if (strcmp(word.c_str(), "opt_type") == 0) {
            getline(lss, read_in);
            update->set_opt_type(read_in);
        }
        else if (strcmp(word.c_str(), "write_plog") == 0) {
            lss >> read_in;
            if (strcmp(read_in.c_str(),"yes") == 0) {
                dto->wplog = true;
                string fname;
                fname = "p" + to_string(me) + "_S";
                fname = fname + to_string(universe->color) + "_k";
                fname = fname + to_string(universe->key) + ".plog";
                dto->plogfname=fname;
                output->createplog(fname);
            }
            else if (strcmp(read_in.c_str(),"no") == 0) {
                dto->wplog = false;
            }
            else {
                err_msg = "ERROR: Illegal write_plog command (yes or no not \""+read_in+"\")";
                error->errsimple(err_msg);
            }
        }
        else if (strcmp(word.c_str(), "write_lmp_log") == 0) {
            lss >> read_in;
            if (strcmp(read_in.c_str(),"yes") == 0) {
                lammpsIO->wllog = true;
            }
            else if (strcmp(read_in.c_str(),"no") == 0) {
                lammpsIO->wllog = false;
                lammpsIO->lammpsdo("log none");
            }
            else {
                err_msg = "ERROR: Illegal write_lmp_log command (yes or no not \""+read_in+"\")";
                error->errsimple(err_msg);
            }
        }
        else if (strcmp(word.c_str(), "write_restart") == 0) {
            output->wrestart = true;
            lss >> read_in;
            output->restart_file = read_in;
            int r_every;
            while(lss >> r_every) {
                output->restart_every = r_every;
            }
        }
        else if (strcmp(word.c_str(), "thermo_style") == 0) {
            if(!optimize->pop_size) {
                err_msg = "ERROR: Must set opt_type before using thermo_style command";
                error->errsimple(err_msg);
            }
            while(lss >> read_in) {
                optimize->thermo_string.push_back(read_in);
            }
        }
        else{
            string msg = "ERROR: command unknown in input or restart file: "+word;
            error->errsimple(msg);

        }
         
        
    }
        
}




// ---------------------------------------------------------------
// Printing info about the inputcprs class (possibly useful for debugging)
void Inputdeto::printall()
{
	fprintf(screen,"\n---------ALL ABOUT inputdeto----------\n");
	//fprintf(screen,"inputcprs filename =  %s\n",fname.c_str());
	fprintf(screen,"---------------------------------------\n\n");
	
}
