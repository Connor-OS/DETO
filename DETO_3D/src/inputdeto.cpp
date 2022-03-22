#include "inputdeto.h"
#include "error.h"
#include "universe.h"
#include "lammpsIO.h"
#include "optimize.h"
#include "simulations.h"
/*

#include "chemistry.h"
#include "solution.h"
#include "fix.h"
#include "krun.h"
#include "output.h"
#include "relax.h"
#include "setconc.h"
#include "store.h"
#ifdef MASKE_WITH_SPECIATION
#include "spec.h"
#endif
*/

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
        std::string msg = "ERROR: no inputdeto file specified\n";
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
// Reading inputmsk file (not from restart)
void Inputdeto::file()
{
	
    if (me == MASTER) fprintf(screen, "\nLooking for \"DETOrestart.dat\" file in current folder...\n");
    
    std::string read_string;
    
    std::ifstream inFile("DETOrestart.dat");
    
    
    
    // Is this a restart or not?
    if (!inFile.is_open()) {
        inFile.open (fname.c_str(), std::ifstream::in);
        if (me == MASTER) fprintf(screen,"\"DETOrestart.dat\" not found. Proceeding to input file \"%s\"...\n",fname.c_str());
    }
    else{ /*
        if (me == MASTER) fprintf(screen,"\"MSKrestart.dat\" file found. Checking if it is active...\n");
        //look for the STATUS line, which must be the first non-comment line
        bool isrestartfound = false;
        while (!inFile.eof() && !isrestartfound) {
            inFile >> read_string;
            std::transform(read_string.begin(), read_string.end(), read_string.begin(), ::tolower);  //convert string lower case
            if (strncmp(read_string.c_str(), "#", 1) == 0) {
                std::getline (inFile, totrash); // lines starting with # are comments, which are disregarded here
            }
            else if (strcmp(read_string.c_str(), "status") == 0) {
                inFile >> read_string;
                std::transform(read_string.begin(), read_string.end(), read_string.begin(), ::tolower);  //convert string lower case
                if (strcmp(read_string.c_str(), "active") == 0){
                    isrestart = true;
                    isrestartfound = true;
                    if (me == MASTER) fprintf(screen,"\"MSKrestart.dat\" file is active; let's go with it...\n");
                }
                else{
                    if (me == MASTER) fprintf(screen,"\"MSKrestart.dat\" is not active (to make it active, the first non-comment line in \"MSKrestart.dat\" should be: \n STATUS ACTIVE \nProceeding to inputcprs file  \"%s\"...\n",fname.c_str());
                    inFile.close();
                    inFile.open (fname.c_str(), std::ifstream::in);
                    isrestartfound = true;
                }
            }
            else {
                if (me == MASTER) fprintf(screen,"\"MSKrestart.dat\" is not active (to make it active, the first non-comment line in \"MSKrestart.dat\" should be: \n STATUS ACTIVE \nProceeding to inputcprs file...\n");
                inFile.close();
                inFile.open (fname.c_str(), std::ifstream::in);
                isrestartfound = true;
            }
        }
           */
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
        std::getline (inFile, read_string);
        if (!read_string.empty()){
            //std::istringstream ss(read_string);
            //ss >> firstWord;
            execline(read_string);
        }
    }
    
    inFile.close();

    if (me == MASTER) fprintf(screen,"DONE Reading input file");
    MPI_Barrier(MPI_COMM_WORLD);
}




// ---------------------------------------------------------------
// reading initial configuration from xyz lammps-dump-like file
void Inputdeto::execline(std::string read_string)
{
    
    bool getout = false;
    
    // each line is read word by word. Some words are keywords after which a number of elements is expected. As soon as the symbol # is found the rest of line is trashed being a comment
    std::istringstream lss(read_string);
    while (lss >> word && !getout) {
        if (strncmp(word.c_str(), "#", 1) == 0) break;
        /*else if (strcmp(word.c_str(), "step") == 0 && isrestart) {
            lss >> msk->step;
            msk->doublestep = (double) msk->step;
            if (me == MASTER) fprintf(screen, "\nlast step in restart was %d", msk->step);
            foundstep = true;
        }
        else if (strcmp(word.c_str(), "tempo") == 0 && isrestart) {
            lss >> msk->tempo;
            if (me == MASTER) fprintf(screen, "\ntempo = %f", msk->tempo);
            foundtempo = true;
        }*/
        else if (strcmp(word.c_str(), "subcomm") == 0) {
            int nsubs;
            lss >> nsubs;
            universe->nsc = nsubs;
            //record names of subcommunicators
            for (int i=0; i<universe->nsc; i++) {
                std::stringstream ss;
                ss << i;
                std::string str = ss.str();
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
                    std::string temp_string;
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
                //fprintf(screen, "In %s LAMMMPS doing: %s\n",(universe->SCnames[universe->color]).c_str(),(lss.str()).c_str());
                //std::getline (inFile, read_string);
                //std::istringstream iss(read_string);
                
                std::string subs, newstring;
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
                std::string err_msg;
                err_msg = "ERROR: lammps not active on a processor";
                error->errsimple(err_msg);
            }
        }
        else if (strcmp(word.c_str(), "opt_map_chi") == 0) {
            std::string mapfname;
            lss >> mapfname;
            optimize->read_chimap(mapfname);
        }
        /*else if (strcmp(word.c_str(), "read_potentials") == 0) {
            std::string potfname;
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
                    std::getline (potFile, read_string2);
                    lammpsIO->lammpsdo(read_string2);
                }
            }
        }*/
        else if (strcmp(word.c_str(), "simulation") == 0) {
            std::string read_string2;
            std::getline(lss, read_string2);
            if (me==MASTER) {
                fprintf(screen,"\n inputdeto -- Adding simulation:  %s \n",read_string2.c_str());
            }
            sims->add(read_string2);
        }
        else if (strcmp(word.c_str(), "add_attribute") == 0) {
            std::string read_string2;
            std::getline(lss, read_string2);
            sims->add_attribute(read_string2);
        }
        else if (strcmp(word.c_str(), "add_objective") == 0) {
            std::string read_string2;
            std::getline(lss, read_string2);
            sims->add_objective(read_string2);
        }
        /*
        else if (strcmp(word.c_str(), "real_types") == 0) {
            int tt;
            while (lss >> tt){  msk->Rtypes.push_back(tt); }
            if (me == MASTER) {
                fprintf(screen, "Recorded real particle types: ");
                for (int i=0; i< msk->Rtypes.size(); i++) fprintf(screen, " %d ",msk->Rtypes[i]);
                fprintf(screen, "\n");
            }
        }
        else if (strcmp(word.c_str(), "trial_types") == 0) {
            int tt;
            while (lss >> tt){  msk->Ttypes.push_back(tt); }
            if (me == MASTER) {
                fprintf(screen, "Recorded trial particle types: ");
                for (int i=0; i< msk->Ttypes.size(); i++) fprintf(screen, " %d ",msk->Ttypes[i]);
                fprintf(screen, "\n");
            }
        }
        else if (strcmp(word.c_str(), "chemDB") == 0){
            lss >> read_string;
            chem->readDB(read_string);
        }
        else if (strcmp(word.c_str(), "kB") == 0) {
            lss >> std::scientific >> msk->kB;
            if (me == MASTER) fprintf(screen, "\n kB = %e", msk->kB);
        }
        else if (strcmp(word.c_str(), "hpl") == 0) {
            lss >> msk->hpl;
            if (me == MASTER) fprintf(screen, "\n hpl = %e", msk->hpl);
        }
        else if (strcmp(word.c_str(), "sol_start") == 0){
            if (me == MASTER) fprintf(screen, "\nFound sol_start\n");
            lss>>read_string;
            if (strcmp(read_string.c_str(), "uniform")== 0){
                solution->soltype = "uniform";
                solution->molins.clear();
                solution->conc.clear();
                solution->concdV.clear();
                int nm;
                lss>>nm;
                for (int i=0; i<nm; i++){
                    double conci;
                    lss>>read_string>>conci;
                    solution->addmol(read_string,conci);
                }
                
                bool comment_found = false;
                while (lss && !comment_found){
                    lss >> read_string;
                    if (strcmp(read_string.c_str(), "Temp")== 0) lss >> solution->Temp;
                    else if (strcmp(read_string.c_str(), "DH_A")== 0) lss >> solution->DH_A;
                    else if (strcmp(read_string.c_str(), "DH_B")== 0) lss >> solution->DH_B;
                    else if (strcmp(read_string.c_str(), "voidV")== 0) lss >> solution->voidV;
                    else if (strcmp(read_string.c_str(), "dV")== 0) lss >> solution->dVtype >> solution->dV;
                    else if (strcmp(read_string.c_str(), "dVvoidV")== 0) lss >> solution->dVvoidV;
                    else if (strcmp(read_string.c_str(), "unitC")== 0) lss >> solution->unitC;
                    else if (strncmp(read_string.c_str(), "#", 1) == 0) {
                        comment_found = true;
                        getout=true;
                    }
                    else {
                        std::string msg = "ERROR: keyword in sol_start  not valid: \""+read_string+"\" \n";
                        error->errsimple(msg);
                    }
                }
                solution->computeNmol();
            }
            else {
                std::string msg = "ERROR: start_sol type not valid. \""+read_string+"\" was found, but only uniform is allowed \n";
                error->errsimple(msg);
            }
            if (me == MASTER) {
                fprintf(screen, "\nDone reading the solution. Temperature is %f and concentrations and numbers for each molecule in box and dV are listed below: \n",solution->Temp);
                for (int i=0; i<chem->Nmol; i++) {
                     fprintf(screen, "  %s   %f  %f %f %f \n",chem->molnames[i].c_str(),chem->mol_cins[i],chem->mol_nins[i],chem->mol_cindV[i],chem->mol_nindV[i]);
                }
            }
        }
        else if (strcmp(word.c_str(), "fix") == 0){
            std::string read_string3;
            std::getline(lss, read_string3);
            if (me==MASTER) {
                fprintf(screen,"\n inputmsk -- Adding fix:  %s \n",read_string3.c_str());
            }
            fix->add(read_string3);
            
            fprintf(screen,"\n inputmsk -- I am processor %d , part of subcomm %s , and I have %d KMC-free fixes and %d Cont fixes defined now \n",me,(universe->SCnames[universe->color]).c_str(),(int)(fix->fKMCtype.size()),(int)(fix->Ctype.size()));
        }
        else if (strcmp(word.c_str(), "store") == 0){
            std::string read_string3;
            std::getline(lss, read_string3);
            if (me==MASTER) {
                fprintf(screen,"\n inputmsk -- Storing this quantity:  %s \n",read_string3.c_str());
            }
            store->add(read_string3);
            
            fprintf(screen,"\n inputmsk -- Processor %d , part of subcomm %s , has stored: \n",me,(universe->SCnames[universe->color]).c_str());
            for (int i=0; i<store->RegNames.size(); i++) fprintf(screen,"Proc %d - %s\n",me,store->RegCmd[i].c_str());
            for (int i=0; i<store->LatNames.size(); i++) fprintf(screen,"Proc %d - %s\n",me,store->LatCmd[i].c_str());
            for (int i=0; i<store->MinNames.size(); i++) fprintf(screen,"Proc %d - %s\n",me,store->MinCmd[i].c_str());
            for (int i=0; i<store->MinNames.size(); i++) fprintf(screen,"Proc %d - %s\n",me,store->MinModCmd[i].c_str());
        }
        else if (strcmp(word.c_str(), "thermo") == 0){
            if (me==MASTER) msk->wthermo = true;
            lss>>output->th_every;
            lss>>msk->th_fname;
            bool comment_found = false;
            lss >> read_string;
            while (lss && !comment_found){
                if (strncmp(read_string.c_str(), "#", 1) != 0){
                    output->add_thqtt(read_string);
                }
                else  {
                    comment_found = true;
                    getout=true;
                }
                lss >> read_string;
            }
            if (me==MASTER) output->createthermo(msk->th_fname);
        }
        else if (strcmp(word.c_str(), "dump") == 0){
            if (lammpsIO->lammps_active) {
                bool comment_found = false;
                std::string subcname;
                lss >> subcname;
                if (strcmp(subcname.c_str(),(universe->SCnames[universe->color]).c_str())==0 || strcasecmp(subcname.c_str(),"all")==0 ) {
                    std::string read_string,dumpID,dump_group,dump_style,dump_fname,dump_string,tdumpID;
                    int dump_every;
                    lss >> dumpID >> dump_every >> dump_group >> dump_style >> dump_fname;
                
                    // add subcomm name to dump name ----
                    dump_fname = dump_fname+"."+universe->SCnames[universe->color];
                    
                    dump_string = dump_group+" "+dump_style+" 1 "+dump_fname;
                    
                    lss >> read_string;
                    while (lss && !comment_found){
                        if (strncmp(read_string.c_str(), "#", 1) != 0){
                            dump_string = dump_string+" "+read_string;
                        }
                        else  {
                            comment_found = true;
                            getout=true;
                        }
                        lss >> read_string;
                    }
                    output->add_dump(dumpID,dump_every,dump_string);
                }
                else {
                    break;
                }
            }
            else {
                getout=true;
            }
        }
        else if (strcmp(word.c_str(), "setconc") == 0) {
            std::string sname,molname;
            double molconc;
            int every;
            lss >> sname >> every >> molname >> molconc;
            
            bool flag_ctr = false;
            std::string counter,ctr_mol,boxdV;
            ctr_mol = "none";
            lss >> counter;
            if (strcmp(counter.c_str(),"counter")==0) {
                flag_ctr = true;
                lss >> ctr_mol;
                lss >> boxdV;
            }
            else boxdV = counter;
            
            if (strcmp(boxdV.c_str(),"box")!=0 && strcmp(boxdV.c_str(),"box+dV")!=0) {
                std::string msg = "ERROR: setconc command must end with either \"box\" or \"box+dV\", instead \""+boxdV+"\" was found\n";
                error->errsimple(msg);
            }
            setconc->add_conc(sname,every,molname,molconc,flag_ctr,ctr_mol,boxdV);
        }
        else if (strcmp(word.c_str(), "relax") == 0) {
            if (lammpsIO->lammps_active) {
                bool comment_found = false;
                int every;
                std::string rid,relaxer,rlx_string,rlx_style,rlx_modify,read_string;
                lss >> rid;
                lss >> every;
                lss >> relaxer;
                if (strcmp(relaxer.c_str(), "minimize") == 0){
                    rlx_string = "minimize";
                    for (int i=0;i<4; i++) {
                        lss >> read_string;
                        rlx_string = rlx_string+" "+read_string;
                    }
                    lss >> rlx_style;
                    rlx_modify = "min_modify";
                    lss >> read_string;
                    while (lss && !comment_found){
                        if (strncmp(read_string.c_str(), "#", 1) != 0){
                            rlx_modify = rlx_modify+" "+read_string;
                        }
                        else  {
                            comment_found = true;
                            getout=true;
                        }
                        lss >> read_string;
                    }
                }
                else {
                    std::string msg = "ERROR: only minimize implemented as a relax mode thus far. No nvt, npt, or anything else yet...";
                    error->errsimple(msg);
                }
                
                relax->add_rlx(rid,every,relaxer,rlx_string,rlx_style,rlx_modify);
            }
        }
        else if (strcmp(word.c_str(), "spec") == 0) {
            #ifdef MASKE_WITH_SPECIATION
                std::string msg = "ERROR: MASKE must be installed with the WITH_SPECIATION option to use this command";
                error->errsimple(msg);
            #endif
            std::string id, type, database, solvent;
            int every, nsolvents;
            double molar_mass;
            lss >> id;
            lss >> type;
            lss >> every;
            lss >> database;
            lss >> nsolvents;
            std::vector<std::string> solvent_names;
            for (int i = 0; i < nsolvents; i++) {
                lss >> solvent;
                solvent_names.push_back(solvent);
            }
            std::vector<double> solvent_molar_masses;
            for (int i = 0; i < nsolvents; i++) {
                lss >> molar_mass; // in g/mol
                solvent_molar_masses.push_back(molar_mass);
            }
            #ifdef MASKE_WITH_SPECIATION
                spec->add_spec(id, type, every, database, solvent_names, solvent_molar_masses);
            #endif
	}
        else if (strcmp(word.c_str(), "Krun") == 0) {
            double deltat;
            lss >> deltat;
            krun->proceed(deltat);
        }
         */
        else{
            std::string msg = "ERROR: command unknown in input or restart file: "+word;
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
