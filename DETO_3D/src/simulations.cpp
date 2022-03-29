#include "simulations.h"
#include <sstream>
#include "error.h"
//#include "universe.h"

//#include "chemistry.h"
//#include "store.h"
//#include "lammpsIO.h"
//#include "error.h"
/*#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
*/

#include <string.h>

using namespace DETO_NS;

Simulations::Simulations(DETO *deto) : Pointers(deto)
{

    // The inputcprs class is run by all processors in COMM_WORLD. This reads the id of the processor
    MPI_Comm_rank(MPI_COMM_WORLD, &me);


}



// ---------------------------------------------------------------
// Class destructor
Simulations::~Simulations()
{
    
}



// ---------------------------------------------------------------
void Simulations::add(std::string read_string)
{
// record name in a vector
    std::string sim_name, sim_type, repstr, repYN;
    
    std::istringstream lss(read_string);
    lss >> sim_name >> sim_type;
    sim_names.push_back(sim_name);
    sim_types.push_back(sim_type);

    // default values for cstgs inputs - will be overwritten if sim is "cstgs"
    cstgs_varname.push_back("NULL");
    cstgs_type.push_back("NULL");
    cstgs_par1.push_back(-1.);
    cstgs_par2.push_back(-1.);
    cstgs_par3.push_back(-1.);
    cstgs_par4.push_back(-1.);
    
    if (strcmp(sim_type.c_str(), "cstgs") == 0) {
        int pos = cstgs_varname.size()-1;
        std::string varname, ttype, incty;
        double par1, par2, par3, par4;
        lss >> varname >> ttype;
        if (strcmp(ttype.c_str(), "binary") == 0) {
            lss >> par1 >> par2 >> par3 >> par4;
            cstgs_par1[pos] = par1;
            cstgs_par2[pos] = par2;
            cstgs_par3[pos] = par3;
            cstgs_par4[pos] = par4;
        }
        if (strcmp(ttype.c_str(), "increment") == 0) {
            lss >> incty;
            if (strcmp(incty.c_str(), "linear") == 0) {
                lss >> par1 >> par2 >> par3;
                par4 = -1.;
            }
            else if (strcmp(incty.c_str(), "power") == 0) {
                lss >> par1 >> par2 >> par3 >> par4;
            }
            else if (strcmp(incty.c_str(), "log") == 0) {
                lss >> par1 >> par2 >> par3;
                par4 = -1.;
            }
            cstgs_par1[pos] = par1;
            cstgs_par2[pos] = par2;
            cstgs_par3[pos] = par3;
            cstgs_par4[pos] = par4;
        }
        
        
    }
    cstgs_crit.push_back("NULL");   // this will be provided separaetly by the user through dedicated "Criterion sim_ID string" command in the input scrupt

    // default values for repeat inputs - will be overwritten if sim is repeat
    sim_rep_vars.push_back(std::vector<std::string>());
    sim_rep_val.push_back(std::vector<std::vector<double>>());
    n_repeats.push_back(1);

    lss >> repstr >> repYN;
    
    if (strcmp(repstr.c_str(), "repeat") != 0) 
    {
        err_msg = "Error: keyword \"repeat\" must follow simulation type in input script\n";
        error->errsimple(err_msg);
    }
    if (strcmp(repYN.c_str(), "yes") == 0)
    {
        std::string repfname;
        lss >> repfname;
        sim_repeat_file.push_back(repfname);
        read_repeat(repfname);
    }
    else if (strcmp(repYN.c_str(), "no") == 0)
    {
        sim_repeat_file.push_back("NULL");
    }
    else 
    {
        err_msg = "Error: repeat must be \"yes\" or \"no\", case sensitive. Instead I found "+repYN+"\n";
        error->errsimple(err_msg);
    }

    // default containers for objectives - will be filled out by seperate add_objective command
    sim_obj_names.push_back(std::vector<std::vector<std::string>>());
    sim_obj_LMPnames.push_back(std::vector<std::vector<std::string>>()); // Should these be initialised here or right at the begining of this function?
    sim_obj_val.push_back(std::vector<std::vector<double>>());
    for (int i=0;i<n_repeats[n_repeats.size()-1];i++) 
    {  
        sim_obj_names[sim_obj_names.size()-1].push_back(std::vector<std::string>());
        sim_obj_LMPnames[sim_obj_LMPnames.size()-1].push_back(std::vector<std::string>());
        sim_obj_val[sim_obj_val.size()-1].push_back(std::vector<double>());
    }

    //default container for attributes - will be filled out by seperate add_attribute command
    sim_attributes.push_back(std::vector<std::string>());
}


// ---------------------------------------------------------------
void Simulations::add_attribute(std::string read_string)
{
// find simulation id matching specified name
    // add string to a vector of vectors of string (one list of attributes per simulation) at the mathcing simulation ID
    std::string sim_ID;
    std::istringstream lss(read_string);
    lss >> sim_ID;
    for(int i = 0; i < sim_names.size(); i++){
        if(strcmp(sim_ID.c_str(),sim_names[i].c_str()) == 0){
            std::string read_string2;
            std::getline(lss, read_string2);
            sim_attributes[i].push_back(read_string2);
            if (me==MASTER){
            // fprintf(screen,"Adding attribute to sim: %s attribute: %s \n",sim_names[i].c_str(),sim_attributes[i][sim_attributes[i].size()-1].c_str());
            }
            return;
        }
    }
    err_msg = "ERROR: attribute assigned to unspecified simulation";
    error->errsimple(err_msg);
    // TODO: Implement error if sim name not found
}


// ---------------------------------------------------------------
void Simulations::add_objective(std::string read_string)
{
    std::string sim_ID, obj_name, obj_LMPname;
    std::istringstream lss(read_string);
    lss >> sim_ID >> obj_name >> obj_LMPname;
    for(int i = 0; i < sim_names.size(); i++){
        if(strcmp(sim_ID.c_str(),sim_names[i].c_str()) == 0){
            for (int j=0; j<n_repeats[i]; j++){
                sim_obj_names[i][j].push_back(obj_name);
                sim_obj_LMPnames[i][j].push_back(obj_LMPname);
                sim_obj_val[i][j].push_back(0.);
            }
            return;
        }
    }
    err_msg = "ERROR: objective assigned to unspecified simulation";
    error->errsimple(err_msg);
}


// ---------------------------------------------------------------
void Simulations::read_repeat(std::string repeatfname)
{
    std::ifstream repeatFile(repeatfname.c_str());
    // std::vector<double>
    bool found_nrep = false;
    int simpos;
    simpos = sim_names.size()-1;
    
    if (!repeatFile.is_open())
    {
        err_msg = "ERROR: cannot read file \"" + repeatfname + "\"";
        error->errsimple(err_msg);
    }
    else
    {
        MPI_Barrier(MPI_COMM_WORLD);
        // READ FILE (all processors need to know this)
        while (!repeatFile.eof())
        {
            MPI_Barrier(MPI_COMM_WORLD);
            std::getline(repeatFile, read_string);
            if (!read_string.empty())
            {
                std::istringstream lss(read_string);
                while (lss >> word)
                {
                    if (strncmp(word.c_str(), "#", 1) == 0) break;
                    else if (strcmp(word.c_str(), "num_rep") == 0)
                    {
                        lss >> n_repeats[simpos];
                        found_nrep = true;
                        break;
                    }
                    else if(strcmp(word.c_str(), "VARS:") == 0)
                    {
                        if (found_nrep == true)
                        {
                            while (lss >> word)
                            {
                                sim_rep_vars[simpos].push_back(word);
                                sim_rep_val[simpos].push_back(std::vector<double>());
                            }
                            for(int i=0; i < n_repeats[simpos]; i++)
                            {
                                std::getline(repeatFile, read_string);
                                if (!read_string.empty())
                                {
                                    std::istringstream lss(read_string);
                                    for (int j = 0; j < sim_rep_vars[simpos].size(); j++)
                                    {
                                        double tmp;
                                        lss >> tmp;
                                        sim_rep_val[simpos][j].push_back(tmp);
                                    }
                                }
                                else
                                {
                                    err_msg = "ERROR: Fewer than specified repeat parameters given";
                                    error->errsimple(err_msg);
                                }
                            }
                        }
                        else
                        {
                            err_msg = "ERROR: num_repeats not found in repeat file before specifying repeat values";
                            error->errsimple(err_msg);
                        }
                    }
                }
            }
        }
    }
}

// ---------------------------------------------------------------
// Printing info about the inputcprs class (possibly useful for debugging)
void Simulations::printall()
{
	fprintf(screen,"\n---------ALL ABOUT Simulations----------\n");
	//fprintf(screen,"inputcprs filename =  %s\n",fname.c_str());
    for(int i=0; i<sim_names.size(); i++) 
    {
        fprintf(screen,"ID: %s\n  Type: %s\n",sim_names[i].c_str(),sim_types[i].c_str());
        // print cstgs parameters
        if(true)    //sim_types[i].compare("cstgs") == 0
        {
            fprintf(screen,"  cstg type: %s\n  cstg variable: %s\n  Parameters: %.2f %.2f\n  Tolerance: %.2f\n  Criterion: %s\n",cstgs_type[i].c_str(),cstgs_varname[i].c_str(),cstgs_par1[i],cstgs_par2[i],cstgs_tol[i],cstgs_crit[i].c_str());
        }
        // print repeat parameters
        if(true)   //n_repeats[i] != 1
        {
            fprintf(screen,"  Repeat: %d\n  Repeat file: %s\n",n_repeats[i],sim_repeat_file[i].c_str()); 
            for (int j=0; j<sim_rep_vars[i].size(); j++)
            {
                fprintf(screen,"  %s  ",sim_rep_vars[i][j].c_str());
            }
            fprintf(screen,"\n");
            for (int j = 0; j < n_repeats[i]; j++)
            {
                for (int k = 0; k < sim_rep_vars[i].size(); k++)
                {
                    fprintf(screen,"  %f ",sim_rep_val[i][k][j]);
                }
                fprintf(screen, "\n");
            }
        }
        else
        {
            fprintf(screen,"Repeat: no\n");
        }
        // print objectives
        fprintf(screen,"\nObjectives:\n");
        for(int j=0; j<sim_obj_names[i][0].size(); j++)
        {
            fprintf(screen,"  Objective: %s\n  Lammps name: %s\n",sim_obj_names[i][0][j].c_str(),sim_obj_LMPnames[i][0][j].c_str());
        }
        // print attributes
        fprintf(screen,"\nAttributes:\n");
        for(int j=0; j<sim_attributes[i].size(); j++)
        {
            fprintf(screen,"  %s\n",sim_attributes[i][j].c_str());
        }
        fprintf(screen,"\n---------------------------------------\n");
    }
	fprintf(screen,"---------------------------------------\n\n");
	
}
