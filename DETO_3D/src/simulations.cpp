#include "simulations.h"
#include <sstream>
#include "error.h"
#include "lammpsIO.h"
#include "universe.h"
#include "optimize.h"
#include "output.h"

#include <string.h>

using namespace DETO_NS;

Simulations::Simulations(DETO *deto) : Pointers(deto)
{

    // The inputcprs class is run by all processors in COMM_WORLD. This reads the id of the processor
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    n_sims = 0;

}



// ---------------------------------------------------------------
// Class destructor
Simulations::~Simulations()
{
    
}



// ---------------------------------------------------------------
void Simulations::add(string read_string)
{
// record name in a vector
    string sim_name, sim_type, repstr, repYN;
    n_sims++;
    
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
    cstgs_crit.push_back("NULL");
    
    if (strcmp(sim_type.c_str(), "cstgs") == 0) {
        int pos = cstgs_varname.size()-1;
        string varname, ttype, incty;
        double par1, par2, par3, par4;
        lss >> varname >> ttype;
        cstgs_varname[pos] = varname;
        cstgs_type[pos] = ttype;
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
    
   

    // default values for repeat inputs - will be overwritten if sim is repeat
    sim_rep_vars.push_back(vector<string>());
    sim_rep_val.push_back(vector<vector<double>>());
    n_repeats.push_back(1);

    lss >> repstr >> repYN;
    
    if (strcmp(repstr.c_str(), "repeat") != 0) 
    {
        err_msg = "Error: keyword \"repeat\" must follow simulation type in input script\n";
        error->errsimple(err_msg);
    }
    if (strcmp(repYN.c_str(), "yes") == 0)
    {
        string repfname;
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
    
    
    cstgs_crit_vnms.push_back(vector<string>());
    if (strcmp(sim_type.c_str(), "cstgs") == 0)  {
        int pos = cstgs_varname.size()-1;
        string kcrit;
        lss >> kcrit;
        if (strcmp(kcrit.c_str(), "crit") == 0){
            std::getline(lss, read_string2);
            //fprintf(screen,"DEBUG: getline records this \"%s\"\n",read_string2.c_str());
            string delimiter = "'";
            size_t spos = 0;
            string critStr;
            spos = read_string2.find(delimiter);
            read_string2.erase(0, spos + delimiter.length());
            spos = read_string2.find(delimiter);
            critStr = read_string2.substr(0, spos);
            cstgs_crit[pos] = critStr;
            read_string2.erase(0, spos + delimiter.length());
        
            std::istringstream lss2(read_string2);
            string vnam;
            while (lss2 >> vnam){
                if  (strncmp(vnam.c_str(), "#", 1) == 0) break;
                else cstgs_crit_vnms[pos].push_back(vnam);
            }
        }
        else {
            // TODO: implement error that crit keyword has not been found for cstgs simulation type (crit must be specified for cstgs simulations)
        }
    }

    // default containers for objectives - will be filled out by seperate add_objective command
    sim_obj_names.push_back(vector<vector<string>>());
    sim_obj_LMPnames.push_back(vector<vector<string>>()); // Should these be initialised here or right at the begining of this function?
    sim_obj_val.push_back(vector<vector<double>>());
    for (int i=0;i<n_repeats[n_repeats.size()-1];i++) 
    {  
        sim_obj_names[sim_obj_names.size()-1].push_back(vector<string>());
        sim_obj_LMPnames[sim_obj_LMPnames.size()-1].push_back(vector<string>());
        sim_obj_val[sim_obj_val.size()-1].push_back(vector<double>());
    }

        // default containers for objectives - will be filled out by seperate add_objective command
    sim_sens_names.push_back(vector<vector<string>>());
    sim_sens_LMPnames.push_back(vector<vector<string>>()); // Should these be initialised here or right at the begining of this function?
    sim_sens_val.push_back(vector<vector<double*>>());
    for (int i=0;i<n_repeats[n_repeats.size()-1];i++) 
    {  
        sim_sens_names[sim_obj_names.size()-1].push_back(vector<string>());
        sim_sens_LMPnames[sim_obj_LMPnames.size()-1].push_back(vector<string>());
        sim_sens_val[sim_obj_val.size()-1].push_back(vector<double*>());
    }


    //default container for attributes - will be filled out by seperate add_attribute command
    sim_attributes.push_back(vector<string>());

    ndim = lammpsIO->extract_setting("dimension");
}


// ---------------------------------------------------------------
void Simulations::add_attribute(string read_string)
{
// find simulation id matching specified name
    // add string to a vector of vectors of string (one list of attributes per simulation) at the mathcing simulation ID
    string sim_ID;
    std::istringstream lss(read_string);
    lss >> sim_ID;
    for(int i = 0; i < sim_names.size(); i++){
        if(strcmp(sim_ID.c_str(),sim_names[i].c_str()) == 0){
            string read_string2;
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
void Simulations::add_objective(string read_string)
{
    string sim_ID, obj_name, obj_LMPname;
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

void Simulations::add_sensitivity(string read_string)
{
    string sim_ID, sens_name, sens_LMPname;
    std::istringstream lss(read_string);
    lss >> sim_ID >> sens_name >> sens_LMPname;
    for(int i = 0; i < sim_names.size(); i++){
        if(strcmp(sim_ID.c_str(),sim_names[i].c_str()) == 0){
            for (int j=0; j<n_repeats[i]; j++){
                sim_sens_names[i][j].push_back(sens_name);
                sim_sens_LMPnames[i][j].push_back(sens_LMPname);
                sim_sens_val[i][j].push_back(NULL); //this needs to be a array of size natoms but for now NULL
            }
            return;
        }
    }
    err_msg = "ERROR: objective assigned to unspecified simulation";
    error->errsimple(err_msg);
}


// ---------------------------------------------------------------
void Simulations::read_repeat(string repeatfname)
{
    std::ifstream repeatFile(repeatfname.c_str());
    // vector<double>
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
                                sim_rep_val[simpos].push_back(vector<double>());
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
// Run simulation of given sim_id
void Simulations::run()
{
    for(int i=0; i<sim_attributes.size(); i++) {
        for(int j=0; j<n_repeats[i]; j++) {
            for(int k=0; k<sim_attributes[i].size(); k++) { //TODO: add while loop for cstgs and for loop for repeat
                lammpsIO->lammpsdo(sim_attributes[i][k]);
            }
            for(int k=0; k<sim_obj_names[i][j].size(); k++) {
                sim_obj_val[i][j][k] = *(double *)lammpsIO->extract_varaiable(sim_obj_LMPnames[i][j][k]);            
            }
            for(int k=0; k<sim_sens_names[i][j].size(); k++) {
                sim_sens_val[i][j][k] = (double *)lammpsIO->gather_atom_varaiable(sim_sens_LMPnames[i][j][k]);
            }
        }
    }
}


// ---------------------------------------------------------------
// Run simulation of given sim_id
void Simulations::run_one(int sim)
{
    for(int j=0; j<n_repeats[sim]; j++) {
        for(int k=0; k<sim_attributes[sim].size(); k++) { //TODO: add while loop for cstgs and for loop for repeat
            lammpsIO->lammpsdo(sim_attributes[sim][k]);
        }
        for(int k=0; k<sim_obj_names[sim][j].size(); k++) {
            sim_obj_val[sim][j][k] = *(double *)lammpsIO->extract_varaiable(sim_obj_LMPnames[sim][j][k]);            
        }
    }
}


// ---------------------------------------------------------------
// Printing info about the inputcprs class (possibly useful for debugging)
void Simulations::printall()
{
	fprintf(screen,"\n---------ALL ABOUT Simulations----------\n");
	//fprintf(screen,"inputcprs filename =  %s\n",fname.c_str());
    fprintf(screen,"Number of simulations %d\n",n_sims);
    for(int i=0; i<n_sims; i++) 
    {
        fprintf(screen,"ID: %s\n  Type: %s\n",sim_names[i].c_str(),sim_types[i].c_str());
        // print cstgs parameters
        fprintf(screen,"  cstg type: %s\n  cstg variable: %s\n  Parameters: %.2f %.2f\n  Criterion: %s\n  Lammps vars: ",cstgs_type[i].c_str(),cstgs_varname[i].c_str(),cstgs_par1[i],cstgs_par2[i],cstgs_crit[i].c_str());
        for (int j=0; j<cstgs_crit_vnms[i].size(); j++){
            fprintf(screen,"%s ",cstgs_crit_vnms[i][j].c_str());
        }
        fprintf(screen,"\n");

        // print repeat parameters
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
                fprintf(screen,"  %.2f ",sim_rep_val[i][k][j]);
            }
            fprintf(screen, "\n");
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
