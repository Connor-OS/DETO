#ifndef UPDATE_H
#define UPDATE_H

#include "pointers.h"
#include "universe.h"
#include <string>
#include <vector>
#include <memory>

#define MASTER 0

namespace DETO_NS {
	using std::vector;
	using std::string;
	using std::make_unique;
	using std::make_shared;
	using std::unique_ptr;
	using std::shared_ptr;
	
	class Update : protected Pointers {
	public:
        Update(class DETO *);
		~Update();

        int me;     // id of the current processor (rank)
		int key;   // rank of current proc in current subcomm
		void evaluate_objective(int id);
		void set_opt_type(string read_string);

		void genetic(const vector<vector<double>>& chi_pop,const vector<vector<int>>& mat_pop,const double* opt_obj_eval, const int* fitness);
		void monte_carlo(const vector<vector<double>>& chi_pop,const vector<vector<int>>& mat_pop,const double* opt_obj_eval);
		void perturbation(const double* chi,const int* mat,const double* opt_obj_eval,vector<double>& dchi,int sim);
		void gradient_descent(const double* chi,const int* mat,vector<double>& dchi);
		void sensitivity_update(const vector<double> dchi, const double* chi, const int* mat);
		void update_chipop(vector<vector<double>>& chi_pop,vector<vector<int>>& mat_pop,const double* opt_obj_eval, const int* fitness);
		void filter_sensitivities(const double* chi, vector<double>& dchi, const int* bonds, const int nbonds);

        void printall();

		std::ofstream sens_file;

	private:

		string err_msg, read_string, word;
		int natoms;

		//Optimisation variables
        vector<string> opt_type; //the type of optimization to be run. i.e genetic, sensitivity etc..
        string opt_style; //the style for the optimisation type to be run. i.e genetic style = tornement or roulette 
		int pop_size; // size of population of solutions
        double opt_par1, opt_par2, opt_par3; //paramaters specifict to optimisation types
		int gen_elitism; // number of best solutions retained during each step of a genetic optimization
        //if "genetic" then par1 = crossover rate, part2 = mutation rate
        //if "sensitivity" then par1 = move limit
		const vector<vector<double>> chi_pop;
        // const double* opt_obj_eval;
		double* chi_next; // vector containing the next chi_pop
		vector<vector<int>> mat_next; // vector containing the next mat_pop

		vector<double> chi;
		vector<int> mat;

		MPI_Status status;

		double* dchips;
		double* dchi_;

		double* update_obj_evalps;
		double* update_obj_eval;
	};
	
}

#endif
