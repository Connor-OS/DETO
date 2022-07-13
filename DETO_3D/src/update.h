#ifndef UPDATE_H
#define UPDATE_H

#include "pointers.h"
#include "universe.h"
// #include "simulations.h"
#include <string>
#include <vector>

#define MASTER 0

namespace DETO_NS {
	
	class Update : protected Pointers {
	public:
        Update(class DETO *);
		~Update();

        int me;     // id of the current processor (rank)
		void evaluate_objective(int id);
		void set_opt_type(std::string read_string);

		void genetic(const std::vector<std::vector<double>>& chi_pop,const std::vector<std::vector<int>>& mat_pop,const double* opt_obj_eval);
		void monte_carlo(const std::vector<std::vector<double>>& chi_pop,const std::vector<std::vector<int>>& mat_pop,const double* opt_obj_eval);
		void sensitivity();
		std::vector<std::vector<double>> update_chipop(const std::vector<std::vector<double>>& chi_pop,const std::vector<std::vector<int>>& mat_pop,const double* opt_obj_eval);

        void printall();


	private:

		std::string err_msg, read_string, word;
		int key;   // rank of current proc in current subcomm
		int natoms;


		//Optimisation variables
        std::string opt_type; //the type of optimization to be run. i.e genetic, sensitivity etc..
        int pop_size; // size of population of solutions
        double opt_par1, opt_par2, opt_par3; //paramaters specifict to optimisation types
        //if "genetic" then par1 = crossover rate, part2 = mutation rate
        //if "sensitivity" then par1 = move limit
		const std::vector<std::vector<double>> chi_pop;
        // const double* opt_obj_eval;
		std::vector<std::vector<double>> chi_next; // vector containing the next chi_pop
		std::vector<std::vector<int>> mat_next; // vector containing the next mat_pop

	};
	
}

#endif
