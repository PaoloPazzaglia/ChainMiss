#ifndef MILP_LET__
#define MILP_LET__

//#include <ilcplex/ilocplex.h>

#include <algorithm>   
#include <vector>  
#include <fstream>   
#include <string>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

#include "milp_data.h"

double MILP_WH_K(std::vector<Task> &taskchain, std::vector<WHconstr> &setofmk, OptTarget mytarget);
//double MILP_WH_K_PATHS(std::vector<Task> &taskchain, std::vector<WHconstr> &setofmk, int num_paths, int chain_D);
//int MILP_WH_SN(std::vector<Task> &taskchain, std::vector<WHconstr> &setofmk);

std::string  convert_to_string(const int Number);

#endif
