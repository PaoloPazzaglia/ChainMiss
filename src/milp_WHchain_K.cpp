#include <ilcplex/ilocplex.h>
#include <string>
#include <algorithm>   
#include <vector>  
#include <fstream>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

#include "milp_data.h"
#include "milp_WHchain.h"

#define __DEBUG_MILP__ 1
#define TOL 0.001
#define TOL_OFFS 0.001

using namespace std;

ILOSTLBEGIN

typedef IloArray<IloNumVarArray>   NumVarMatrix;
typedef IloArray<IloIntVarArray>   IntVarMatrix;
typedef IloArray<IntVarMatrix>     IntVar3Matrix;


double MILP_WH_K(vector<Task> &taskchain, vector<WHconstr> &setofmk, OptTarget mytarget)
{
	//-----------------------------------------------------------------------------
	// PROBLEM PARAMETERS
	//-----------------------------------------------------------------------------

	// Inputs from user
	const int NUMBER_OF_PATHS = 2; 

	// Number of tasks in a chain
	const int NUMBER_OF_TASKS_IN_CHAIN = taskchain.size();

	// Big-M (to represent infinity)
	const double BIGM = INT_MAX;

	//-----------------------------------------------------------------------------
	// START MILP DESIGN
	//-----------------------------------------------------------------------------

	IloEnv env;

	// MILP output
	double MILP_out;


#ifdef __DEBUG_MILP__
	cout << "[MILP] Setting up variables...";
#endif
	try
	{
		IloModel model(env);

		//----------------------------------------------------------------------------
		// VARIABLES DEFINITION
		//----------------------------------------------------------------------------

		// Release offset of a task 
		IloNumVarArray OFFS(env, NUMBER_OF_TASKS_IN_CHAIN);
		for (int t = 0; t < NUMBER_OF_TASKS_IN_CHAIN; t++) {
			string name = "OFFS" + convert_to_string(t);
			int T = taskchain.at(t).period;
			OFFS[t] = IloNumVar(env, 0.0, T - TOL_OFFS, name.c_str());
		}

		// Index of effective job
		IntVarMatrix EFFECTIVEJOB(env, NUMBER_OF_TASKS_IN_CHAIN);
		for (int t = 0; t < NUMBER_OF_TASKS_IN_CHAIN; t++) {
			EFFECTIVEJOB[t] = IloIntVarArray(env, NUMBER_OF_PATHS);
			for (unsigned int l = 0; l < NUMBER_OF_PATHS; l++) {
				string name = "EID" + convert_to_string(t) + convert_to_string(l);
				EFFECTIVEJOB[t][l] = IloIntVar(env, 0.0, UINT16_MAX, name.c_str());
			}
		}

		// Number of redundant jobs until a new input is available (R jobs)
		IntVarMatrix REDUNDHITS(env, NUMBER_OF_TASKS_IN_CHAIN);
		for (int t = 0; t < NUMBER_OF_TASKS_IN_CHAIN; t++) {
			REDUNDHITS[t] = IloIntVarArray(env, NUMBER_OF_PATHS - 1);
			for (unsigned int l = 0; l < NUMBER_OF_PATHS - 1; l++) {
				string name = "nRED" + convert_to_string(t) + convert_to_string(l);
				REDUNDHITS[t][l] = IloIntVar(env, 0.0, UINT16_MAX, name.c_str());
			}
		}

		// Number of missed jobs after the first input overwrite (M1 jobs)
		IntVarMatrix MISSWNEWINPUT(env, NUMBER_OF_TASKS_IN_CHAIN);
		for (int t = 0; t < NUMBER_OF_TASKS_IN_CHAIN; t++) {
			MISSWNEWINPUT[t] = IloIntVarArray(env, NUMBER_OF_PATHS);
			for (unsigned int l = 0; l < NUMBER_OF_PATHS; l++) {
				string name = "nMISSni" + convert_to_string(t) + convert_to_string(l);
				MISSWNEWINPUT[t][l] = IloIntVar(env, 0.0, UINT16_MAX, name.c_str());
			}
		}

		// Number of V jobs
		IntVarMatrix VOIDHITS(env, NUMBER_OF_TASKS_IN_CHAIN);
		for (int t = 0; t < NUMBER_OF_TASKS_IN_CHAIN; t++) {
			VOIDHITS[t] = IloIntVarArray(env, NUMBER_OF_PATHS - 1);
			for (unsigned int l = 0; l < NUMBER_OF_PATHS - 1; l++) {
				string name = "nINC" + convert_to_string(t) + convert_to_string(l);
				VOIDHITS[t][l] = IloIntVar(env, 0.0, UINT16_MAX, name.c_str());
			}
		}

		// Number of missed jobs after completion of next effective job of producer task (M2 jobs)
		IntVarMatrix MISSAFTEREFFECTIVE(env, NUMBER_OF_TASKS_IN_CHAIN);
		for (int t = 0; t < NUMBER_OF_TASKS_IN_CHAIN; t++) {
			MISSAFTEREFFECTIVE[t] = IloIntVarArray(env, NUMBER_OF_PATHS);
			for (unsigned int l = 0; l < NUMBER_OF_PATHS; l++) {
				string name = "nMISSV" + convert_to_string(t) + convert_to_string(l);
				MISSAFTEREFFECTIVE[t][l] = IloIntVar(env, 0.0, UINT16_MAX, name.c_str());
			}
		}

		// Auxiliary variable: there is at least one V job of task t between the jobs of paths p and p+1
		IntVarMatrix boolVOIDJOBS(env, NUMBER_OF_TASKS_IN_CHAIN);
		for (int t = 0; t < NUMBER_OF_TASKS_IN_CHAIN; t++) {
			boolVOIDJOBS[t] = IloIntVarArray(env, NUMBER_OF_PATHS - 1);
			for (unsigned int l = 0; l < NUMBER_OF_PATHS - 1; l++) {
				string name = "boolHV" + convert_to_string(t) + convert_to_string(l);
				boolVOIDJOBS[t][l] = IloIntVar(env, 0.0, 1.0, name.c_str());
			}
		}

		// Auxiliary variable: length of subsequence between two blocks of misses is <= k
		IntVar3Matrix boolLENGTHK(env, NUMBER_OF_TASKS_IN_CHAIN);
		for (int t = 0; t < NUMBER_OF_TASKS_IN_CHAIN; t++) {
			boolLENGTHK[t] = IntVarMatrix(env, 2 * NUMBER_OF_PATHS);
			for (unsigned int l = 0; l < 2 * NUMBER_OF_PATHS; l++) {
				boolLENGTHK[t][l] = IloIntVarArray(env, 2 * NUMBER_OF_PATHS);
				for (unsigned int p = 0; p < 2 * NUMBER_OF_PATHS; p++) {
					string name = "boolLEN" + convert_to_string(t) + convert_to_string(l) + convert_to_string(p);
					boolLENGTHK[t][l][p] = IloIntVar(env, 0.0, 1.0, name.c_str());
				}
			}
		}


#ifdef __DEBUG_MILP__
		cout << "DONE!" << endl;
		//-----------------------------------------------------------------------------------------------------------------------------------------  
		cout << "[MILP] Starting constraints" << endl;
#endif


		//----------------------------------------------------------------------------
		// CONSTRAINT 1
		// Constraining variables for head task of the chain		
		// Offset of head task is 0
		model.add(OFFS[0] == 0);
		// First effective job of head task has index 0
		model.add(EFFECTIVEJOB[0][0] == 0);


		//----------------------------------------------------------------------------
		// CONSTRAINT 2
		// Encoded in definition of OFFS


		//----------------------------------------------------------------------------
		// CONSTRAINT 3
		// Head task cannot have redundant jobs
		for (int p = 0; p < NUMBER_OF_PATHS - 1; p++) {
			model.add(REDUNDHITS[0][p] == 0);
		}
		// Head task cannot have "misses after effective job of producer task" (it has no producer!)
		for (int p = 0; p < NUMBER_OF_PATHS; p++) {
			model.add(MISSAFTEREFFECTIVE[0][p] == 0);
		}
		// Tail task cannot have void jobs
		for (int p = 0; p < NUMBER_OF_PATHS - 1; p++) {
			model.add(VOIDHITS[NUMBER_OF_TASKS_IN_CHAIN - 1][p] == 0);
			model.add(boolVOIDJOBS[NUMBER_OF_TASKS_IN_CHAIN - 1][p] == 0);
		}


		//----------------------------------------------------------------------------
		// CONSTRAINT 4
		// Checking if there exist void hits at level of task t
		// Note that task tail cannot have void hits
		for (int t = 0; t < NUMBER_OF_TASKS_IN_CHAIN - 1; t++) {
			for (int p = 0; p < NUMBER_OF_PATHS - 1; p++) {

				int Tt = taskchain.at(t).period;
				int Dt = taskchain.at(t).deadline;

				model.add(VOIDHITS[t][p] <= boolVOIDJOBS[t][p] * BIGM);
				model.add(VOIDHITS[t][p] >= boolVOIDJOBS[t][p]);
			}
		}


		//----------------------------------------------------------------------------
		// CONSTRAINT 5
		// Effective job of i-th task of the chain must start before the
		// next hit job of (i-1)th task of the chain completes with a different output 
		for (int t = 1; t < NUMBER_OF_TASKS_IN_CHAIN; t++) {
			for (int p = 0; p < NUMBER_OF_PATHS - 1; p++) {

				int Tt = taskchain.at(t).period;
				int Dt = taskchain.at(t).deadline;

				int Tt1 = taskchain.at(t - 1).period;
				int Dt1 = taskchain.at(t - 1).deadline;

				// The activation of EFFECTIVEJOB_tp occurs after (or at) the completion of EFFECTIVEJOB_(t-1)p
				model.add(OFFS[t] + Tt * EFFECTIVEJOB[t][p] >= 
					OFFS[t - 1] + Tt1 * EFFECTIVEJOB[t - 1][p] + Dt1);

				// The activation of EFFECTIVEJOB_tp occurs before the completion of EFFECTIVEJOB_(t-1)(p+1)
				model.add(OFFS[t] + Tt * EFFECTIVEJOB[t][p] <=
					OFFS[t - 1] + Tt1 * EFFECTIVEJOB[t - 1][p + 1] + Dt1 - TOL);
			}
		}


		//----------------------------------------------------------------------------
		// CONSTRAINT 6
		// Between the end of the effective job of task t-1 of path p, and the beginning of the 
		// effective job of task t of the same path p, task t-1 cannot have INCOMPLHITS
		for (int t = 1; t < NUMBER_OF_TASKS_IN_CHAIN; t++) {
			for (int p = 0; p < NUMBER_OF_PATHS - 1; p++) {

				int Tt = taskchain.at(t).period;
				int Dt = taskchain.at(t).deadline;

				int Tt1 = taskchain.at(t - 1).period;
				int Dt1 = taskchain.at(t - 1).deadline;

				model.add(OFFS[t - 1] + (EFFECTIVEJOB[t - 1][p] + REDUNDHITS[t - 1][p]
					+ MISSWNEWINPUT[t - 1][p] + 1) * Tt1 + Dt1 - TOL >=
					OFFS[t] + EFFECTIVEJOB[t][p] * Tt - (1 - boolVOIDJOBS[t - 1][p]) * BIGM);
			}
		}


		//----------------------------------------------------------------------------
		// CONSTRAINT 7
		// A task t may have redundant hits after a effective job, until task t-1 has produced a new output
		// Head task has no redundant hits (already defined above)
		for (int t = 1; t < NUMBER_OF_TASKS_IN_CHAIN; t++) {
			for (int p = 0; p < NUMBER_OF_PATHS - 1; p++) {

				int Tt = taskchain.at(t).period;
				int Dt = taskchain.at(t).deadline;

				int Tt1 = taskchain.at(t - 1).period;
				int Dt1 = taskchain.at(t - 1).deadline;

				// The activation of the last redundant hit must occur before the end of the first job of
				// producer task that works on new data and that successfully completes
				model.add(OFFS[t] + Tt * (EFFECTIVEJOB[t][p] + REDUNDHITS[t][p]) <=
					OFFS[t - 1] + Tt1 * (EFFECTIVEJOB[t - 1][p] + REDUNDHITS[t - 1][p]
						+ MISSWNEWINPUT[t - 1][p] + 1)
					+ Dt1 - TOL + (1 - boolVOIDJOBS[t - 1][p]) * BIGM);

				model.add(OFFS[t] + Tt * (EFFECTIVEJOB[t][p] + REDUNDHITS[t][p]) <=
					OFFS[t - 1] + Tt1 * (EFFECTIVEJOB[t - 1][p] + REDUNDHITS[t - 1][p]
						+ MISSWNEWINPUT[t - 1][p] + MISSAFTEREFFECTIVE[t - 1][p + 1] + 1)
					+ Dt1 - TOL + boolVOIDJOBS[t - 1][p] * BIGM);
			}
		}


		//----------------------------------------------------------------------------
		// CONSTRAINT 8
		// Between the end of the effective job of task t-1 of path p, and the beginning of the 
		// effective job of task t of the same path p, task t may miss MISSAFTEREFFECTIVE jobs
		for (int t = 1; t < NUMBER_OF_TASKS_IN_CHAIN; t++) {
			for (int p = 0; p < NUMBER_OF_PATHS; p++) {

				int Tt = taskchain.at(t).period;
				int Dt = taskchain.at(t).deadline;

				int Tt1 = taskchain.at(t - 1).period;
				int Dt1 = taskchain.at(t - 1).deadline;

				// floor function of the number of instances of task t between the end of 
				// EFFECTIVEJOB_(t-1)p and the activation of EFFECTIVEJOB_tp
				model.add(MISSAFTEREFFECTIVE[t][p] >=
					(OFFS[t] + EFFECTIVEJOB[t][p] * Tt -
					(OFFS[t - 1] + EFFECTIVEJOB[t - 1][p] * Tt1 + Dt1)) / Tt - 1 + TOL);

				model.add(MISSAFTEREFFECTIVE[t][p] <=
					(OFFS[t] + EFFECTIVEJOB[t][p] * Tt -
					(OFFS[t - 1] + EFFECTIVEJOB[t - 1][p] * Tt1 + Dt1)) / Tt);
			}
		}


		//----------------------------------------------------------------------------
		// CONSTRAINT 9
		// The index of the effective job of path p+1 is equal to the index of the effective job of path p
		// plus REDUNDHITS + MISSWNEWINPUT + VOIDHITS + MISSAFTEREFFECTIVE + 1
		for (int t = 0; t < NUMBER_OF_TASKS_IN_CHAIN; t++) {
			for (int p = 0; p < NUMBER_OF_PATHS - 1; p++) {

				model.add(EFFECTIVEJOB[t][p + 1] ==
					EFFECTIVEJOB[t][p] + REDUNDHITS[t][p] + MISSWNEWINPUT[t][p]
					+ VOIDHITS[t][p] + MISSAFTEREFFECTIVE[t][p + 1] + 1);
			}
		}

		// Index of effective job of path p+1 is greater than index of effective job of path p
		for (int t = 0; t < NUMBER_OF_TASKS_IN_CHAIN; t++) {
			for (int p = 0; p < NUMBER_OF_PATHS - 1; p++) {
				model.add(EFFECTIVEJOB[t][p + 1] >= 1 + EFFECTIVEJOB[t][p]);
			}
		}


		//----------------------------------------------------------------------------
		// CONSTRAINT 10
		// A task cannot miss more than m_consec misses
		for (int t = 0; t < NUMBER_OF_TASKS_IN_CHAIN; t++) {
			for (int p = 0; p < NUMBER_OF_PATHS - 1; p++) {
				model.add(MISSWNEWINPUT[t][p] <= setofmk.at(t).mconsec);
				model.add(MISSAFTEREFFECTIVE[t][p] <= setofmk.at(t).mconsec);
			}
			model.add(MISSAFTEREFFECTIVE[t][NUMBER_OF_PATHS - 1] <= setofmk.at(t).mconsec);
		}


		//----------------------------------------------------------------------------
		// CONSTRAINT 11
		// If there are no void jobs, MISSWNEWINPUT and MISSAFTEREFFECTIVE occur side by side
		// thus mconsec must be enforced for the whole sequence
		for (int t = 0; t < NUMBER_OF_TASKS_IN_CHAIN; t++) {
			for (int p = 0; p < NUMBER_OF_PATHS - 1; p++) {

				int Tt = taskchain.at(t).period;
				int Dt = taskchain.at(t).deadline;

				model.add(MISSWNEWINPUT[t][p] + MISSAFTEREFFECTIVE[t][p + 1] <=
					setofmk.at(t).mconsec + boolVOIDJOBS[t][p] * BIGM);
			}
		}
		

		//----------------------------------------------------------------------------
		// CONSTRAINT 12 & CONSTRAINT 13
		// Miss and hit patterns must satisfy also the (m,k) constraints
		for (int t = 0; t < NUMBER_OF_TASKS_IN_CHAIN; t++) {

			// Checking sequences starting from MISSAFTEREFFECTIVE
			for (int s = 0; s < NUMBER_OF_PATHS - 1; s++) {

				IloExpr LENGTHSEQ(env);
				IloExpr NUMMISSES(env);

				LENGTHSEQ += MISSAFTEREFFECTIVE[t][s];
				NUMMISSES += MISSAFTEREFFECTIVE[t][s];

				for (int p = s; p < NUMBER_OF_PATHS - 1; p++) {

					// Check until MISSNEWINPUT
					LENGTHSEQ += 1 + REDUNDHITS[t][p] + MISSWNEWINPUT[t][p];
					NUMMISSES += MISSWNEWINPUT[t][p];

					for (int i = 0; i < setofmk.at(t).mk.size(); i++) {

						int m = setofmk.at(t).mk.at(i).m;
						int k = setofmk.at(t).mk.at(i).k;

						model.add(LENGTHSEQ <= k + (1 - boolLENGTHK[t][2 * s][2 * p]) * BIGM);
						model.add(LENGTHSEQ >= k + 1 - boolLENGTHK[t][2 * s][2 * p] * BIGM);

						// Adding necessary and sufficient constraints to check (m,k)
						model.add(NUMMISSES <= m + (1 - boolLENGTHK[t][2 * s][2 * p]) * BIGM);
						model.add(LENGTHSEQ - NUMMISSES >= k - m - boolLENGTHK[t][2 * s][2 * p] * BIGM);
					}

					// Continue checking until MISSAFTEREFFECTIVE
					LENGTHSEQ += VOIDHITS[t][p] + MISSAFTEREFFECTIVE[t][p + 1];
					NUMMISSES += MISSAFTEREFFECTIVE[t][p + 1];

					for (int i = 0; i < setofmk.at(t).mk.size(); i++) {

						int m = setofmk.at(t).mk.at(i).m;
						int k = setofmk.at(t).mk.at(i).k;

						model.add(LENGTHSEQ <= k + (1 - boolLENGTHK[t][2 * s][2 * p + 1]) * BIGM);
						model.add(LENGTHSEQ >= k + 1 - boolLENGTHK[t][2 * s][2 * p + 1] * BIGM);

						// Adding necessary and sufficient constraints to check (m,k)
						model.add(NUMMISSES <= m + (1 - boolLENGTHK[t][2 * s][2 * p + 1]) * BIGM);
						model.add(LENGTHSEQ - NUMMISSES >= k - m - boolLENGTHK[t][2 * s][2 * p + 1] * BIGM);
					}
				}
			}

			// Checking sequences starting from MISSNEWINPUT
			for (int s = 0; s < NUMBER_OF_PATHS - 1; s++) {

				IloExpr LENGTHSEQ(env);
				IloExpr NUMMISSES(env);

				LENGTHSEQ += MISSWNEWINPUT[t][s];
				NUMMISSES += MISSWNEWINPUT[t][s];

				for (int p = s; p < NUMBER_OF_PATHS - 1; p++) {

					// Check until MISSAFTEREFFECTIVE
					LENGTHSEQ += VOIDHITS[t][p] + MISSAFTEREFFECTIVE[t][p + 1];
					NUMMISSES += MISSAFTEREFFECTIVE[t][p + 1];

					for (int i = 0; i < setofmk.at(t).mk.size(); i++) {

						int m = setofmk.at(t).mk.at(i).m;
						int k = setofmk.at(t).mk.at(i).k;

						model.add(LENGTHSEQ <= k + (1 - boolLENGTHK[t][2 * s + 1][2 * p]) * BIGM);
						model.add(LENGTHSEQ >= k + 1 - boolLENGTHK[t][2 * s + 1][2 * p] * BIGM);

						// Adding necessary and sufficient constraints to check (m,k)
						model.add(NUMMISSES <= m + (1 - boolLENGTHK[t][2 * s + 1][2 * p]) * BIGM);
						model.add(LENGTHSEQ - NUMMISSES >= k - m - boolLENGTHK[t][2 * s + 1][2 * p] * BIGM);
					}
					
					if (p < NUMBER_OF_PATHS - 2) {
						// Check until MISSNEWINPUT
						LENGTHSEQ += 1 + REDUNDHITS[t][p+1] + MISSWNEWINPUT[t][p+1];
						NUMMISSES += MISSWNEWINPUT[t][p+1];

						for (int i = 0; i < setofmk.at(t).mk.size(); i++) {

							int m = setofmk.at(t).mk.at(i).m;
							int k = setofmk.at(t).mk.at(i).k;

							model.add(LENGTHSEQ <= k + (1 - boolLENGTHK[t][2 * s + 1][2 * p + 1]) * BIGM);
							model.add(LENGTHSEQ >= k + 1 - boolLENGTHK[t][2 * s + 1][2 * p + 1] * BIGM);

							// Adding necessary and sufficient constraints to check (m,k)
							model.add(NUMMISSES <= m + (1 - boolLENGTHK[t][2 * s + 1][2 * p + 1]) * BIGM);
							model.add(LENGTHSEQ - NUMMISSES >= k - m - boolLENGTHK[t][2 * s + 1][2 * p + 1] * BIGM);
						}
					}					
				}
			}
		}


#ifdef __DEBUG_MILP__
		cout << "Constraints DONE." << endl;
		//-----------------------------------------------------------------------------------------------------------------------------------------  
		cout << "[MILP] Objective Funtion:" << endl;
#endif

		//-----------------------------------------------------------------------------
		// OBJECTIVE FUNCTION
		//-----------------------------------------------------------------------------

		int tailtask_id = NUMBER_OF_TASKS_IN_CHAIN - 1;
		int Tt = taskchain.at(tailtask_id).period;
		int Dt = taskchain.at(tailtask_id).deadline;

		IloNumVar OBJ(env, -INT_MAX, INT_MAX);

		switch (mytarget) {

		case MAXIMIZE_LATENCY: // Maximize end-to-end latency of effective path
			model.add(OBJ <= OFFS[tailtask_id] + EFFECTIVEJOB[tailtask_id][0] * Tt + Dt);
			break;

		case MAXIMIZE_DATAAGE: // Maximize data age
			model.add(OBJ <= OFFS[tailtask_id] + EFFECTIVEJOB[tailtask_id][1] * Tt);
			break;

		case MAXIMIZE_UPDATE_INT: // Maximize update interval
			model.add(OBJ <= (EFFECTIVEJOB[tailtask_id][1] - EFFECTIVEJOB[tailtask_id][0]) * Tt);
			break;

		case MINIMIZE_UPDATE_INT: // Minimize update interval
			model.add(OBJ <= -(EFFECTIVEJOB[tailtask_id][1] - EFFECTIVEJOB[tailtask_id][0]) * Tt);
			break;

		default:
			env.error() << "Unknown optimization target" << endl;
			throw(-1);
		}

		model.add(IloMaximize(env, OBJ));


#ifdef __DEBUG_MILP__
		cout << "DONE." << endl;
		//-----------------------------------------------------------------------------------------------------------------------------------------  
#endif

		env.out() << "Rows populated" << endl;

		//-----------------------------------------------------------------------------
		// SOLVER PARAMETERS
		//-----------------------------------------------------------------------------

		IloCplex cplex(model);
		cplex.exportModel("qcpex1.lp");
		// Optimize the problem and obtain solution.

		// Set minimum GAP to 1%
		cplex.setParam(IloCplex::EpGap, 1e-2);

		// Stop at the first feasibile solution
		//cplex.setParam(IloCplex::IntSolLim, 1);

		// Stop after reaching the time limit 2 hours
		cplex.setParam(IloCplex::TiLim, 7200);

		// Set maximum number of threads 
		cplex.setParam(IloCplex::Threads, 4);

		if (!cplex.solve()) {
			env.error() << "Failed to optimize LP" << endl;
			env.out() << "Solution status = " << cplex.getStatus() << endl;
			throw(-1);
		}

		env.out() << "Solution status = " << cplex.getStatus() << endl;
		env.out() << "Solution value  = " << cplex.getObjValue() << endl;


		//-----------------------------------------------------------------------------
		// SAVE EVERYTHING
		//-----------------------------------------------------------------------------

		// Save objective function output
		if (mytarget == MINIMIZE_UPDATE_INT)
			MILP_out = cplex.getValue(-OBJ);
		else
			MILP_out = cplex.getValue(OBJ);

		ofstream results;
		results.open("results.txt");

		results << "Task" << "\t" << "Period" << "\t" << "Offs";

		for (int p = 1; p < NUMBER_OF_PATHS; p++) {
			results << "\t" << "VMiss" << p;
			results << "\t" << "Vhit" << p; 
			results << "\t" << "RedHs" << p;
			results << "\t" << "Miss" << p;
			results << "\t" << "IncHs" << p;
			
		}
		results << "\t" << "VMiss" << NUMBER_OF_PATHS;
		results << "\t" << "Vhit" << NUMBER_OF_PATHS;
		results << endl;

		for (int t = 0; t < NUMBER_OF_TASKS_IN_CHAIN; t++) {

			results << taskchain.at(t).id << "\t" << taskchain.at(t).period << "\t";
			results << (cplex.getValue(OFFS[t]));
			
			for (int p = 0; p < NUMBER_OF_PATHS - 1; p++) {

				int Mv = round(cplex.getValue(MISSAFTEREFFECTIVE[t][p]));
				results << "\t" << Mv;
				int V = round(cplex.getValue(EFFECTIVEJOB[t][p]));
				results << "\t" << V ;
				int Hr = round(cplex.getValue(REDUNDHITS[t][p]));
				results << "\t" << Hr;
				int NM = round(cplex.getValue(MISSWNEWINPUT[t][p]));
				results << "\t" << NM;
				int Hv = round(cplex.getValue(VOIDHITS[t][p]));
				results << "\t" << Hv;

			}
			int Mv = round(cplex.getValue(MISSAFTEREFFECTIVE[t][NUMBER_OF_PATHS - 1]));
			results << "\t" << Mv;
			int V = round(cplex.getValue(EFFECTIVEJOB[t][NUMBER_OF_PATHS - 1]));
			results << "\t"  << V;

			results << endl;
		}
		results.close();

	} // End of try

	//-----------------------------------------------------------------------------
	// EXCEPTIONS
	//-----------------------------------------------------------------------------

	catch (IloAlgorithm::CannotExtractException &e) {
		IloExtractableArray &failed = e.getExtractables();
		std::cerr << "Failed to extract:" << std::endl;
		for (IloInt i = 0; i < failed.getSize(); ++i)
			std::cerr << "\t" << failed[i] << std::endl;
	}
	catch (IloException& e) {
		cerr << "Concert exception caught: " << e << endl;
	}

	env.end();

	return MILP_out;
} 
