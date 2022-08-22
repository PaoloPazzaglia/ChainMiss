#include "milp_WHchain.h"
#include <random>
#include <iostream>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <chrono>

using namespace std;

#define INPUT_FILE true

#define NUM_TESTS 67

#define NUM_TASKS 25
#define PERIODS_IN_BUCKET false

// Select k for (m,k) constraint
#define CHOSEN_K 10		

#define MYPERRULE 2
#define MYMKTASKS 3


int main()
{
	// Initialize random engine
	std::random_device rd;     // only used once to initialise (seed) engine
	std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)


	// Extract period rule (oversampling, undersampling, mixed)
	PeriodsRule myPerRule = static_cast<PeriodsRule>(MYPERRULE);

	// Extract weakly hard task rule (top, mid, bottom, mixed);
	mkTasks mymkTasks = static_cast<mkTasks>(MYMKTASKS);


	// Taskchain initialization
	std::vector<Task> taskchain;
	std::vector<WHconstr> setofmk;


	if (INPUT_FILE) { // Input file

		// Iterate for each metrics
		for (int i = 0; i < 4; i++) {

			// Create and open output file
			ofstream all_out;
			string nameoutputfile;

			switch (i) {
			case 0:
				nameoutputfile = "data_maxIOL";
				break;
			case 1:
				nameoutputfile = "data_maxDA";
				break;
			case 2:
				nameoutputfile = "data_maxUI";
				break;
			case 3:
				nameoutputfile = "data_minUI";
				break;
			default:
				return 1;
				break;
			}

			all_out.open(nameoutputfile + ".csv");

			// For each mk values for chosen task
			for (int m = 0; m <= 30; m++) {

				all_out << m;

				// For each chain
				for (int c = 1; c <= 5; c++) {

					// Read input file
					std::string namefile = "perceptin";
					ifstream infile(namefile + std::to_string(c) + ".txt");

					if (!infile.is_open()) {
						exit(EXIT_FAILURE);
					}

					// Initialize appo variables
					int id, period, deadline;
					bool mktask;
					vector<int> mktaskid;

					taskchain.clear();
					setofmk.clear();

					// Start reading data from input file
					string str;
					getline(infile, str); // skip the first line

					while (infile >> id >> period >> deadline >> mktask) {

						// Create task chain
						Task t;
						t.id = id;
						t.period = period;
						t.deadline = deadline;
						taskchain.push_back(t);

						// Create mk model
						MKconstr mkc;
						WHconstr whc;

						// Initialize everything hard deadline
						mkc.m = 0;
						mkc.k = 1;

						// Store id of task with weakly-hard behavior
						if (mktask)
							mktaskid.push_back(id);

						// Add mk parameters to list
						whc.taskid = id;
						whc.mconsec = mkc.m;
						whc.mk.push_back(mkc);
						setofmk.push_back(whc);

					} // end reading file

					  					  
					for (int j = 0; j < mktaskid.size(); j++) {
						// Update mk values
						setofmk.at(mktaskid.at(j)).mconsec = m; // + j //for test VIII.C
						setofmk.at(mktaskid.at(j)).mk.at(0).k = 50;
						setofmk.at(mktaskid.at(j)).mk.at(0).m = m;  // + j //for test VIII.C
					}
					  

					// Perform the test
					OptTarget mytarget = static_cast<OptTarget>(i);
					double output = MILP_WH_K(taskchain, setofmk, mytarget);
					all_out << ',' << output;

				} // end test of chain

				all_out << endl;

			} // end mk value 

			  // Close output file
			all_out.close();
		} // end metric test
	}


	//---------------------------------------------------------

	else { // Pseudo-random task chain

		// Initialize structures to save the results
		vector<double> max_m, min_m, av_m;
		vector<double> tmax, tmin, tav;

		for (int i = 0; i < 4; i++) {
			max_m.push_back(0);
			min_m.push_back(INT_MAX);
			av_m.push_back(0);
		}

		for (int i = 0; i < 4; i++) {
			tmax.push_back(0);
			tmin.push_back(INT_MAX);
			tav.push_back(0);
		}

		// Start tests
		for (int s = 0; s < NUM_TESTS; s++) {


			// Building the task set

			//-------------------------------------------------------------
			// Store chosen periods
			vector<int> chosen_periods;

			// Build random set of periods
			if (PERIODS_IN_BUCKET) {
				// Periods to choose from (automotive standard, ms)
				vector<int> period_bucket = { 1, 2, 5, 10, 20, 25, 50, 100 };

				int TOT_NUM_PERIODS = period_bucket.size();
				std::uniform_int_distribution<int> uni(0, TOT_NUM_PERIODS - 1); // guaranteed unbiased

				for (int i = 0; i < NUM_TASKS; i++) {
					// Choose random period id					
					int appo_id = uni(rng);
					chosen_periods.push_back(period_bucket.at(appo_id));
				}
			}

			else { // Random periods

				// Maximum value for period (ms)
				int MAX_T = 100;
				std::uniform_int_distribution<int> uni(1, MAX_T); // guaranteed unbiased

				for (int i = 0; i < NUM_TASKS; i++) {
					// Choose random period

					int rand_period = uni(rng);
					chosen_periods.push_back(rand_period);
				}
			}


			// Reorder periods if required
			switch (myPerRule) {
			case UN:
				sort(chosen_periods.begin(), chosen_periods.end());
				break;
			case OV:
				sort(chosen_periods.rbegin(), chosen_periods.rend());
				break;
			}

			//-------------------------------------------------------------
			// Create tasks in chain

			taskchain.clear();

			for (int i = 0; i < NUM_TASKS; i++) {
				Task t;
				t.id = i;
				t.period = chosen_periods.at(i);
				t.deadline = t.period;
				taskchain.push_back(t);
			}

			//-------------------------------------------------------------
			// Assign (m,k) values

			setofmk.clear();

			std::uniform_int_distribution<int> uni(3, CHOSEN_K); // guaranteed unbiased

			for (int i = 0; i < NUM_TASKS; i++) {
				MKconstr mkc;
				WHconstr whc;

				mkc.k = uni(rng);
				mkc.m = 0;

				// Select m depending on the chosen value K and assignment rule
				switch (mymkTasks) {
				case TOP:
					if (i < NUM_TASKS / 3)
						mkc.m = floor(static_cast<double>(CHOSEN_K) / 3.0);
					break;
				case MID:
					if (i >= NUM_TASKS / 3 && i < 2 * NUM_TASKS / 3)
						mkc.m = floor(static_cast<double>(CHOSEN_K) / 3.0);
					break;
				case BOT:
					if (i >= 2 * NUM_TASKS / 3)
						mkc.m = floor(static_cast<double>(CHOSEN_K) / 3.0);
					break;
				case MIX:
					mkc.m = floor(static_cast<double>(CHOSEN_K) / 3.0);
					break;
				}

				// Consecutive deadline misses
				whc.mconsec = mkc.m;

				// Add weakly-hard constraint to set
				whc.mk.push_back(mkc);
				setofmk.push_back(whc);
			}

			// Perform the test for all the metrics
			for (int i = 0; i < 1; i++) {
				OptTarget mytarget = static_cast<OptTarget>(i);

				cout << "********TEST NUMBER " << s << endl << endl;

				auto start_time = chrono::steady_clock::now();
				double output = MILP_WH_K(taskchain, setofmk, mytarget);
				auto end_time = chrono::steady_clock::now();

				/*
				all_out << output << ',';

				if (output > max_m.at(i))
					max_m.at(i) = output;
				if (output < min_m.at(i))
					min_m.at(i) = output;
				av_m.at(i) += (double)(output / NUM_TESTS);
				*/

				double runtime = chrono::duration_cast<chrono::seconds>(end_time - start_time).count();

				ofstream time_out;
				string nameout = "exec_time_" + std::to_string(NUM_TASKS) + ".csv";
				time_out.open(nameout, std::ios_base::app);

				time_out << runtime << ',' << endl;

				time_out.close();

			}
		}
		

		ofstream summary_out;
		summary_out.open("summary_data.csv");

		summary_out << NUM_TASKS;
		for (int i = 0; i < 4; i++) {
			summary_out << max_m.at(i) << ',';
		}
		summary_out << endl;
		summary_out << NUM_TASKS;
		for (int i = 0; i < 4; i++) {
			summary_out << min_m.at(i) << ',';
		}
		summary_out << endl;
		summary_out << NUM_TASKS;
		for (int i = 0; i < 4; i++) {
			summary_out << av_m.at(i) << ',';
		}
		summary_out << endl;

		summary_out.close();


		

	}

	
	std::fflush(stdout);
#ifdef _WIN32
	std::system("pause");
#endif
	return 0;
}
