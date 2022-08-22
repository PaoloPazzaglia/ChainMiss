#ifndef MILP_DATA_H__
#define MILP_DATA_H__

#include <vector>
#include <string>
#include <cstdint>
#include <iostream>

struct Task {
	int id;
	int deadline;
	int period;
	std::string name;
	int core_id = 0;
};

struct MKconstr {
	int m;
	int k;
};

struct WHconstr {
	int taskid;
	int mconsec;
	std::vector<MKconstr> mk;
};

enum OptTarget {
	MAXIMIZE_LATENCY = 0,
	MAXIMIZE_DATAAGE = 1,
	MAXIMIZE_UPDATE_INT = 2,
	MINIMIZE_UPDATE_INT = 3
};

enum DMstrat {
	KILL = 0,
	SKIP_NEXT = 1
};

enum PeriodsRule {
	OV = 0,
	UN = 1,
	MX = 2
};

enum mkTasks {
	TOP = 0,
	MID = 1,
	BOT = 2,
	MIX = 3
};

#endif
