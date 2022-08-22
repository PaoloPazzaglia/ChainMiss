#include "milp_WHchain.h"
#include <sstream>

std::string  convert_to_string(const int Number) {
	std::stringstream convert;

	convert << Number;

	std::string Result = "_" + convert.str();
	return Result;
}