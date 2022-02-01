#ifndef UTILS_H
#define UTILS_H

#include "pointers.h"
#include<vector>
#include<string>

class Utils:Pointers
{
	public:
// 		Utils(int, double, double);
		Utils(Main*);
		~Utils();
		std::vector<int> Sec2HHMMSS(double);
		std::string ToString(double, int);
		void ChangeLogName(std::string, std::string);
};

#endif

