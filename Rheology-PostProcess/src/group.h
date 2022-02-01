#ifndef GROUP_H
#define GROUP_H

#include "pointers.h"
#include<string>
#include<vector>

class Group:Pointers
{
	public:
		Group(Main*);
		~Group(){}
		int nall, nactive, npassive;
		double tauP, f0, pc;
		std::vector<int> all;
		std::vector<int> nselect;
		std::vector<int> active;
		std::vector<int> passive;
		std::vector<int> neigh;

		void SelectAtoms(std::string);
};

#endif
