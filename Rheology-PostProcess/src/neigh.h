#ifndef NEIGH_H
#define NEIGH_H

#include "pointers.h"
// #include<map>

class Neigh:Pointers
{
	public:
		Neigh(Main*,int);
		~Neigh();

		int nbuild;
		double neigh_time;

// 		std::vector<std::multimap<double,int>> mapnlist;
		int **nlist;
		int *nlistsize;

		void Build(double);
		void VerletList(double);
};

#endif
