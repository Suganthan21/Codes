#ifndef THERMO_H
#define THERMO_H

#include "pointers.h"
#include<vector>

class Thermo:Pointers
{
	public:
		Thermo(Main*, int);
		~Thermo();
		
		int *step;
		double *pe, *etot, *sxy, *fmax, *fnorm, *pol;
		double **ke, **temp, **press;

		// to read a per atom property from data file instead of dump. 
		// length of array must be specified in constructor.
		double **custom; 
		
		int keflag, peflag, pressflag;
		void PE(const std::vector<int>&);
		void KE(const std::vector<int>&);
		void ETOT(const std::vector<int>&);
		void PRESS(const std::vector<int>&);

		void Reset();
};

#endif

