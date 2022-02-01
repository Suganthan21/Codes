#ifndef PAIR_H
#define PAIR_H

#include "pointers.h"
#include<vector>

class Pair:Pointers
{
	public:
		Pair(Main*);
		~Pair();

		double pair_time;

		//LJ
		float epsilon[3][3];
		float sigma[3][3];
		float rc2[3][3];
		double lj1[3][3];
		double lj2[3][3];
		double lj3[3][3];
		double lj4[3][3];
		double offset[3][3];

		//DPD
		float gamma_dpd;
		float cutoff_dpd;
		float cutoffsq_dpd;

		//LD
		float gamma_ld;
		
		void Compute(const std::vector<int>&);
		void PerAtomProps(int,int,double,double,double,double, double);
		void LJ(const std::vector<int>&);
		void LJPolydisperse(const std::vector<int>&);
		void DragDPD(const std::vector<int>&);
		void DragLD(const std::vector<int>&);
		void TotalForce(const std::vector<int>&);

};

#endif

