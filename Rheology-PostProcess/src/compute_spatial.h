#ifndef COMPUTE_SPATIAL_H
#define COMPUTE_SPATIAL_H

#include "pointers.h"
#include<fstream>
#include<vector>
#include<map>

using namespace std;

class Compute_spatial:Pointers
{
	public:
		Compute_spatial(Main*, int);
		~Compute_spatial();

		void InitSquares(const std::vector<int>&);
		void FilledSquares(const std::vector<int>&);
		void Correlation(const std::vector<int>&);

	private:
		class Square{
			public:
			double xmid,ymid, xlo, xhi, ylo, yhi, initxlo, initxhi, initxmid, len;
			int fill;
		};
		int nsquare;
		std::vector<Square> sq;
		std::vector<multimap <double,double>> map_coarse_corr;
		map<int,double> adder, counter;

		bool flag_spacecorr,flag_filledsquares;
		
		void CoarsenCorrelation(const multimap<double,double>&);
		void FinalWrite();

		FILE* fp1;
};

#endif
