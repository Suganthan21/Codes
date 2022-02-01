#ifndef COMPUTE_TEMPORAL_H
#define COMPUTE_TEMPORAL_H

using namespace std;

#include "pointers.h"

#include<vector>
#include<fstream>

class Compute_temporal:Pointers
{
	public:
		Compute_temporal(Main*);
		~Compute_temporal();
		
		void TimeOriginAvg(const std::vector<int>&);
		void MSD();
	
	private:
		ifstream* CalibFile;
		string calibratefile;
		int ntimes;
		std::vector<int> itstep;
		//OverlapFunc and MSD
		std::vector<int> tdiff;
		double **q, **msd; 	// First index is tdiff position, q[][1] is q value, q[][2] is counter
		double **msd_theta, **msd_phi;
		int tavgrows;
		static int tdiffcount;

		void CalibrateTdiff();
		void WriteTimeOriginAvg(int);
		void PositionMSD(int, int, const std::vector<int>&);
		void AngularMSD(int, int, const std::vector<int>&);
		void EmptyArrays();
		
		FILE* fp1;

};
#endif
