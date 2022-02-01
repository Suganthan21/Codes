#ifndef MAIN_H
#define MAIN_H

#include<cstdio>
//singleton class - instantiate only once

class Main
{
	public:
// 		static Main* Instantiate();
// 		void Destruct();
		Main(){}
		~Main();
		
		int n; //current ntimes - post process
		FILE *log; //logfile
		class FileHandle* filehandle;
		class Memory* memory;
		class Particle* atom;
		class SimBox* simbox;
		class Group* group;
		class Neigh* neigh;
		class Pair* pair;
		class Thermo* thermo;
		class Compute_spatial* c_spatial;
		class Compute_temporal* c_temporal;
		class Run* run;
		class Utils* utils;

		void InitProcess();
		void PostProcess();
		void Simulation();
		void TimingStatistics(double);
// 	private:
// 		static Main* instance;
// 		Main(){}
// 		Main(const Main& copy);
// 		Main& operator=(const Main& copy);
};

#endif
