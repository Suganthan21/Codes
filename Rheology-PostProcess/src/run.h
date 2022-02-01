#ifndef RUN_H
#define RUN_H

#include "pointers.h"
#include "fix_swap_mc.h"

class Run:Pointers{
	public:
		Run(Main*);
		~Run();
		int nsteps, ntimestep;
		double neigh_cutoff;
		void Simulation();
		void OneStep();

	private:
		class FixSwapMC *swapmc;
		void TimingStatistics();
		void ThermoOutput();
		void Init();
		void Test();
};

#endif
