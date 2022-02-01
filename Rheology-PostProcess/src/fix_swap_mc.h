#ifndef FIX_SWAP_MC_H
#define FIX_SWAP_MC_H

#include "pointers.h"
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>

class FixSwapMC:Pointers
{
	public:
		FixSwapMC(Main*);
		~FixSwapMC();
		void OneStep();
		void Statistics();
		double T, pswap; // temperature, probability to swap particles

	private:
		gsl_rng *r1, *r2, *r3, *r4;
		double delta;  // translational distance to move
		double Ti, Tf, deltaT;
		double dispsq; // for efficiency check
		int nswap_attempt, nswap, ntrans_attempt, ntrans;     // attempted, accepted swap moves
		int config_change; // has config changed in last MC step

		void SwapRadius(int, int);
		void Efficiency();
		int Metropolis(double, double);

};

#endif
