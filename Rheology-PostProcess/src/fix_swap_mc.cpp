#include<cmath>
#include "fix_swap_mc.h"
#include "particle.h"
#include "simbox.h"
#include "group.h"
#include "pair.h"
#include "run.h"
#include "thermo.h"

FixSwapMC::FixSwapMC(Main* ptr):Pointers(ptr) {
	r1 =  gsl_rng_alloc(gsl_rng_mt19937); // to check pswap
	r2 =  gsl_rng_alloc(gsl_rng_mt19937); // to pick random particle
	r3 =  gsl_rng_alloc(gsl_rng_mt19937); // to displace particles
	r4 =  gsl_rng_alloc(gsl_rng_mt19937); // Metropolis

	// seed
	gsl_rng_set(r1,12322);
	gsl_rng_set(r2,14595);
	gsl_rng_set(r3,579826);
	gsl_rng_set(r4,98687);
	

	nswap_attempt = 0;
	ntrans_attempt = 0;
	nswap = 0;
	ntrans = 0;
	config_change = 1;
	dispsq = 0.0;

	pswap = 0.0;
	Ti = 5.0;
	Tf = 0.01;
	deltaT = 0.02;
	T = Ti;
	delta = 0.45*atom->MaxRadius();

	if(simbox->shearflag){
		fprintf(log,"ERROR: FIXSWAPMC - SIMULATION BOX MUST BE ORTHOGONAL");
		throw std::exception();
	}
	fprintf(log,"Fix Swap MC setttings:\n"
				"\tPswap = %g, Ti = %g, Tf = %g\n", pswap, Ti, Tf);
	fprintf(log,"\tTranslate particles by length = %g\n\n", delta);
}

FixSwapMC::~FixSwapMC(){	
	gsl_rng_free(r1);
	gsl_rng_free(r2);
	gsl_rng_free(r3);
	gsl_rng_free(r4);
}

void FixSwapMC::OneStep(){
	if(run->ntimestep%10000 == 0)
		T -= deltaT;
// 	T = Ti - (Ti-Tf)*(run->ntimestep/run->nsteps);
	double randnum = gsl_rng_uniform(r1);
	double Ei, Ef;
	
	// Swap Particles
	if(randnum < pswap){
		nswap_attempt++;
		int id1, id2;
		id1 = gsl_rng_uniform_int(r2, simbox->natoms); // pick a particle btwn [0,natoms-1]
		id2 = gsl_rng_uniform_int(r2, simbox->natoms);

		if(config_change){
			pair->Compute(group->all);
			thermo->PE(group->all);
		}

		Ei = thermo->pe[n]*simbox->natoms;
		
		SwapRadius(id1, id2);

		pair->Compute(group->all);
		thermo->PE(group->all);
		Ef = thermo->pe[n]*simbox->natoms;
		
		config_change = Metropolis(Ei,Ef);

		if(config_change)
			nswap++;
		else{
			SwapRadius(id1,id2);  //swap back
			thermo->pe[n] = Ei/simbox->natoms;
		}
	}
	// Translational move
	else{	
		ntrans_attempt++;
		int id = gsl_rng_uniform_int(r2, simbox->natoms);
		if(config_change){
			pair->Compute(group->all);
			thermo->PE(group->all);
		}
		Ei = thermo->pe[n]*simbox->natoms;

		double randx = gsl_rng_uniform(r3)-0.5;
		double randy = gsl_rng_uniform(r3)-0.5;
		double randz = gsl_rng_uniform(r3)-0.5;

		double ***r = atom->r;
		r[0][id][n] += randx*delta;
		r[1][id][n] += randy*delta;
		r[2][id][n] += randz*delta;
		
		pair->Compute(group->all);
		thermo->PE(group->all);
		Ef = thermo->pe[n]*simbox->natoms;

		config_change = Metropolis(Ei,Ef);
		
		if(config_change){
			ntrans++;
			dispsq += pow(randx*delta,2) + pow(randy*delta,2) + pow(randz*delta, 2);
		}
		else{	
			r[0][id][n] -= randx*delta;
			r[1][id][n] -= randy*delta;
			r[2][id][n] -= randz*delta;
			thermo->pe[n] = Ei/simbox->natoms;
		}
	}
}

void FixSwapMC::SwapRadius(int id1, int id2){
	double *rad = atom->rad;
	double tmp = rad[id1];
	rad[id1] = rad[id2];
	rad[id2] = tmp;
}

int FixSwapMC::Metropolis(double Ei, double Ef){
	double randnum = gsl_rng_uniform(r4);
	if( Ef < Ei || randnum < exp(-(Ef-Ei)/ T)  )
		return 1;
	else
		return 0;
}

void FixSwapMC::Statistics(){
	fprintf(log,"\nMonte carlo statistics:\n"
				"\tNumber of attempted swap moves = %d\n"
				"\tNumber of successful swap moves = %d\n"
				"\tAcceptance ratio swap moves = %g\n"
				"\tNumber of attempted translation moves = %d\n"
				"\tNumber of successful translation moves = %d\n"
				"\tAcceptance ratio translation moves = %g\n\n",
				nswap_attempt, nswap, double(nswap)/double(nswap_attempt),
				ntrans_attempt, ntrans, double(ntrans)/double(ntrans_attempt) );
}
