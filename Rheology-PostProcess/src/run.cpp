#include "run.h"
#include "filehandle.h"
#include "memory.h"
#include "simbox.h"
#include "pair.h"
#include "neigh.h"
#include "group.h"
#include "particle.h"
#include "thermo.h"
#include "fix_swap_mc.h"
#include "utils.h"

using std::vector;
using std::string;
using std::to_string;

Run::Run(Main *ptr):Pointers(ptr){
	ntimestep = 0;
	swapmc = new FixSwapMC(ptr);
// 	neigh_cutoff = (2*atom->MaxRadius())*1.122462048309373 + 0.3; //(rad[i] + rad[j])*2^(1/6) + skin
	neigh_cutoff = 2.8;
}

Run::~Run(){
	delete swapmc;
}

void Run::Init(){
	filehandle->writedump = 1;
	string str = "_T"+utils->ToString(swapmc->T,1)+"-p"+utils->ToString(swapmc->pswap,1);
	string logname = "log"+str+".log";
// 	fclose(log);
// 	utils->ChangeLogName("log.dat",logname);
// 	log = fopen(logname.c_str(),"a");
	filehandle->odump.push_back("dump-data/config"+str+".dump");
	filehandle->OpenClose("open",0);
	neigh->Build(neigh_cutoff);
	pair->Compute(group->all);
	filehandle->WriteDump(group->all,ntimestep);
	fprintf(log,"\nThermodynamic Output\nStep PotEng MCTemp Press\n");
	ThermoOutput();
}

void Run::ThermoOutput(){
	pair->Compute(group->all);
	thermo->PE(group->all);
	thermo->PRESS(group->all);

	fprintf(log,"%d %.8f %g %.8f\n", ntimestep, thermo->pe[n], swapmc->T, thermo->press[n][0]);	
}

void Run::Simulation(){
	nsteps = 11000;
	Init();

	for(int i=1; i<=nsteps; i++){
		ntimestep++;
		if(ntimestep%15 == 0)
			neigh->Build(neigh_cutoff);

		swapmc->OneStep();

		if(ntimestep%1000 == 0)
			filehandle->WriteDump(group->all,ntimestep);

		if(ntimestep%10 == 0)
			ThermoOutput();

// 		if(ntimestep%10 == 0)
// 			Test();
// 	 	thermo->Reset();
	}
	filehandle->OpenClose("close",0);
	swapmc->Statistics();
	TimingStatistics();
}

void Run::TimingStatistics(){
	vector<int> pair_time = utils->Sec2HHMMSS(pair->pair_time);
	vector<int> neigh_time = utils->Sec2HHMMSS(neigh->neigh_time);

	fprintf(log,"Timing breakdown:\n"
				"\tPair potential = %d:%d:%d\n",
				pair_time[0],pair_time[1], pair_time[2]);
	fprintf(log,"\t%d Verlet Lists built in time = %d:%d:%d\n",
				neigh->nbuild,neigh_time[0],neigh_time[1], neigh_time[2]);

}

void Run::Test(){
	double ***r = atom->r;
	double ***f = atom->f;
	fprintf(log,"\n\n");
	for(int i=0; i<10; i++)
		fprintf(log,"%g %g %g %g %g %g\n", r[0][i][n], r[1][i][n], r[2][i][n], f[0][i][n], f[1][i][n], f[2][i][n]);
	fprintf(log,"\n\n");
}
