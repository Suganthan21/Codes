#include "main.h"
#include "filehandle.h"
#include "memory.h"
#include "particle.h"
#include "simbox.h"
#include "group.h"
#include "pair.h"
#include "neigh.h"
#include "thermo.h"
#include "compute_spatial.h"
#include "compute_temporal.h"
#include "run.h"
#include "utils.h"

using std::vector;

// Main* Main::Instantiate(){
// 	if(instance==nullptr)
// 		instance=new Main; 
// 	return instance;
// }

// Main* Main::instance=nullptr;

Main::~Main(){
// 	Destruct();
// }
// 
// void Main::Destruct(){
// 	delete instance;
// 	instance=nullptr;
	delete memory;
	memory = nullptr;
	delete simbox;
	simbox=nullptr;
	delete group;
	group=nullptr;
	delete filehandle;
	filehandle=nullptr;
	delete atom;
	atom=nullptr;
	delete pair;
	pair=nullptr;
	delete thermo;
	thermo= nullptr;
	delete c_spatial;
	c_spatial = nullptr;
	delete c_temporal;
	c_temporal = nullptr;
	delete run;
	run = nullptr;
	delete utils;
	utils = nullptr;
}


void Main::Simulation(){
	log = stderr;
// 	log = fopen("log.dat","w");
	memory = new Memory(this);
	filehandle = new FileHandle(this);
	filehandle->ntimes = 1; //only current timestep per atom info is stored
	filehandle->filetype = "dump"; //dump or dump,data
	simbox = new SimBox(this);
	simbox->dt = 0.001;
	simbox->dim = 3;
// 	string fname="samples/VF0.70-N864-S1-CG-vseed11111-Minimize.int";
	string fname="samples/VF0.70-N864-Equilibrate-T5.int";
	fprintf(log,"INPUT FILE = %s\n", fname.c_str());
	simbox->BoxData(fname,"init");

	atom = new Particle(this,filehandle->ntimes,simbox->natoms);

	simbox->Initialise(fname,"Atoms,Radius,Masses");
	neigh = new Neigh(this,simbox->natoms);
	pair = new Pair(this);

	string groups = "all";
	group = new Group(this);
	group  ->  SelectAtoms(groups);

	thermo = new Thermo(this,filehandle->ntimes);
	c_spatial = new Compute_spatial(this,filehandle->ntimes);
	c_temporal = new Compute_temporal(this);
	utils = new Utils(this);
	
	fprintf(log,"\nSetup complete\n");
	memory->Allocation();
	fprintf(log,"\nStarting simulation...\n");

	run = new Run(this);
	run->Simulation();
}

void Main::InitProcess(){
	log = stderr;
	memory = new Memory(this);
	filehandle = new FileHandle(this);

	filehandle->filetype = "dump"; //dump or dump,data
	string ft = filehandle->filetype;

	simbox = new SimBox(this);
	group = new Group(this);
	
	filehandle->Init();
	simbox->dt = 0.001;
	simbox->dim=3;
	simbox->natomtypes = 1;
	
	if(ft == "dump" || ft == "dump,data" ){
		simbox->BoxData(filehandle->idump[0], ft);
		atom = new Particle(this,filehandle->ntimes,simbox->natoms);
		neigh = new Neigh(this,simbox->natoms);
		pair = new Pair(this);
		string groups = "all";
		group  ->  SelectAtoms(groups);
	}

	thermo = new Thermo(this,filehandle->ntimes);
	c_spatial = new Compute_spatial(this,filehandle->ntimes);
	c_temporal = new Compute_temporal(this);
	utils = new Utils(this);

	fprintf(log,"\nSetup complete\n");
	memory->Allocation();
	fprintf(log,"\nStarting to process...\n");
}


void Main::PostProcess(){
	InitProcess();

	int len = filehandle->idump.size();
	len = 1;

	for(int i=0; i < len; i++){
		n=0;
		filehandle ->  OpenClose("open",i);
// 		simbox->Properties(filehandle->idump[i]);

// 		filehandle	-> StoreData();
// 		c_temporal	-> MSD();

		double time;
		for(auto step:filehandle->ststep){
			fprintf(log,"Timestep = %s\n", step.c_str());
			time = double(stod(step))*simbox->dt;
			filehandle   ->   StoreDump(step);
			filehandle	 ->	  Dump2Dat(step,group->all);

// 			atom		 ->	  Unwrap(group->all);
// 			c_temporal	 ->   TimeOriginAvg(group->all);
// 			c_spatial	 ->   Correlation(group->all);
// 			filehandle   ->   WriteDump(step,group->all);
			n++;
		}

		filehandle->OpenClose("close",i);
	}
}

void Main::TimingStatistics(double time_sec){
	
	std::vector<int> total_time = utils->Sec2HHMMSS(time_sec);

	fprintf(log,"\tTotal Time = %d:%d:%d\n",
				total_time[0],total_time[1], total_time[2]);
}

int main(){
// 	Main *main = Main::Instantiate();
	Main *main = new Main();
	main->n=0;

	clock_t clk_start;
	clk_start = clock();

	try{
// 		main->PostProcess();
		main->Simulation();
	}
	catch(...){
		fprintf(main->log,"ENDING EXECUTION\n");
	}
	//will result in error if class is singleton

	double time_sec = double(clock() - clk_start)/(CLOCKS_PER_SEC);
	main->TimingStatistics(time_sec);
	
	fclose(main->log);
	delete main;

	return 0;
}
