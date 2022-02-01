#ifndef POINTERS_H
#define POINTERS_H

#include "main.h"

class Pointers
{
	public:
		Pointers(Main* ptr):
		n(ptr->n),
		log(ptr->log),
		filehandle(ptr->filehandle),
		memory(ptr->memory),
		atom(ptr->atom),
		simbox(ptr->simbox),
		group(ptr->group),
		neigh(ptr->neigh),
		pair(ptr->pair),
		thermo(ptr->thermo),
		c_spatial(ptr->c_spatial),
		c_temporal(ptr->c_temporal),
		run(ptr->run),
		utils(ptr->utils){}
	
		int &n;
		FILE* &log;
		FileHandle* &filehandle;
		Memory* &memory;
		Particle* &atom;
		SimBox* &simbox;
		Group* &group;
		Neigh* &neigh;
		Pair* &pair;
		Thermo* &thermo;
		Compute_spatial* &c_spatial;
		Compute_temporal* &c_temporal;
		Run* &run;
		Utils* &utils;
};

#endif
