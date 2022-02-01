#include<cmath>
#include "neigh.h"
#include "memory.h"
#include "particle.h"
#include "group.h"
#include "simbox.h"

Neigh::Neigh(Main* ptr, int natoms):Pointers(ptr){
	neigh_time = 0.0;
	nbuild = 0;
	int size = 500;
	nlist = memory->Create2DArray<int>(natoms,size); //second entry - max no.of neighbours one atom can have
	nlistsize = memory->Create1DArray<int>(natoms);
// 	mapnlist.resize(natoms);
	fprintf(log,"\nNeigh list settings: One atom can have %d neighbours\n", size);
}

Neigh::~Neigh(){
	memory->Delete2DArray(nlist);
	memory->Delete1DArray(nlistsize);
}

void Neigh::Build(double cut){
	clock_t clk_start = clock();
	nbuild++;

// 	mapnlist.clear();
	for(int i=0;i<simbox->natoms;i++)
		nlistsize[i] = 0;

	VerletList(cut);
	neigh_time += double(clock() - clk_start)/(CLOCKS_PER_SEC);
}

void Neigh::VerletList(double cut){

	double L = simbox->L;
	double cutsq = cut*cut;
	double ***r = atom->r;
	double r2, delx, dely,delz;

	int i,j,pos;
	for(i=0; i<simbox->natoms; i++){
		for(j=i+1; j<simbox->natoms; j++){

			delx = r[0][i][n] - r[0][j][n]; delx -= round(delx/L)*L; 
			dely = r[1][i][n] - r[1][j][n];	dely -= round(dely/L)*L;
			delz = r[2][i][n] - r[2][j][n];	delz -= round(delz/L)*L;

			r2 = delx*delx + dely*dely + delz*delz;

			if(r2 < cutsq){
				pos = nlistsize[i]; 	
				nlist[i][pos] = j;
				nlistsize[i]++;
// 				mapnlist[i].insert({r2,j});
			}			

		}
	}
		
	
}

