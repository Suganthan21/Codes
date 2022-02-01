#ifndef PARTICLE_H
#define PARTICLE_H

#include "pointers.h"
#include<vector>

class Particle:Pointers
{
	public:
		Particle(Main*,int,int);
		~Particle();

		double ***r,***v,***f,***lj,***fdrag, ***press;

		double **F, **LJ, **FDRAG;
		double **pe, **ke, **sxx, **syy, **sxy, **sxz, **syz, **szz;
		double ***pol;

		int **ix, **iy, **iz;
		int *id, *type, *mol;
		double *rad, *mass;

		void Unwrap(const std::vector<int>&);
		double MaxRadius();
};

#endif

