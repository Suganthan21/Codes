#include<cmath>
#include "particle.h"
#include "memory.h"
#include "simbox.h"

Particle::Particle(Main* ptr,int ntimes,int natoms):Pointers(ptr){
	r = nullptr;
	v = nullptr;
	f = nullptr;
	lj = nullptr;
	fdrag = nullptr;
	press = nullptr;

	F = nullptr;
	LJ = nullptr;
	FDRAG = nullptr;
	pe = nullptr;
	ke = nullptr;
	sxx = nullptr;
	syy = nullptr;
	sxy = nullptr;
	szz = nullptr;
	sxz = nullptr;
	syz = nullptr;

	pol = nullptr;

	ix = nullptr;
	iy = nullptr;
	iz = nullptr;
	id = nullptr;
	type = nullptr;
	mol = nullptr;
	rad = nullptr;
	mass = nullptr;

	r = memory->Create3DArray<double>(3,natoms,ntimes);
	v = memory->Create3DArray<double>(3,natoms,ntimes);
	f = memory->Create3DArray<double>(3,natoms,ntimes);
	lj = memory->Create3DArray<double>(3,natoms,ntimes);
	fdrag = memory->Create3DArray<double>(3,natoms,ntimes);
	press = memory->Create3DArray<double>(6,natoms,ntimes);

	F = memory->Create2DArray<double>(natoms,ntimes);
	LJ = memory->Create2DArray<double>(natoms,ntimes);
	FDRAG = memory->Create2DArray<double>(natoms,ntimes);
	pe = memory->Create2DArray<double>(natoms,ntimes);
	ke = memory->Create2DArray<double>(natoms,ntimes);
	sxx = memory->Create2DArray<double>(natoms,ntimes);
	syy = memory->Create2DArray<double>(natoms,ntimes);
	sxy = memory->Create2DArray<double>(natoms,ntimes);
	sxz = memory->Create2DArray<double>(natoms,ntimes);
	syz = memory->Create2DArray<double>(natoms,ntimes);
	szz = memory->Create2DArray<double>(natoms,ntimes);

	pol = memory->Create3DArray<double>(7,natoms,ntimes);

	ix = memory->Create2DArray<int>(natoms,ntimes);
	iy = memory->Create2DArray<int>(natoms,ntimes);
	iz = memory->Create2DArray<int>(natoms,ntimes);
	id = memory->Create1DArray<int>(natoms);
	type = memory->Create1DArray<int>(natoms);
	mol = memory->Create1DArray<int>(natoms);

	rad = memory->Create1DArray<double>(natoms);
	mass = memory->Create1DArray<double>(natoms);

	fprintf(log,"\nMemory allocated for particle properties \n");

}

Particle::~Particle(){

	memory->Delete3DArray(r);
	memory->Delete3DArray(v);
	memory->Delete3DArray(f);
	memory->Delete3DArray(lj);
	memory->Delete3DArray(fdrag);
	memory->Delete3DArray(press);

	memory->Delete2DArray(F);
	memory->Delete2DArray(LJ);
	memory->Delete2DArray(FDRAG);
	memory->Delete2DArray(pe);
	memory->Delete2DArray(ke);
	memory->Delete2DArray(sxx);
	memory->Delete2DArray(sxy);
	memory->Delete2DArray(syy);
	memory->Delete2DArray(sxz);
	memory->Delete2DArray(syz);
	memory->Delete2DArray(szz);

	memory->Delete3DArray(pol);

	memory->Delete2DArray(ix);
	memory->Delete2DArray(iy);
	memory->Delete2DArray(iz);

	memory->Delete1DArray(id);
	memory->Delete1DArray(type);
	memory->Delete1DArray(mol);
	memory->Delete1DArray(rad);
	memory->Delete1DArray(mass);

}


void Particle::Unwrap(const std::vector<int>& group){
	for(auto i:group)
		pol[5][i][n] = pol[5][i][n] + pol[6][i][n]*M_PI;
}

double Particle::MaxRadius(){
	double *rad = atom->rad;
	double maxrad = 0.0;
	for(int i=0; i<=simbox->natoms; i++)
		if(rad[i]>maxrad)
			maxrad = rad[i];

	return maxrad;
}
