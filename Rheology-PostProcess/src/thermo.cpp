#include<cmath>
#include "thermo.h"
#include "memory.h"
#include "filehandle.h"
#include "particle.h"
#include "simbox.h"
#include "pair.h"

using std::vector;

Thermo::Thermo(Main* ptr, int ntimes):Pointers(ptr){
	step = memory->Create1DArray<int>(ntimes);
	pe = memory->Create1DArray<double>(ntimes);
	sxy = memory->Create1DArray<double>(ntimes);
	fmax = memory->Create1DArray<double>(ntimes);
	fnorm = memory->Create1DArray<double>(ntimes);
	pol = memory->Create1DArray<double>(ntimes);

	ke = memory->Create2DArray<double>(7,ntimes);
	temp = memory->Create2DArray<double>(7,ntimes);
	press = memory->Create2DArray<double>(7,ntimes);
	custom = memory->Create2DArray<double>(40,ntimes);
}

Thermo::~Thermo(){
	memory->Delete1DArray(step);
	memory->Delete1DArray(pe);
	memory->Delete1DArray(sxy);
	memory->Delete1DArray(fmax);
	memory->Delete1DArray(fnorm);
	memory->Delete1DArray(pol);

	memory->Delete2DArray(ke);
	memory->Delete2DArray(temp);
	memory->Delete2DArray(press);
	memory->Delete2DArray(custom);
}

// void Thermo::Reset(){
// 	int ntimes = filehandle->ntimes;
// 	ke = memory->Empty2DArray(thermo->ke,7,ntimes); 
// 	temp = memory->Empty2DArray(thermo->temp,7,ntimes); 
// 	press = memory->Empty2DArray(thermo->press,7,ntimes);
// 	pe = memory->Empty1DArray(thermo->pe,ntimes);
// 
// 	peflag = 0;
// 	keflag = 0;
// 	pressflag = 0;
// }

void Thermo::PE(const vector<int>& group){
	pe[n] = 0.0;
	for(auto i:group)
		pe[n] += atom->pe[i][n];
	pe[n] /= double(simbox->natoms);
}

void Thermo::KE(const vector<int>& group){
	for(int i=0;i<7;i++){
		ke[i][n] = 0.0;
		temp[i][n] = 0.0;
	}
	double*** v = atom->v;

	for(auto i:group){
		ke[0][n] += 0.5*(v[0][i][n]*v[0][i][n] + v[1][i][n]*v[1][i][n] + v[2][i][n]*v[2][i][n]); //KE
		ke[1][n] += 0.5*(v[0][i][n]*v[0][i][n]); //KE - xx 
		ke[2][n] += 0.5*(v[1][i][n]*v[1][i][n]); //KE - yy 
		ke[3][n] += 0.5*(v[2][i][n]*v[2][i][n]); //KE - zz 
		ke[4][n] += 0.5*(v[0][i][n]*v[1][i][n]); //KE - xy 
		ke[5][n] += 0.5*(v[0][i][n]*v[2][i][n]); //KE - xz 
		ke[6][n] += 0.5*(v[1][i][n]*v[2][i][n]); //KE - yz 

	}
	for(int i=0;i<7;i++){
		temp[i][n] = 2.0*ke[i][n]/double(simbox->dim*simbox->natoms);
		ke[i][n] /= double(simbox->natoms);
	}
}

void Thermo::ETOT(const vector<int>& group){
	etot[n] = pe[n] + ke[0][n]; 
}

void Thermo::PRESS(const vector<int>& group){
	KE(group);

	for(int i=0;i<7;i++)
		press[i][n] = 0.0;
	

	for(int i=1;i<7;i++)
		for(auto j:group)
			press[i][n] += atom->press[i-1][j][n];
	
	double vol = pow(simbox->L,3.0);
	for(int i=1;i<7;i++)
		press[i][n] *= 1.0/vol;

	temp[0][n] = 5.0;
	double ideal = double(simbox->natoms)*temp[0][n]/vol;
	press[0][n] = (-1.0/double(simbox->dim))*(press[1][n] + press[2][n] + press[3][n]);
	press[0][n] += ideal;

}
