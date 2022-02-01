#include<cmath>
#include "pair.h"
#include "memory.h"
#include "filehandle.h"
#include "particle.h"
#include "simbox.h"
#include "neigh.h"

using std::vector;

Pair::Pair(Main* ptr):Pointers(ptr){
	pair_time = 0.0;
	gamma_dpd = 1.0;
	cutoff_dpd = 2.5;
	gamma_ld = 1.0;

	epsilon[1][1] = 1.0;
	epsilon[1][2] = 1.5;
	epsilon[2][1] = 1.5;
	epsilon[2][2] = 0.5;

	sigma[1][1] = 1.0;
	sigma[1][2] = 0.8;
	sigma[2][1] = 0.8;
	sigma[2][2] = 0.88;

	rc2[1][1] = 2.5*2.5;
	rc2[1][2] = 2.0*2.0;
	rc2[2][1] = 2.0*2.0;
	rc2[2][2] = 2.2*2.2;

	cutoffsq_dpd = cutoff_dpd*cutoff_dpd;

	for(int i=1; i<=2; i++){
		for(int j=1; j<=2; j++){
			lj1[i][j] = 48.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
			lj2[i][j] = 24.0 * epsilon[i][j] * pow(sigma[i][j],6.0);
			lj3[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
			lj4[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],6.0);
			// 		  printf("%d %d %g %g",i,j,lj1[i][j],lj2[i][j]);

			double ratio = sigma[i][j] / sqrt(rc2[i][j]);
			offset[i][j] = 4.0 * epsilon[i][j] * (pow(ratio,12.0) - pow(ratio,6.0));
		}
	}
	
}

Pair::~Pair(){}

// list all the pair potentials that must be computed
void Pair::Compute(const vector<int>& group){
	clock_t clk_start = clock();

	int natoms = simbox->natoms;
	int ntimes = filehandle->ntimes;

	atom->f = memory->Empty3DArray(atom->f,3,natoms,ntimes);
	atom->press = memory->Empty3DArray(atom->press,6,natoms,ntimes);
	atom->pe = memory->Empty2DArray(atom->pe,natoms,ntimes);

	LJPolydisperse(group);

	pair_time += double(clock() - clk_start)/(CLOCKS_PER_SEC);
}

void Pair::LJ(const vector<int>& group){
	double delx, dely, delz, forcelj, r2, r2i, r6i, energy;
	int itype, jtype;

	// define array f as lj to output lj force alone
	double ***f = atom->f;
	double ***r = atom->r;
	double L = simbox->L;

	int j,jj;
	for(auto i:group){
		itype = atom->type[i];

		for(jj=0; jj < neigh->nlistsize[i]; jj++){
			j = neigh->nlist[i][jj];
			jtype = atom->type[j];

			delx = r[0][i][n] - r[0][j][n]; delx -= round(delx/L)*L; 
			dely = r[1][i][n] - r[1][j][n];	dely -= round(dely/L)*L;
			delz = r[2][i][n] - r[2][j][n];	delz -= round(delz/L)*L;

			if(simbox->dim == 2) delz =0.0;

			r2 = delx*delx + dely*dely + delz*delz;

			if(r2 < rc2[itype][jtype]){
				r2i = 1/r2;
				r6i = r2i*r2i*r2i;
				forcelj = r6i * r2i *(lj1[itype][jtype]*r6i - lj2[itype][jtype]);

				f[0][i][n] += forcelj*delx;	  f[0][j][n] -= forcelj*delx;
				f[1][i][n] += forcelj*dely;	  f[1][j][n] -= forcelj*dely;
				f[2][i][n] += forcelj*delz;	  f[2][j][n] -= forcelj*delz;

				energy = r6i*(lj3[itype][jtype]*r6i - lj4[itype][jtype]) - offset[itype][jtype];

				PerAtomProps(i,j,delx,dely,delz,forcelj,energy);
			}
		}
	}
// 	for (int i=0;i < simbox->natoms; i++)
// 		atom->LJ[i][n] = sqrt(lj[0][i][n]*lj[0][i][n]  +  lj[1][i][n]*lj[1][i][n]  + lj[2][i][n]*lj[2][i][n]);
	
}

//assumes epsilon_ij = 1
void Pair::LJPolydisperse(const vector<int>& group){
	double delx, dely, delz, r2, r2i, r6i, forcelj;

	double ***f = atom->f;
	double ***r = atom->r;
	double *rad = atom->rad;
	double L = simbox->L;
	double cutsq, sij, sij2, sij6, energy;

	int j,jj;
	for(auto i:group){
		for(jj=0; jj < neigh->nlistsize[i]; jj++){
			j = neigh->nlist[i][jj];

			delx = r[0][i][n] - r[0][j][n]; delx -= round(delx/L)*L; 
			dely = r[1][i][n] - r[1][j][n];	dely -= round(dely/L)*L;
			delz = r[2][i][n] - r[2][j][n];	delz -= round(delz/L)*L;

			if(simbox->dim == 2) delz = 0.0;

			r2 = delx*delx + dely*dely + delz*delz;
			sij = rad[i] + rad[j];
			sij2 = sij*sij;
			cutsq = sij2*1.25992104989487;

			if (r2 < cutsq ) {

				r2i = 1.0/r2;
				r6i = r2i*r2i*r2i;

				sij6 = sij2*sij2*sij2;

				forcelj = r6i*r2i*sij6*(48.0*r6i*sij6 - 24.0);

				f[0][i][n] += forcelj*delx;	  f[0][j][n] -= forcelj*delx;
				f[1][i][n] += forcelj*dely;	  f[1][j][n] -= forcelj*dely;
				f[2][i][n] += forcelj*delz;	  f[2][j][n] -= forcelj*delz;

				// epsilon = 1.0
				energy = (4.0*sij6*r6i*(sij6*r6i - 1.0) + 1.0);

				PerAtomProps(i,j,delx,dely,delz,forcelj,energy);
			}
		}
	}
}

void Pair::PerAtomProps(int i, int j, double delx, double dely, 
								double delz, double fpair, double energy)
{
	double ***press = atom->press;	
	vector<int> vec = {i,j};
	for(auto id:vec){
		press[0][id][n] += -0.5*delx*delx*fpair;
		press[1][id][n] += -0.5*dely*dely*fpair;
		press[2][id][n] += -0.5*delz*delz*fpair;
		press[3][id][n] += -0.5*delx*dely*fpair;
		press[4][id][n] += -0.5*delx*delz*fpair;
		press[5][id][n] += -0.5*dely*delz*fpair;
	}

	double **pe = atom->pe;
	pe[i][n] += 0.5*energy;
	pe[j][n] += 0.5*energy;
}

//assumes DPD cutoff = lj cutoff. If not, request neighbour list for the desired cutoff
void Pair::DragDPD(const vector<int>& group){
	double delx, dely, delz, delvx, delvy, delvz;
	double r2, dot, r2inv;

	double ***f = atom->f;
	double ***r = atom->r;
	double ***v = atom->v;
	double L = simbox->L;

	for(int i=0; i < 3; i++)
		for(int j=0; j < simbox->natoms; j++)
			f[i][j][n] = 0.0;
	int j,jj;
	for(auto i:group){
		for(jj=0; jj < neigh->nlistsize[i]; jj++){
			j = neigh->nlist[i][jj];

			delx = r[0][i][n] - r[0][j][n]; delx -= round(delx/L)*L; 
			dely = r[1][i][n] - r[1][j][n];	dely -= round(dely/L)*L;
			delz = r[2][i][n] - r[2][j][n];	delz -= round(delz/L)*L;

			if(simbox->dim == 2) delz =0.0;

			r2 = delx*delx + dely*dely + delz*delz;

			if(r2 < cutoffsq_dpd){
				r2inv = 1.0/r2;
				delvx = v[0][i][n] - v[0][j][n];
				delvy = v[1][i][n] - v[1][j][n];
				delvz = v[2][i][n] - v[2][j][n];

				dot = delx*delvx + dely*delvy + delz*delvz;
				f[0][i][n] += -gamma_dpd*dot*r2inv*delx;		f[0][j][n] -= -gamma_dpd*dot*r2inv*delx;
				f[1][i][n] += -gamma_dpd*dot*r2inv*dely;		f[1][j][n] -= -gamma_dpd*dot*r2inv*dely;
				f[2][i][n] += -gamma_dpd*dot*r2inv*delz;		f[2][j][n] -= -gamma_dpd*dot*r2inv*delz;

// 				cout << r2 << "\t" << delx << "\t" << delvx << "\t" << dely << "\t" << delvy << "\t" << dot << "\n";
			}
		}
	}
	//total drag
	for (int i=0;i <= simbox->natoms; i++)
		atom->FDRAG[i][n] = sqrt(f[0][i][n]*f[0][i][n]  +  f[1][i][n]*f[1][i][n]  + f[2][i][n]*f[2][i][n]);
	
}


void Pair::DragLD(const vector<int>& group){
	double ***f = atom->f;
	double ***v = atom->v;
	for(auto i:group){
		f[0][i][n] = -gamma_ld*v[0][i][n];
		f[1][i][n] = -gamma_ld*v[1][i][n];
		f[2][i][n] = -gamma_ld*v[2][i][n];

		if(simbox->dim == 2) v[2][i][n] =0.0;
		
	}
	//total drag
	for(auto i:group)
		atom->FDRAG[i][n] = sqrt(f[0][i][n]*f[0][i][n]  +  f[1][i][n]*f[1][i][n]  + f[2][i][n]*f[2][i][n]);
}

void Pair::TotalForce(const vector<int>& group){
	for(auto i:group)
		atom->F[i][n] = atom->LJ[i][n] + atom->LJ[i][n] + simbox->f0;
}

// 	 	double** pol1 = atom->pol1;
// 		atom[i].facx = f0 * cos(pol1[i]);
// 		atom[i].facy = f0 * sin(pol1[i]);
// 		atom[i].fac = sqrt(atom[i].facx*atom[i].facx  +  atom[i].facy*atom[i].facy);
// 
// 		atom[i].f = sqrt(atom[i].fx*atom[i].fx  +  atom[i].fy*atom[i].fy );
// 
// 		atom[i].FX = atom[i].ljx + atom[i].fx + atom[i].facx;
// 		atom[i].FY = atom[i].ljy + atom[i].fy + atom[i].facy;
