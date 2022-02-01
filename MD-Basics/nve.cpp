#include<iostream>
#include<stdio.h>
#include <stdlib.h>
#include<cstdlib>
#include<cmath>
#include <iomanip>
#include<fstream>
#include<vector>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>

const double L = 9, nx = 6;   // nx - no of particles on one axis. chosen such that number density = 0.3
const double T = 0.6,  total_time = 10, dt = 0.005;
const double rc = L/2, rc2 = pow(rc,2) ,  ecut = 4*(  pow( (1/rc), 12) - pow( (1/rc), 6)  ), rv = (rc+0.25), rv2 = pow(rv,2); 

using namespace std;

class particle{

	public:
		double x, y, z, xp, yp, zp, fx, fy, fz, fxp, fyp, fzp, vx, vy, vz, dispx, dispy, dispz;  //xp - previous position
		int nlist[300]; int listsize;

		particle(){
			x = y = z = xp = yp = zp = vx = vy = vz = fx = fy = fz = fxp = fyp = fzp = dispx = dispy = dispz = 0.0;
		}

};


//--------------------------------------------------------------------Initialise-simple cubic-------------------------------------------------------------------


vector<double> initialise(particle atom[]){

	double pos[10];
	double dist = (double)(L-1)/(double)(nx-1);
	cout << "dist =" << dist<<"\n";
	pos[0] = 0.5;
	for(int i = 1; i < nx; i++) {
		pos[i] = pos[i-1] + dist ; 
		//cout << pos[i] << "\n";
	}

	//simple cubic
	int n=0, q=0;
	double svx=0, svx2=0, svy=0,svy2=0, svz=0, svz2=0, sv2=0, fs;

	gsl_rng *r;
	r =  gsl_rng_alloc(gsl_rng_mt19937);

	for(int i=0; i < nx; i++){
		for(int j=0; j < nx; j++){
			for(int k=0; k < nx; k++){

				atom[n].x = pos[i];   atom[n].vx = gsl_ran_gaussian(r,1.0); 		svx += atom[n].vx;
				atom[n].y = pos[j];   atom[n].vy = gsl_ran_gaussian(r,1.0);		svy += atom[n].vy;
				atom[n].z = pos[k];   atom[n].vz = gsl_ran_gaussian(r,1.0); 		svz += atom[n].vz;
				n++; q += 3;
			}
		}
	}

	int npart = n;

	cout<< "total particles = "<< npart <<"\n";

	svx /= npart; 		svy /= npart;			svz /= npart;

	double vcm1 = sqrt(svx*svx + svy*svy + svz*svz);

	cout << svx << "\t" << svy << "\t" << svz << "\t" << "\n vcm =" << vcm1 << "\n" ;

	// vcm = 0 and temperature is maintained, previous positions recorded.


	svx2 = svy2 = svz2 = sv2 =0;

	for(int i=0; i < npart; i++){
		atom[i].vx = (atom[i].vx - svx);			svx2 += atom[i].vx*atom[i].vx;
		atom[i].vy = (atom[i].vy - svy);			svy2 += atom[i].vy*atom[i].vy;
		atom[i].vz = (atom[i].vz - svz);			svz2 += atom[i].vz*atom[i].vz;

	}


	sv2 = (svx2 + svy2 + svz2);
	fs = sqrt(3*npart*T/sv2); 
	cout << "fs = " << fs << "\n";

	double d1,d2,d3; 
	svx2 = svy2 = svz2 = sv2 = 0;

	for(int i=0; i<npart; i++){
		atom[i].vx *= fs; 		svx2 += atom[i].vx*atom[i].vx;
		atom[i].vy *= fs; 		svy2 += atom[i].vy*atom[i].vy;
		atom[i].vz *= fs; 		svz2 += atom[i].vz*atom[i].vz;

		d1 += atom[i].vx; 		d2 += atom[i].vy;		d3 += atom[i].vz;
	}

	vcm1 = sqrt(d1*d1 + d2*d2 + d3*d3);
	cout << "New vcm =" << vcm1 << "\n";

	sv2 = (svx2 + svy2+ svz2);
	cout << "initial temp = \t" << (sv2/(3*npart)) << "\n";

	vector<double> initial;
	initial.push_back(npart);
	initial.push_back(sv2);
	initial.push_back(vcm1);

	return initial;
}

//--------------------------------------------------------------------Neighbour list----------------------------------------------------------------------------

void verletlist(particle* atom, int npart){
	
	for(int i=0;i<npart; i++){
		atom[i].listsize =0;
	}
	int k1, k2;
	for(int i=0; i<npart-1; i++){
		for(int j=i+1; j<npart; j++){
			double dx, dy, dz, r2;
		
			dx = atom[i].x - atom[j].x; dx -= round(dx/L)*L; 
			dy = atom[i].y - atom[j].y;	dy -= round(dy/L)*L;
			dz = atom[i].z - atom[j].z;	dz -= round(dz/L)*L;	

			r2 = dx*dx + dy*dy + dz*dz;	
			
			if(r2 < rv2){
				k1 = atom[i].listsize; 	
				atom[i].nlist[k1] = j;	atom[i].listsize++;
						
			}			

		}
	}
}

bool checkdisp(particle* atom, int npart){
	
	double dx, dy, dz,r2, maxdr=0.0;
	bool check=false;
	
	for(int i=0; i<npart; i++){	

		dx = atom[i].x - atom[i].xp;	dx -= round(dx/L)*L; 
		dy = atom[i].y - atom[i].yp;	dy -= round(dy/L)*L;
		dz = atom[i].z - atom[i].zp;	dz -= round(dz/L)*L;

		atom[i].dispx += dx;
		atom[i].dispy += dy;
		atom[i].dispz += dz;

		r2 = atom[i].dispx*atom[i].dispx + atom[i].dispy*atom[i].dispy + atom[i].dispz*atom[i].dispz;	

		if(r2 > maxdr){
			maxdr = r2;
		}

	}

	if(sqrt(maxdr) > (rv-rc)/2){
		check = true;
		for(int i=0; i<npart; i++){
			atom[i].dispx = atom[i].dispy = atom[i].dispz =0;
		}
	}

	return check;

}

//--------------------------------------------------------------------Leonard Jones Force-------------------------------------------------------------------


double force(particle* atom, int npart){
	double dx, dy, dz, ff, r2, r2i, r6i, en=0;
	int neighbour;
	//force on every particle

	for(int i=0; i < npart; i++){
		atom[i].fx = 0; atom[i].fy = 0; atom[i].fz = 0;
	}

	for(int i=0; i<npart; i++){

		for(int j=0; j<(atom[i].listsize); j++){

			neighbour = atom[i].nlist[j];

			dx = atom[i].x - atom[neighbour].x; dx -= round(dx/L)*L; 
			dy = atom[i].y - atom[neighbour].y;	dy -= round(dy/L)*L;
			dz = atom[i].z - atom[neighbour].z;	dz -= round(dz/L)*L;

			r2 = dx*dx + dy*dy + dz*dz;
			//cout << dx << "\n";
			if(r2 < rc2){
				r2i = 1/r2;
				r6i = r2i*r2i*r2i;
				ff = 48*r2i*r6i*(r6i-0.5);

				atom[i].fx += ff*(dx); 	atom[neighbour].fx -= ff*(dx);
				atom[i].fy += ff*(dy); 	atom[neighbour].fy -= ff*(dy);
				atom[i].fz += ff*(dz); 	atom[neighbour].fz -= ff*(dz);
				en = en + (4*r6i*(r6i-1)-ecut);
			}
		}
		//cout<<"force\t"<< i <<"\t" << en <<"\t"<< atom[i].fx<<"\t"<< atom[i].fy << "\t" << atom[i].fz<<  "\n";
		// cout << atom[i].x << "\t" << atom[i].y << "\t" << atom[i].z << "\n";
	}
	return en;
}


//-----------------------------------------------------------------------verlet-----------------------------------------------------------------------------------

vector<double> verlet(particle* atom, int npart){

	double svx=0.0, svy=0.0, svz=0.0, svx2, svy2, svz2, sv2, vcm;
	vector<double> output;

	for(int i=0; i<npart; i++){

		//previous positions
		atom[i].xp = atom[i].x;
		atom[i].yp = atom[i].y;
		atom[i].zp = atom[i].z;

		//velocity-verlet algorithm
		atom[i].x += atom[i].vx*dt + 0.5*atom[i].fx*dt*dt;
		atom[i].y += atom[i].vy*dt + 0.5*atom[i].fy*dt*dt;
		atom[i].z += atom[i].vz*dt + 0.5*atom[i].fz*dt*dt;

		//boundary conditions
		if(atom[i].x > L) atom[i].x -= L; 	if(atom[i].x < 0) atom[i].x += L; 
		if(atom[i].y > L) atom[i].y -= L; 	if(atom[i].y < 0) atom[i].y += L; 
		if(atom[i].z > L) atom[i].z -= L; 	if(atom[i].z < 0) atom[i].z += L;	


	}

	for(int i=0; i<npart; i++) {
		atom[i].fxp = atom[i].fx;
		atom[i].fyp = atom[i].fy;
		atom[i].fzp = atom[i].fz;
	}

	double pot = force(atom, npart);

	for(int i=0; i<npart; i++){

		atom[i].vx += 0.5*dt*(atom[i].fx + atom[i].fxp);
		atom[i].vy += 0.5*dt*(atom[i].fy + atom[i].fyp);
		atom[i].vz += 0.5*dt*(atom[i].fz + atom[i].fzp);

		//cout << atom[i].vx << "\t" <<atom[i].vy << "\t" <<atom[i].vz << "\n";
		svx += atom[i].vx; svy += atom[i].vy; svz += atom[i].vz;

		svx2 +=  (atom[i].vx)*(atom[i].vx); 	svy2 +=  (atom[i].vy)*(atom[i].vy); 	svz2 +=  (atom[i].vz)*(atom[i].vz);

	}


	//cout << svx2 << "\t" << svy2 << "\t" <<svz2 << "\n";

	svx /= npart;       svy /= npart;         svz /= npart;
	sv2 = (svx2 + svy2 + svz2);

	vcm = sqrt(svx*svx + svy*svy + svz*svz);

	//cout << "2*ke =\t"<< sv2 << "\tvcm = \t" << vcm << "\n";

	output.push_back(pot);
	output.push_back(sv2);
	output.push_back(vcm);

	return output;

}



int main(){	

	// 	int files = 3;
	// 	string filename[files]={"output.dat","vel_initial.dat","velocity.dat"};
	// 
	// 	ofstream outputFile[files];
	// 	for (int i = 0; i < files; i++){
	// 		outputFile[i].open(filename[i].c_str());
	// 	}	
	// 
	// 	string positionfile[1000];
	// 	for (int i = 0; i < 1000; i++){
	// 		positionfile[i]="time"+std::to_string(i)+".xyz";
	// 	}
	// 
	// 	ofstream outputfile[1000];
	// 	for (int i = 0; i < 1000; i++){
	// 		outputfile[i].open(positionfile[i].c_str());
	// 		outputfile[i]<<"216"<<"\n"<<"#Simple cubic - velocity verlet"<<"\n";
	// 	}
	// 


	int size = round(total_time/dt); 
	cout << "size =\t" << size << "\n";

	particle* atom = new particle[600];


	//-------------------------------------------------------------------calculate quantities-------------------------------------------------------------------

	double* pe = new double[size];  double* ke = new double[size];   double* etot = new double[size];	double* temp = new double[size];	
	double* vcm = new double[size];	double* current_time = new double[size]; 
	
	bool check;
	vector<double> initial, output; int npart;

	initial = initialise(atom);

	npart = static_cast<int>(initial[0]);	
	ke[0] = 0.5*initial[1];
	vcm[0] = initial[2];
	temp[0] = 2*ke[0]/(3*npart);

	//Find nearest neighbours
	verletlist(atom, npart);	

	pe[0] = force(atom, npart); 
	etot[0] = (pe[0]+ke[0])/npart;

	cout << ke[0] <<"\t" << vcm[0] << "\n";


	current_time[0]=0; 	
	int time=0;	
	clock_t clk_start, clk_end;
	clk_start = clock();

	for(int t=0; t < size; t++){

		cout << setprecision(5) << fixed << current_time[t] << "\t" << pe[t] << "\t" << ke[t] << "\t" << etot[t] << "\t" << vcm[t] << "\t" << temp[t] <<"\n";

		current_time[t+1] = current_time[t] + dt;

		// 		if(t<500 || t>(size-500)){
		// 
		// 			for(int i=0; i< npart; i++){
		// 				outputfile[time] << setprecision(5) << fixed << atom[i].x <<"\t" << atom[i].y <<"\t"<< atom[i].z << "\n";
		// 			}
		// 			time++;
		// 		}
		// 
		// 		if(t==0){  //initial velocities are stored
		// 			for(int i=0; i < npart; i++)
		// 				outputFile[1] << setprecision(5) << fixed << atom[i].vx <<"\t" << atom[i].vy <<"\t"<< atom[i].vz << "\n";
		// 		} 
		// 
		// 
		// 		if(t > (size-1000)){		//final 1000 timestep velocities are stored	
		// 			for(int i=0; i < npart; i++)				
		// 				outputFile[2] << setprecision(5) << fixed << atom[i].vx <<"\t" << atom[i].vy <<"\t"<< atom[i].vz << "\n";	
		// 		}

	
		output = verlet(atom, npart);	
		pe[t+1] = output[0];						//potential energy
		ke[t+1] = 0.5*output[1];					//kinetic energy	
		vcm[t+1] = output[2];						//velocity centre of mass
		etot[t+1] = (pe[t+1] + ke[t+1])/npart;		//total energy

		temp[t+1] = 2*ke[t+1]/(3*npart);			//instantaneous temperature

		check = checkdisp(atom, npart);

		if(check){
// 			cout << t << "\n";
			verletlist(atom, npart);
		}


//      outputFile[0]<<setprecision(5)<<fixed<< current_time[t] << "\t" << pe[t] << "\t" << ke[t] << "\t" << etot[t] << "\t" << vcm[t] << "\t" << temp[t] <<"\n";

	}
	clk_end = clock();

	cout << setprecision(5) << "Time taken = " << double(clk_end - clk_start)/(CLOCKS_PER_SEC) << "\n";


	// 	for (int i = 0; i < files; i++){  
	// 		outputFile[i].close();
	// 	}	
	// 
	// 	for (int i = 0;  i<time; i++){                   // i is time index and j is particle
	// 		outputfile[i].close();
	// 	}

	return 0;
}



