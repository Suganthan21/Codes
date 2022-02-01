#include<cmath>
#include "compute_temporal.h"
#include "memory.h"
#include "filehandle.h"
#include "simbox.h"
#include "thermo.h"
#include "particle.h"

Compute_temporal::Compute_temporal(Main* ptr):Pointers(ptr){
	CalibFile = new ifstream;
	fp1 = nullptr;
	
	string nstart = to_string(filehandle->nstart);
	string nevery = to_string(filehandle->nevery);
	string nstop  = to_string(filehandle->nstop);
	
	calibratefile = "calibrate_"+nstart+ "-"+ nevery+"-"+nstop+".dat";

	tavgrows = 120000;

	itstep = filehandle->itstep;
	ntimes = itstep.size();

	q = memory->Create2DArray<double>(tavgrows,2);
	msd = memory->Create2DArray<double>(tavgrows,2);
	msd_theta = memory->Create2DArray<double>(tavgrows,2);
	msd_phi = memory->Create2DArray<double>(tavgrows,2);
}

Compute_temporal::~Compute_temporal(){
	delete fp1;
	delete CalibFile;

	memory->Delete2DArray(q);
	memory->Delete2DArray(msd);
	memory->Delete2DArray(msd_theta);
	memory->Delete2DArray(msd_phi);

}

void Compute_temporal::EmptyArrays(){
	msd       = memory->Empty2DArray(msd, tavgrows, 2);
	msd_theta = memory->Empty2DArray(msd_theta, tavgrows, 2);
	msd_phi   = memory->Empty2DArray(msd_phi, tavgrows, 2);
	q         = memory->Empty2DArray(q, tavgrows, 2);
}

void Compute_temporal::CalibrateTdiff(){

	printf("Calibrating Time Difference\n");
	
	FILE* fp;
	fp = fopen(calibratefile.c_str(),"w");

	int dt, count=0;
	bool add;
	
	vector<int> t, pos;

	int i, ii;
	for(i=0; i < ntimes; i++){
		ii=i;
		while(ii >= 0){
			add = 1;
			dt = itstep[i]-itstep[ii];

			int len = t.size();
			//Check if the tdiff is already added
			for(int k=0; k < len; k++){
				if (dt == t[k]){
					add=0;
					pos.push_back(k);
				}
			}
			if(add){
				t.push_back(dt);
				pos.push_back(count);
				count++;
			}

			fprintf(fp,"%d\n",pos.back());

// 			printf("%d  %d  %d  %d  %d\n",itstep[i],itstep[ii],dt,pos.back(),add)

			ii--;
		}
	}
	printf("Unique tdiff = %d, Total tdiff = %lu\n",count,pos.size());
	fclose(fp);

}

void Compute_temporal::TimeOriginAvg(const vector<int>& group){
	
	if(n==0){
		std::ifstream infile(calibratefile);
		if( !infile.good() )
			CalibrateTdiff();
		CalibFile->open(calibratefile);
		EmptyArrays();
		fprintf(log,"\nComputing Time Origin Average\n");
	}
	string str1;

	double dt;
	int nn=n,l;

	while(nn >= 0){
		getline(*CalibFile, str1);
		l = stoi(str1);

		q[l][2] += 1.0;
		msd[l][2] += 1.0;
		msd_theta[l][2] += 1.0;
		msd_phi[l][2] += 1.0;

		if(q[l][2] < 1.1){  //double
			dt = itstep[n-1] - itstep[nn-1];
			tdiff.push_back(dt);
		}
		
// 		PositionMSD(l, nn, group);
		AngularMSD(l, nn, group);

		nn--;
	}

	if(n==ntimes-1){
		WriteTimeOriginAvg(group.size());
		CalibFile->close();
	}
}

void Compute_temporal::PositionMSD(int l, int nn, const vector<int>& group){
	double ***r = atom->r;
	double L = simbox->L;
	double delx, dely, delz, dr;
	for(auto i:group){

		delx = r[0][i][n] - r[0][i][nn];	delx -= round(delx/L)*L; 
		dely = r[1][i][n] - r[1][i][nn];	dely -= round(dely/L)*L;
		delz = r[2][i][n] - r[2][i][nn];	delz -= round(delz/L)*L;

		if(simbox->dim ==2) delz = 0.0;

		dr = sqrt(delx*delx + dely*dely + delz*delz);

		msd[l][1] += dr*dr;

		if(dr <= 0.3){
			q[l][1] += 1.0;	//qvalue
		}
	}
	// 		printf("%d  %d  %d  %g  %g  %g\n",n, itstep[n-1],itstep[nn-1],dr,q[l][1],l,q[l][2])
}


void Compute_temporal::AngularMSD(int l, int nn, const vector<int>& group){

	double dtheta, dphi;
	double ***pol = atom->pol;
	for(auto i:group){
		dphi = pol[0][i][nn] - pol[0][i][n];
		dtheta = pol[5][i][nn] - pol[5][i][n];
		msd_phi[l][1] += dphi*dphi;
		msd_theta[l][1] += dtheta*dtheta;
	}

}


void Compute_temporal::WriteTimeOriginAvg(int count){

	double **arr1 = msd_theta; 
	double **arr2 = msd_phi; 
	double dt = simbox->dt;

	string fname = "msd_taup"+simbox->TAUP+".dat";

	FILE* fp;
	fp = fopen(fname.c_str(),"w");

	for(unsigned int i=1; i<tdiff.size(); i++){
		if(arr1[i][2] >= 10){
			fprintf(fp,"%g %g %g %g\n",tdiff[i]*dt, arr1[i][1]/(double(count)*arr1[i][2]), 
					arr2[i][1]/(double(count)*arr2[i][2]), arr1[i][2]);
		}
	}
	fclose(fp);
}
