#include<cmath>
#include<map>
#include "compute_spatial.h"
#include "simbox.h"
#include "particle.h"
#include "neigh.h"
#include "filehandle.h"

Compute_spatial::Compute_spatial(Main* ptr, int ntimes):Pointers(ptr){
	fp1 = NULL;
	map_coarse_corr.resize(ntimes);
	flag_spacecorr = 0;
	flag_filledsquares = 0;
}

Compute_spatial::~Compute_spatial(){
	delete fp1;
}

void Compute_spatial::InitSquares(const vector<int>& group){

	flag_filledsquares = 1;
	nsquare = group.size();
	double spacing = (simbox->L-1.0)/round(sqrt(nsquare));
	sq.clear();
	sq.resize(nsquare+1);

	int j;
	double diff=spacing/2.0;
	double lty=-spacing;

	for(int i=0; i<nsquare; i++){
		if(i<10) j=i;
		else j=i%10;
		if(j==0) lty += spacing;
		sq[i].initxhi = diff + j*spacing;
		sq[i].xhi = sq[i].initxhi;
		sq[i].yhi = diff+lty;

		sq[i].initxlo = sq[i].initxhi - spacing;
		sq[i].xlo = sq[i].initxlo;
		sq[i].ylo = sq[i].yhi - spacing;

		sq[i].initxmid = 0.5*(sq[i].initxlo + sq[i].initxhi) ;
		sq[i].ymid = 0.5*(sq[i].ylo + sq[i].yhi);
		sq[i].xmid = 0.5*(sq[i].xlo + sq[i].xhi);
	}
}

void Compute_spatial::FilledSquares(const vector<int>& group){
	if(n==0) InitSquares(group);
	
	double boxtheta = simbox->theta;
	double L = simbox->L;
	double*** r = atom->r;

	double excess;
	double tiltsign = simbox->tiltsign;

	// boundary of squares is adjusted as box is tilted
	if(simbox->shearflag){
		for(int i=0; i<nsquare; i++){
			excess = abs(tan(boxtheta)*sq[i].yhi);
			sq[i].xlo = sq[i].initxlo + tiltsign*excess;
			sq[i].xhi = sq[i].initxhi + tiltsign*excess;

			excess = abs(tan(boxtheta)*sq[i].ymid);
			sq[i].xmid = sq[i].initxmid + tiltsign*excess;

// 			OutfileDumpLattice << i << "\t" << sq[i].xhi << "\t" << sq[i].yhi << "\t0.0\n";
		}
	}

	int nfill = 0;
	double delx, dely;

	for(int j=0; j<nsquare; j++)
		sq[j].fill =0;

	double limitlo, limithi, sqxmid;
	int j;
	for(auto i:group){
		for(j=0; j<nsquare; j++){
			dely = r[1][i][n] - sq[j].ymid;	dely -= round(dely/L)*L;
			dely = abs(dely);

			if(sq[j].ymid+dely > sq[j].yhi || sq[j].ymid-dely < sq[j].ylo )
				continue;

			excess = tiltsign*abs(tan(boxtheta))*( sq[j].yhi - r[1][i][n]);
			limitlo = sq[j].xlo - excess;
			limithi = sq[j].xhi - excess;
			sqxmid = 0.5*(limitlo + limithi);

			delx = r[0][i][n] - sqxmid;   delx -= round(delx/L)*L; 
			delx = abs(delx);
			
			if( sqxmid + delx  < limithi &&  sqxmid - delx > limitlo)
				sq[j].fill += 1;
		}
	}

	for(int j=0; j<nsquare; j++){
		if(sq[j].fill >= 1)
			nfill++;
	}
	fprintf(fp1,"%g\n",(double)nfill/(double)nsquare);
}

void Compute_spatial::Correlation(const vector<int>& group){
	
	printf("Computing Spatial Correlation...\n");
	flag_spacecorr = 1;

	double ***q = atom->v; 
	double ***r = atom->r;
	double r2, delx, dely,delz;
	double L = simbox->L;

	double Q1, Q2, dist, corr;
	multimap <double,double> map_corr;

	adder.clear();
	counter.clear();
	
	//compute spatial correlation till cutoff of nlist
	neigh->Build(simbox->L); 
	int j,jj;
	for(auto i:group){
		for(jj=0; jj < neigh->nlistsize[i]; jj++){
			j = neigh->nlist[i][jj];
			delx = r[0][i][n] - r[0][j][n]; delx -= round(delx/L)*L; 
			dely = r[1][i][n] - r[1][j][n];	dely -= round(dely/L)*L;
			delz = r[2][i][n] - r[2][j][n];	delz -= round(delz/L)*L;

			r2 = delx*delx + dely*dely + delz*delz;
			dist = sqrt(r2);

			Q1 = sqrt(q[0][i][n]*q[0][i][n] + q[1][i][n]*q[1][i][n] + q[2][i][n]*q[2][i][n]);
			Q2 = sqrt(q[0][j][n]*q[0][j][n] + q[1][j][n]*q[1][j][n] + q[2][j][n]*q[2][j][n]);

			corr = Q1*Q2;

			map_corr.insert({dist,corr});
		}
		CoarsenCorrelation(map_corr);
		map_corr.clear();
	}
	
	double avg=0.0;
	for(auto i:group)
		avg += sqrt(q[0][i][n]*q[0][i][n] + q[1][i][n]*q[1][i][n] + q[2][i][n]*q[2][i][n]);
	avg /= group.size();
	
	for(unsigned int i=1; i<=counter.size(); i++){ //should not start from zero
		dist = adder[i];
		corr = adder[-i];
		dist /= counter[i];
		corr /= (counter[i]*avg*avg);
		map_coarse_corr[n].insert({dist,corr}); 
	}
	
	if(n==filehandle->ntimes-1)
		FinalWrite();
	
}

void Compute_spatial::CoarsenCorrelation(const multimap<double,double>& map_corr){

	if(map_corr.size() <= 0)
		return;

	double max_rad, binwidth = 1.0; //binwidth = LJ dia
	multimap <double,double>::const_iterator it = map_corr.begin();
	max_rad = binwidth;
	
	int bin = 0;
	while(max_rad < map_corr.rbegin()->first){
		bin++;
		// adder : positive indices - distance, negative indices - correlation
		while(it->first < max_rad){
			adder[bin] += it->first;
			adder[-bin] += it->second; 
			counter[bin] += 1;
			it++;
		}
		max_rad += binwidth;
	}

}


void Compute_spatial::FinalWrite(){
	FILE* fp;
	int ntimes = filehandle->ntimes;

	if(flag_spacecorr){
		fp = fopen("spacecorr1.dat","w");
		multimap <double,double>::iterator it[ntimes];
		
		fprintf(fp,"# Timesteps = ");
		for(int i=0;i<ntimes;i++){
			it[i] = map_coarse_corr[i].begin();
			fprintf(fp,"%d\t",filehandle->itstep[i]);
		}
		fprintf(fp,"\n");
		
		int len = map_coarse_corr[n].size();
		for(int i=0; i<len; i++){
			for(int j=0; j<ntimes; j++){
				fprintf(fp,"%g %g\t",it[j]->first,it[j]->second);	
				it[j]++;
			}
			fprintf(fp,"\n");
		}
		fclose(fp);
	}

}
