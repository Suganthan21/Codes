#include<iomanip>
#include<cmath>
#include<sstream>
#include "filehandle.h"
#include "simbox.h"
#include "group.h"
#include "thermo.h"
#include "particle.h"

using namespace std;

FileHandle::FileHandle(Main* ptr):Pointers(ptr){
	readdump = 0;
	writedump = 0;
	Idump = nullptr;
	Odump = nullptr;
	Idata = nullptr;
}

FileHandle::~FileHandle(){
	delete Idump;
	delete Idata;
}

void FileHandle::Init(){
// 	log = fopen("log.process","w");
	readdump = 1;
	nstart = 0;
	nstop = 5000000;
	nevery = 1000;
	
	simbox->srate = 0.0003;
	simbox->SRATE = "0.0003";

	FileName();

//  strain   1.5 2.0 2.5 3.0 3.5 4.0 4.5
//  sr0.01  150000,200000,250000,300000,350000,400000,450000
//  sr0.003  500010,666680,833350,1000020,1166690,1333360,1500030
//  sr0.001   1500015,2000020,2500025,3000030,3500035,4000040,4500045
//  sr0.0003  5000010,6666680,8333350,10000020,11666690,13333360,15000030
	itstep.insert(itstep.end(),{5000010,6666680,8333350,10000020,11666690,13333360,15000030});

	DumpSteps("custom"); //file or equispaced or custom(itstep should be inserted above)

	if(idump.size()<1 && idata.size()<1){
		fprintf(log,"ERROR: NO INPUT FILES\n");
		throw exception();
	}
}

void FileHandle::FileName(){

	string path = "/home/suganth/lammps-3Mar20/residual/config/";
	
	idump.push_back(path+"VF0.68068-N97556-FCC-seed111-Rate0.0005-Minimize_SRate"+simbox->SRATE+"_PBC-1.dump");
	
// 	odump.push_back("config1.dump");

// 	for(int i=0;i<4; i++){
// 		idump.push_back(path+"config-3d_sp_taup"+tp[i]+"_f1.dump");
// 		idata.push_back("");
// 		odump.push_back("");
// 	}
}

void FileHandle::DumpSteps(string type){
	fprintf(log,"Storing timesteps to read...");

	if(type == "file"){
		FileGood("steps.dat");
		ifstream dstepfile;
		dstepfile.open("steps.dat");
		string str1;
		const char *ch1; 
		while(!dstepfile.eof()){
			getline(dstepfile, str1);
			ch1 = str1.c_str();
			if(ch1[0] != '#')
				itstep.push_back(atoi(ch1));
		}
		dstepfile.close();
	}
	else if(type == "equispaced"){
		if( (nstop - nstart)%nevery != 0 ){
			fprintf(log,"ERROR: CHECK VARIABLES FOR TIMESTEPS TO READ\n");
			throw exception();
		}
		int currstep;
		itstep.push_back(nstart);
		currstep= nstart;
		while (currstep <= nstop-nevery){
			currstep += nevery;
			itstep.push_back(currstep);
		}
	}
	else{
		if(type != "custom"){
			fprintf(log,"ERROR: STEPS TO READ CANNOT BE STORED\n");
			throw exception();
		}
	}
	for(auto i:itstep)
		ststep.push_back(std::to_string(i));

	ntimes=ststep.size();
	fprintf(log,"Done\n");
}


//Check if the ifstream file is present in folder
void FileHandle::FileGood(string filename){
	ifstream checkfile(filename);
	if(!checkfile.good()){
		fprintf(log,"ERROR: FILE %s NOT FOUND\n", filename.c_str());
		throw exception();
	}
	checkfile.close();
}

void FileHandle::OpenClose(string action, int i){


	if(filetype != "dump" && filetype != "data" && filetype != "dump,data"){
		fprintf(log,"ERROR: CHECK FILETYPE\n");
		throw exception();
	}

	if(action=="open"){
		if(filetype == "dump" || filetype == "dump,data"){
			if(readdump){
				Idump = new ifstream;
				Idump->open(idump[i]);
				fprintf(log,"\nReading from\n\t%s\n",idump[i].c_str());
				fprintf(log,"\tTotal number of steps to read = %d\n",ntimes) ;
			}
			if(writedump){
				fprintf(log,"Writing to \n");
				Odump=fopen(odump[i].c_str(),"w");
				fprintf(log,"\t%s\n", odump[i].c_str());
			}
		}
		else if(filetype == "data" || filetype == "dump,data"){
			Idata = new ifstream;
			Idata->open(idata[i]);
			fprintf(log,"\nReading from\n\t%s\n", idata[i].c_str());
			FileGood(idata[i]);
		}
	}
	else if(action=="close"){
		if(filetype == "dump" || filetype == "dump,data"){
			if(readdump){
				Idump->close();
				delete Idump;
			}
			if(writedump)
				fclose(Odump);
		}
		else if(filetype == "data" || filetype == "dump,data"){
			Idata->close();
			delete Idata;
		}
	}
	else{
		fprintf(log,"ERROR: CHECK ARGS TO OPEN/CLOSE FILES\n");
		throw exception();
	}
}


//Take the read pointer(seekg) to the first per-atom info line of the desired timestep.
void FileHandle::ParseDump(string gotostep){
	string str1;
	bool tcheck=0;
	string strtimestep="ITEM: TIMESTEP";

	while(!tcheck && ! Idump->eof() ){
		getline(*Idump,str1);
		if(str1 == strtimestep){
			getline(*Idump,str1);
			if( str1 == gotostep ){
				tcheck=1;
				for(int i=1; i<=7;i++){ 
					getline(*Idump,str1);
					stringstream ss(str1);
					if(i==4)ss >> simbox->xlo >> simbox->xhi >> simbox->xy; 
					if(i==5)ss >> simbox->ylo >> simbox->yhi >> simbox->xz; 
					if(i==6)ss >> simbox->zlo >> simbox->zhi >> simbox->yz; 
					if(i==7)dumpfields=str1;
				}
			}
		}
		simbox->theta = atan(simbox->xy/simbox->L);

		if(Idump->eof()){
			fprintf(log,"ERR : TIMESTEP %s NOT FOUND\n",gotostep.c_str());
			throw exception();
		}
	}
	simbox->tiltsign = (simbox->xy < 0.0) ? -1.0 : 1.0;
}

void FileHandle::StoreDump(string dumpstep){

	int id;
	string str1;
	ParseDump(dumpstep);
	int counter=0;

	while(counter < simbox->natoms){
		counter++;
		getline(*Idump,str1);
		stringstream ss(str1);
// 		ss >> id >>  atom->type[id-1] >>  atom->mol[id-1][n] >>  atom->r[0][id-1][n] >>  atom->r[1][id-1][n] 
// 		>> atom->r[2][id-1][n] >> atom->ix[id-1][n] >> atom->iy[id-1][n] >> atom->v[0][id-1][n] >> atom->v[1][id-1][n]
// 		>> atom->f[0][id-1][n] >> atom->f[1][id-1][n] >> atom->sxx[id-1][n] >>  atom->syy[id-1][n] >> atom->sxy[id-1][n]
// 		>> atom->pol1[id-1][n] >> atom->pe[id-1][n] >> atom->ke[id-1][n];
		
        ss >> id >>  atom->type[id-1] >>  atom->r[0][id-1][n] >>  atom->r[1][id-1][n] >> atom->r[2][id-1][n]
           >> atom->ix[id-1][n] >> atom->iy[id-1][n] >> atom->iz[id-1][n] >> atom->rad[id-1]
           >> atom->v[0][id-1][n] >> atom->v[1][id-1][n]>> atom->v[2][id-1][n]
           >> atom->sxx[id-1][n] >> atom->syy[id-1][n] >> atom->szz[id-1][n] >> atom->sxy[id-1][n] >> atom->sxz[id-1][n] >> atom->syz[id-1][n];
	}

}

void FileHandle::StoreData(){
// Step Time KinEng PotEng Pressure
	double d1;
	string str1;
	string comment = "#";
	int i=0;
	
	while(!Idata->eof()){
		getline(*Idata,str1);
		stringstream ss(str1);
		if(str1.compare(0,1,comment) != 0){
			ss >> thermo->step[i] >> d1 >> thermo->ke[0][i] >> thermo->pe[i] >> thermo->press[0][i] >> thermo->custom[0][i];
			i++;
		}
	}
	
	fprintf(log,"Datafile Stored\n");
}

void FileHandle::WriteDumpHeader(FILE* Outfile, int timestep){
	if(simbox->shearflag){
		fprintf(Outfile, "ITEM: TIMESTEP\n%d\nITEM: NUMBER OF ATOMS\n%d\n",timestep, simbox->natoms);
		fprintf(Outfile, "ITEM: BOX BOUNDS xy xz yz pp pp pp\n" 
						 "%.16g %.16g %.16g\n%.16g %.16g %.16g\n%.16g %.16g %.16g\n",
						 simbox->xlo,simbox->xhi,simbox->xy,
						 simbox->ylo,simbox->yhi,simbox->xz,
						 simbox->zlo,simbox->zhi,simbox->yz);
	}
	else{
		fprintf(Outfile,"ITEM: TIMESTEP\n%d\nITEM: NUMBER OF ATOMS\n%d\n",timestep, simbox->natoms);
		fprintf(Outfile,"ITEM: BOX BOUNDS pp pp pp\n" 
						 "%.16g %.16g\n%.16g %.16g\n%.16g %.16g\n",
						simbox->xlo,simbox->xhi,
                        simbox->ylo,simbox->yhi,
                        simbox->zlo,simbox->zhi);
	}

}


void FileHandle::WriteDump(const vector<int>&group, int timestep){
	WriteDumpHeader(Odump, timestep);
// 	fprintf(Odump,"%s\n", dumpfields.c_str());
	fprintf(Odump,"ITEM: ATOMS id x y z fx fy fz atompe\n");

	for(auto i:group){
// 		 i+1, atom->type[id], atom->mol[id][n], atom->r[0][id][n], atom->r[1][id][n],
// 		 atom->r[2][id][n], atom->ix[id][n], atom->iy[id][n], atom->v[0][id][n], atom->v[1][id][n],
// 		 atom->f[0][id][n], atom->f[0][id][n], atom->sxx[id][n],  atom->syy[id][n], atom->sxy[id][n],
// 		 atom->pol1[id][n], atom->pe[id][n], atom->ke[id][n];

		fprintf(Odump,"%d %.8f %.8f %.8f %1.5e %1.5e %1.5e\n",
						i+1, atom->r[0][i][n], atom->r[1][i][n], atom->r[2][i][n],
						  atom->f[0][i][n], atom->f[1][i][n], atom->f[2][i][n]);

	}
}

void FileHandle::Dump2Dat(string step, const vector<int>& group){
	FILE *fp;
	double strainval = stod(step)*simbox->dt*simbox->srate;
	char str[40];
	sprintf(str,"%.1f",strainval);
	string strain(str);
	string fname = "samples/VF0.68068-N97556-FCC-seed111-Rate0.0005-Minimize_SRate"+simbox->SRATE+"-Strain"+strain+".dat";
	fp = fopen(fname.c_str(),"w");

	// header
	fprintf(fp, "Strained Sample\n\n%d atoms\n%d atom types\n\n",simbox->natoms,simbox->natomtypes);
	fprintf(fp,	"%.16g %.16g xlo xhi\n%.16g %.16g ylo yhi\n%.16g %.16g zlo zhi\n",
				simbox->ylo,simbox->xhi,
                simbox->ylo,simbox->yhi,
                simbox->zlo,simbox->zhi);
	
	fprintf(fp, "\n%.16g %.16g %.16g xy xz yz\n",simbox->xy, simbox->xz, simbox->yz);

	fprintf(fp,"\nMasses\n\n1 1.0\n");
	
	// positions
	fprintf(fp,"\nAtoms\n\n");
	for(auto i:group)
		fprintf(fp,"%d %d %.16f %.16f %.16f\n",i+1,atom->type[i],atom->r[0][i][n], atom->r[1][i][n], atom->r[2][i][n]);
	
	//velocities 
	fprintf(fp,"\nVelocities\n\n");
	for(auto i:group)
		fprintf(fp,"%d %.16f %.16f %.16f\n",i+1,atom->v[0][i][n], atom->v[1][i][n], atom->v[2][i][n]);

	//radius
	fprintf(fp,"\nRadius\n\n");
	for(auto i:group)
		fprintf(fp,"%d %.16f\n",i+1,atom->rad[i]);
	
	fclose(fp);
}
