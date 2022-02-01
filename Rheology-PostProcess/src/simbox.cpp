#include<regex>
#include<sstream>
#include<cmath>
#include "simbox.h"
#include "filehandle.h"
#include "particle.h"

using namespace std;

SimBox::SimBox(Main* ptr):Pointers(ptr){
	xy = 0.0;
	xz = 0.0;
	yz = 0.0;
	theta = 0.0;
	shearflag = 0;
}

void SimBox::BoxData(string filename, string ft){

	fprintf(log,"\nAcquiring simulation box settings\n");
	ifstream boxdata(filename);
	filehandle->FileGood(filename);
	string str1;

	if(ft=="dump" || ft=="dump,data"){
		string boxshear("ITEM: BOX BOUNDS xy xz yz pp pp pp");

		//First 9 lines contain all info in sim box
		for(int i=1;i<=9;i++){
			getline(boxdata, str1);
			stringstream ss(str1);
			if(i==4)natoms = stoi(str1);
			if(i==5)
				if(str1==boxshear)
					shearflag=1;
			if(i==6)ss >> xlo >> xhi >> xy; 
			if(i==7)ss >> ylo >> yhi >> xz; 
			if(i==8)ss >> zlo >> zhi >> yz; 
		}
	}
	else if(ft=="init"){

		regex r1("(.*) atoms");
		regex r2("(.*) atom types");
		regex r3("(.*) xlo xhi");
		regex r4("(.*) ylo yhi");
		regex r5("(.*) zlo zhi");
		regex r6("(.*) xy xz yz");

		//store header info
		getline(boxdata, str1);
		while (str1.find("Atoms") == string::npos){ 
			smatch match;
			//remove leading and trailing whitespace
			str1 = regex_replace(str1, regex("^ +| +$|( ) +"), "$1"); 
			const string s=str1;
			if (regex_search(s.begin(), s.end(), match, r1))
				natoms = stod(match[1]);
			if (regex_search(s.begin(), s.end(), match, r2))
				natomtypes = stod(match[1]);
			if (regex_search(s.begin(), s.end(), match, r3)){
				stringstream ss(match[1]);
				ss >> xlo >> xhi;
			}
			if (regex_search(s.begin(), s.end(), match, r4)){
				stringstream ss(match[1]);
				ss >> ylo >> yhi;
			}
			if (regex_search(s.begin(), s.end(), match, r5)){
				stringstream ss(match[1]);
				ss >> zlo >> zhi;
			}
			if (regex_search(s.begin(), s.end(), match, r6)){
				shearflag=1;
				stringstream ss(match[1]);
				ss >> xy >> xz >> yz;
			}
			getline(boxdata, str1);
			if(boxdata.eof()){
				fprintf(log,"ERROR: CHECK INIT FILE\n");
				throw exception();
			}
		}
	}
	else{
		fprintf(log,"ERROR: INVALID FILETYPE\n");
		throw exception();
	}

	L = yhi - ylo;
	theta = atan(xy/L);

	fprintf(log,"Number of atoms = %d\nBox Dimensions:\n"
	"\tX  %.16g\t%.16g \n\tY  %.16g\t%.16g   \n\tZ  %.16g\t%.16g\n",
	natoms, xlo, xhi, ylo, yhi, zlo, zhi);
	fprintf(log,"Length of box = %g\n", L);
	fprintf(log,"Box tilt factors:\n"
	"\tXY = %g, XZ = %g, YZ = %g\n", xy, xz, yz);

	boxdata.close();
}

void SimBox::Initialise(string filename,string props){
	fprintf(log, "\nReading initial configuration...\n");
	stringstream ss(props);
	string tmp, str1;
	vector<string> substring;

	while(getline(ss,tmp,','))
		substring.push_back(tmp);

	ifstream initfile(filename);
	int id;

	for(auto word:substring){
		fprintf(log,"Storing %s info...",word.c_str());
		while (str1.find(word) == string::npos){
			getline(initfile, str1);
			if(initfile.eof()){
				fprintf(log,"ERROR: CHECK INIT FILE. ATTRIBUTE %s NOT FOUND\n", word.c_str());
				throw exception();
			}
		}
		getline(initfile, str1);

		if(word=="Atoms"){
			for(int i=0; i<natoms; i++){
				getline(initfile,str1);
				stringstream ss(str1);
				ss >> id >>  atom->type[id-1] >>  atom->r[0][id-1][n] >>  atom->r[1][id-1][n] >> atom->r[2][id-1][n];
			}
			fprintf(log,"done\n");
		}
		else if(word=="Velocities"){
			for(int i=0; i<natoms; i++){
				getline(initfile,str1);
				stringstream ss(str1);
				ss >> id >> atom->v[0][id-1][n] >>  atom->v[1][id-1][n] >> atom->v[2][id-1][n];
			}
			fprintf(log,"done\n");
		}
		else if(word=="Radius"){
			for(int i=0; i<natoms; i++){
				getline(initfile,str1);
				stringstream ss(str1);
				ss >> id >> atom->rad[id-1];
			}
			fprintf(log,"done\n");
		}
		else if(word=="Masses"){
			double mass;
			for(int j=0; j<natomtypes; j++){
				getline(initfile,str1);
				stringstream ss(str1);
				ss >> id >> mass; //here id denotes type id
				for(int i=0; i<natoms; i++){
					if(atom->type[i] == id)
						atom->mass[i] = mass;
				}
			}
			fprintf(log,"done\n");
		}
		else{
			fprintf(log,"\nERROR: CHECK INIT FILE. ATTRIBUTE %s NOT DEFINED\n", word.c_str());
			throw exception();
		}
		initfile.seekg(0);
	}
	
		
}

void SimBox::Properties(const string s){

	regex rgx(".*sp_taup(.*)_f(.*)\\..*");
    smatch match;
    if (regex_search(s.begin(), s.end(), match, rgx)){
		TAUP = match[1]; 
		F0 = match[2]; 

		taup = stod(match[1]);
		f0 =stod(match[2]);
// 		fprintf(log,"Active Particles = %d\t, TauP = %g\t, f0=%g\n",nactive,tauP,f0);
	}
}



void SimBox::RemapAtoms(){
	double excess;
	
	double*** r = atom->r;
	double L = simbox->L;

	for(int i=1; i<=natoms; i++){
		excess = abs(tan(theta)*r[1][i][n]);
		r[0][i][n] += -tiltsign*excess;
	}
	for(int i=1; i<=natoms; i++){
		if(r[0][i][n] < 0) r[0][i][n] += L;
		if(r[0][i][n] > L) r[0][i][n] -= L;
	}

}
