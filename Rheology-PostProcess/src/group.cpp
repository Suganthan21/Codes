#include<sstream>
#include "particle.h"
#include "simbox.h"
#include "group.h"

using std::string;

Group::Group(Main* ptr):Pointers(ptr){}

void Group::SelectAtoms(string groups){
	fprintf(log, "\nMaking atom groups\n");
	all.clear();
	nselect.clear();
	active.clear();
	passive.clear();

	std::stringstream ss(groups);
	string tmp;
	std::vector<string> substring;

	while(getline(ss,tmp,','))
		substring.push_back(tmp);

	int ngroup;

	for(auto word:substring){
		if(word == "all"){
			for(int i=0; i<simbox->natoms; i++)
				all.push_back(i);
			nall = all.size();
			ngroup = nall;
		}
		else if(word == "active"){
			for(int i=0; i <simbox->natoms; i++)
				if(atom->mol[i] == 1)
					active.push_back(i);
			nactive = active.size();
			ngroup = nactive;
		}
		else if(word=="passive"){
			for(int i=0; i <simbox->natoms; i++)
				if(atom->mol[i] == 0)
					passive.push_back(i);
			npassive = passive.size();
			ngroup = npassive;
		}
		else{
			fprintf(log,"ERROR: UNKNOWN GROUP\n");
			throw std::exception();
		}
		fprintf(log,"%d atoms in group \"%s\"\n",ngroup, word.c_str());
	}

	
}
