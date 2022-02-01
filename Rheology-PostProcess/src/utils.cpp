#include<cmath>
#include "utils.h"
#include "filehandle.h"

using namespace std;

Utils::Utils(Main* ptr):Pointers(ptr){
}

Utils::~Utils(){}

vector<int> Utils::Sec2HHMMSS(double sec){
	int seconds, minutes, hours;
	seconds = round(sec);
	minutes = seconds / 60;
	hours = minutes / 60;
	vector<int> vec = {hours,(minutes%60),(seconds%60)};
	return vec;
}

void Utils::ChangeLogName(string iname, string oname){
	ifstream infile;
	infile.open(iname);
	filehandle->FileGood(iname);
	string str1;
	FILE *fp;
	fp = fopen(oname.c_str(),"w");
	while(!infile.eof()){
		getline(infile, str1);
		fprintf(fp,"%s\n",str1.c_str());
	}
	infile.close();
	fclose(fp);
}

string Utils::ToString(double value, int n){
	char arr[64];
	snprintf(arr,sizeof(arr),"%.*f", n, value);
	string str(arr);
	return str;
}
