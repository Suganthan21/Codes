#ifndef DUMP_HANDLE_H
#define DUMP_HANDLE_H

#include "pointers.h"
#include<vector>
#include<string>
#include<fstream>

class FileHandle:Pointers
{
	public:
		FileHandle(Main*);
		~FileHandle();
		std::string filetype;
		std::vector<std::string> idump;
		std::vector<std::string> odump;
		std::vector<std::string> idata;
		std::vector<std::string> ststep;
		std::vector<int> itstep;
		int nstart, nevery, nstop;
		int ntimes;

		bool readdump, writedump;
		
		void Init();
		void StoreDump(std::string);
		void StoreData();
		void FileGood(std::string filename);
		void OpenClose(std::string,int);
		void WriteDump(const std::vector<int>&, int);
		void Dump2Dat(std::string,const std::vector<int>&);

	private:

		std::string dumpfields;

		std::ifstream* Idata;
		std::ifstream* Idump;
		FILE* Odump;

		void FileName();
		void DumpSteps(std::string);
		void ParseDump(std::string gotostep);
		void WriteDumpHeader(FILE*, int);
		void WriteDataHeader();
		void WriteData();
};

#endif
