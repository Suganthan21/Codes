#ifndef SIM_BOX_H
#define SIM_BOX_H

#include<string>
#include "pointers.h"

class SimBox:Pointers
{
	public: 
		SimBox(Main*);
		~SimBox(){}
		double dt;
		double L, xlo, ylo, zlo, xhi, yhi, zhi;
		double xy, xz, yz, theta;
		bool shearflag;
		int dim, natoms, natomtypes;
		int tiltsign;

		double pc,taup,f0, srate;
		std::string PC, TAUP, F0, SRATE;

		void BoxData(std::string, std::string);
		void Initialise(std::string, std::string);
		void Properties(const std::string);
		void RemapAtoms();
};

#endif

