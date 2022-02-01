#include "memory.h"

Memory::Memory(Main* ptr):Pointers(ptr){}

void Memory::Allocation(){
	double double_size, int_size;
	double_size = ndouble*8*1e-06;
	int_size = nint*4*1e-06;

	fprintf(log,"Memory Allocation:\n");
	fprintf(log,"\t%d doubles = %g MB\n", ndouble, double_size);
	fprintf(log,"\t%d integers = %g MB\n", nint, int_size);
	fprintf(log,"\tTotal = %g MB\n", double_size+int_size);

}
