#ifndef MEMORY_H 
#define MEMORY_H 

#include "pointers.h"
#include<type_traits>
#include<stdexcept>

class Memory:Pointers
{
	public:
		Memory(Main* );
		~Memory(){}

		template <typename T>
		T*** Create3DArray(int dim1, int dim2, int dim3)
		{
			if(dim1<=0 || dim2<=0 || dim3<=0){
				fprintf(log,"3D array creation - Invalid dimensions");
				throw std::exception();
			}

			T*** ptr = new T**[dim1];
			T** minipool = new T*[dim1*dim2];
			T* pool = new T[dim1*dim2*dim3];

			for (int i = 0; i < dim1; ++i, minipool += dim2){
				ptr[i] = minipool;
				for (int j = 0; j < dim2; ++j, pool += dim3){
					ptr[i][j] = pool;
				}
			}
			
			for (int i=0; i<dim1; i++)
				for(int j=0; j<dim2; j++)
					for(int k=0; k<dim3; k++)
						ptr[i][j][k] = T(0.0);

			if(std::is_same<T,double>::value)
				ndouble += dim1*dim2*dim3;
			if(std::is_same<T,int>::value)
				nint += dim1*dim2*dim3;
				
			return ptr;
		}


		template <typename T>
		T** Create2DArray(int dim1, int dim2)
		{
			if(dim1<=0 || dim2<=0 ){
				fprintf(log,"2D array creation - Invalid dimensions");
				throw std::exception();
			}

		   T** ptr = new T*[dim1];  
		   T* pool = new T[dim1*dim2];

		   for (int i = 0; i < dim1; ++i, pool += dim2 )
			   ptr[i] = pool;

		   for (int i=0; i<dim1; i++)
			   for(int j=0; j<dim2; j++)
				   ptr[i][j] = T(0.0);

		   if(std::is_same<T,double>::value)
			   ndouble += dim1*dim2;
		   if(std::is_same<T,int>::value)
			   nint += dim1*dim2;

		   return ptr;
		}

		template <typename T>
			T* Create1DArray(int dim1)
		{
			if(dim1<=0){
				fprintf(log,"1D array creation - Invalid dimensions");
				throw std::exception();
			}

		   T* ptr = new T[dim1];

		   for (int i=0; i<dim1; i++)
			   ptr[i] = T(0.0);

		   if(std::is_same<T,double>::value)
			   ndouble += dim1;
		   if(std::is_same<T,int>::value)
			   nint += dim1;

		   return ptr;
		}


		template <typename T>
		T*** Empty3DArray(T*** arr, int dim1, int dim2, int dim3)
		{
			for (int i=0; i<dim1; i++)
				for(int j=0; j<dim2; j++)
					for(int k=0; k<dim3; k++)
						arr[i][j][k] = T(0.0);
			return arr;
		}

		template <typename T>
		T** Empty2DArray(T** arr, int dim1, int dim2)
		{
			for (int i=0; i<dim1; i++)
				for(int j=0; j<dim2; j++)
					arr[i][j] = T(0.0);
			
			return arr;
		}

		template <typename T>
		T* Empty1DArray(T* arr, int dim1)
		{
			for (int i=0; i<dim1; i++)
				arr[i] = T(0.0);
			
			return arr;
		}

		template <typename T>
		void Delete3DArray(T*** arr)
		{
			delete [] arr[0][0];  
			delete [] arr[0];     
			delete [] arr;
			arr = nullptr;
		}

		template <typename T>
		void Delete2DArray(T** arr)
		{
			delete [] arr[0];
			delete [] arr;     
			arr = nullptr;
		}

		template <typename T>
		void Delete1DArray(T* arr)
		{
			delete[] arr;
			arr = nullptr;
		}
	void Allocation();

	private:
		int ndouble=0, nint=0;

};

#endif
