#include<iostream>
#include<stdio.h>
#include <stdlib.h>
#include<cstdlib>
#include<cmath>
#include <iomanip>
#include<fstream>
using namespace std;

int main()
{
	srand48(time(0));
	int n=256,cycles=70000;//lattice dimensions and number of mc cycles
	int i,j,k,l,a[n+2][n+2],E1,E2,eqE[40][cycles],up,down,m,avgE[40],eqcount,acc;

	// float temp[]={0.01, 0.05, 1, 1.2, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.1, 2.15, 2.2, 2.23, 2.235, 2.24, 2.245, 2.25, 2.255, 2.258, 2.26, 2.263, 2.265, 2.268, 2.27, 2.275, 2.28, 2.285, 2.3, 2.4, 2.5, 3, 4, 6};
	// float eqbmcyc[]={6.5, 6.5, 6.5, 8, 8, 8, 9, 9, 10, 10, 10.3, 10.5, 10.5, 10.5, 10.5, 10.5, 10.5, 10.5, 10.5, 11.5, 11.5, 11.5, 11.5, 11.5, 11.5, 11, 11, 11, 11, 9.5, 9, 9, 8, 8 , 8};

	float temp[]={1.5};
	float eqbmcyc[20]={10};


	long sumE, sumAcc;
	float M[40][cycles],sumM,avgM[40];
	float avgAcc[40];
	float x,y,T;
	int length = sizeof(temp)/sizeof(temp[0]);

	cout<<length<<"\n";
	for(k=0;k<length;k++){   //For different temperatures
		cout<<k<<"\n";

		//Defining the elements
		for(i=0;i<=n+1;i++){
			for(j=0;j<=n+1;j++){
				x=drand48();
				if(x<0.5){
					a[i][j]=+1;
				}
				else{
					a[i][j]=-1;
				}


			}
		}
		// To display the matrix 
		// for(i=0;i<=n+1;i++){
		//   for(j=0;j<=n+1;j++){
		//     cout<<"    "<<a[i][j];
		//     if (j==n+1){
		// 	cout<<"\n";
		//     }
		//   }
		// }

		//sumAcc=0;
		eqcount=0;
		sumM=0;
		sumE=0;
		T=temp[k];
		for(int l=0;l<cycles;l++){
			//acc=0;
			m=0;
			up=0;
			down=0;
			M[k][l]=0;
			eqE[k][l]=0;
			for(int i=1;i<=n;i++){
				for(int j=1;j<=n;j++){

					E1=-(a[i][j]*a[i][j+1]  +  a[i][j]*a[i+1][j]  +  a[i][j]*a[i-1][j] +  a[i][j]*a[i][j-1]);

					a[i][j]=-a[i][j];

					E2=-(a[i][j]*a[i][j+1]  +  a[i][j]*a[i+1][j]  +  a[i][j]*a[i-1][j] +  a[i][j]*a[i][j-1]);

					y=drand48();

					if(E2<=E1 || y<exp(-(E2-E1)/T) ){
						//if(l>exp(eqbmcyc[k])){
						// acc=acc+1;
						//}	  
						a[i][j]=a[i][j];
					}

					else{
						a[i][j]=-a[i][j];	
					}

					//Boundary conditions 
					if(i==1){
						a[n+1][j]=a[i][j];
					}
					if(j==1){
						a[i][n+1]=a[i][j];
					}

					if(i==n){
						a[0][j]=a[i][j];
					}      

					if(j==n){
						a[i][0]=a[i][j];
					}

				}   
			}

			//sumAcc=acc+sumAcc;


			if(l>exp(eqbmcyc[k])){
				eqcount=eqcount+1;
				for(int g=1;g<n;g++){
					for(int h=1;h<n;h++){
						eqE[k][l]=eqE[k][l]-a[g][h]*a[g][h+1]-a[g][h]*a[g+1][h];
						if(a[g][h]==1){
							up=up+1;
						}
						else{
							down=down+1;
						}

					}
				}
				m=abs(up-down);
				M[k][l]=(float)m/float(65536);
				sumM=sumM+M[k][l];
				sumE=sumE+eqE[k][l];
			}



		}//end of l

		avgM[k]=sumM/eqcount;
		avgE[k]=sumE/eqcount;

		//avgAcc[k]=(float)sumAcc/((float)65536*(float)eqcount);
		cout<<setprecision(8)<<fixed<<T<<"\t"<<avgE[k]<<"\t"<<avgM[k]<<"\n";
		// cout<<setprecision(8)<<fixed<<T<<"\t"<<avgM[k]<<"\n";


	}//end of k

	int files=2;
	string filename[files]={"E vs time.dat", "M vs time.dat"};

	ofstream outputFile[files];
	for (int i = 0; i < files; i++){
		outputFile[i].open(filename[i].c_str());
	}

	for (l=0 ; l<cycles ; l++ )  { //start from exp(eqbmcyc[k])
		outputFile[0] << setprecision(5) << fixed << log(l+1) << "\t";
		outputFile[1] << setprecision(5) << fixed << log(l+1) << "\t";
		for(k=0;k<length;k++){     
			outputFile[0] << setprecision(5) << fixed << eqE[k][l] << "\t";
			outputFile[1] << setprecision(5) << fixed << M[k][l] << "\t";
			if(k==length-1){
				outputFile[0]<<"\n";  outputFile[1] << "\n";
			}
		}
	}

	for (int i = 0; i < files; i++){  
		outputFile[i].close();
	}

	return 0;
}


