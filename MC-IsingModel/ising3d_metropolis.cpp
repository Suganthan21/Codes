#include<iostream>
#include<stdio.h>
#include <stdlib.h>
#include<cmath>
#include <iomanip>
#include<fstream>
using namespace std;

int main()
{
  srand48(time(0));
  int n=10,cycles=20000; //lattice dimensions and number of mc cycles
  int i,j,k,l,t,a[n+2][n+2][n+2],E1,E2,E[40][cycles],up,down,m;
  //float temp[]={0.01, 0.05, 1, 1.1, 1.2, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.1, 2.15, 2.2, 2.3, 2.4, 2.5, 3, 4, 6};
  // float temp[]={2.23, 2.235, 2.24, 2.245, 2.25, 2.255, 2.258, 2.26, 2.263, 2.265, 2.268, 2.27, 2.275, 2.28, 2.285};
  float temp[]={2,3,4,5,6};
  float x,y,T,M[40][cycles];
  int length = sizeof(temp)/sizeof(temp[0]);
  cout<<length<<"\n";

  for(t=0;t<length;t++){          //For different temperatures
    cout<<t<<"\n";
      
    for(i=0;i<=n+1;i++){
      for(j=0;j<=n+1;j++){
	for(k=0;k<=n+1;k++){
	  x=drand48();
	  if(x<0.5){
	    a[i][j][k]=+1;
	  }
	  else{
	    a[i][j][k]=-1;
	  }
	}	
      }
    }

    T=temp[k];
 
    for(int l=0;l<cycles;l++){
      cout<<l<<"\t";
      E[t][l]=0; m=0;      up=0;      down=0;      M[t][l]=0;
      
      for(int i=1;i<=n;i++){
	for(int j=1;j<=n;j++){
	  for(int k=1;k<=n;k++){
	    	    
	    //Metropolis
	  
	    E1 = -a[i][j][k]* ( a[i][j][k+1]  +  a[i][j][k-1]  +  a[i][j+1][k]  +  a[i][j-1][k]  +  a[i-1][j][k]  +  a[i+1][j][k] );

	    a[i][j][k] =-a[i][j][k];

	    E2 = -a[i][j][k]* ( a[i][j][k+1]  +  a[i][j][k-1]  +  a[i][j+1][k]  +  a[i][j-1][k]  +  a[i-1][j][k]  +  a[i+1][j][k] );
	  
	    y=drand48();
	  
	    // cout<<E2<<"\n";
	    // cout<<exp(-(E2-E1)/T)<<"\n\n\n";
	  
	    if(E2<=E1 ||y<exp(-(E2-E1)/T) ){
	      a[i][j][k] = a[i][j][k];
	    }
	      
	    else{
	      a[i][j][k] = -a[i][j][k];	
	    }

	    //Boundary conditions
	    if(j==1){
	      a[i][n+1][k] = a[i][j][k];
	    }
	    if(k==1){
	      a[i][j][n+1] = a[i][j][k];
	    }

	    if(j==n){
	      a[i][0][k] = a[i][j][k];
	    }      

	    if(k==n){
	      a[i][j][0] = a[i][j][k];
	    }
	    if(i==1){
	      a[n+1][j][k] = a[i][j][k];
	    }
	    if(i==n){
	      a[0][j][k] = a[i][j][k];
	    }
	    
	  }   
	}
      } //end of i
      
      for(int f=1;f<=n;f++){
	for(int g=1;g<=n;g++){
	  for(int h=1;h<=n;h++){
	  
	    E[t][l] = E[t][l] - a[f][g][h]* ( a[f][g][h+1]  +  a[f][g][h-1]  +  a[f][g+1][h]  +  a[f][g-1][h]  +  a[f+1][g][h]  +  a[f-1][g][h] );
	    
	    if(a[f][g][h]==1){
	      up=up+1;
	    }
	    else{
	      down=down+1;
	    }
	  }
	}
      }//end of f
      m = abs(up-down);
      M[t][l] = (float)m/float(n*n*n);
      cout<<E[t][l]<<"\t"<<M[t][l]<<"\t";
    }//end of l
  
  }//end of k


  int files=2;
  string filename[files]={"E vs time.dat", "M vs time.dat"};
 
  ofstream outputFile[files];
  for (int i = 0; i < files; i++){
    outputFile[i].open(filename[i].c_str());
  }

  for (l=0 ; l<cycles ; l++ )  {
    outputFile[0] << setprecision(5) << fixed << log(l+1) << "\t";
    outputFile[1] << setprecision(5) << fixed << log(l+1) << "\t";
    for(k=0;k<length;k++){     
      outputFile[0] << setprecision(5) << fixed << E[t][l] << "\t";
      outputFile[1] << setprecision(5) << fixed << M[t][l] << "\t";
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
    
      
