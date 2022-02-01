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
  int n = 256, cycles = 70000; //lattice dimensions and number of mc cycles
  
  int i,j,k,l,a[n+1][n+1], count, p, q, r, s, m, up, down ,dummy;
  
  //float temp[]={0.01, 0.05, 1, 1.2, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.1, 2.2, 2.23, 2.24, 2.245, 2.25, 2.255, 2.26, 2.265, 2.269, 2.27, 2.275, 2.28, 2.285, 2.29, 2.3, 2.35,     2.4, 2.45, 2.5,    2.6, 2.8, 3, 3.5, 4, 4.5, 5};

  float temp[] = { 2.3, 2.35,     2.4, 2.45, 2.5,    2.6, 2.8, 3};
  
  int length = sizeof(temp)/sizeof(temp[0]);
  float x,y1,y2,z,T,M[100][cycles], E[100][cycles];

  
  //Defining the elements

  for(k=0;k<length;k++){          //For different temperatures
      
    for(i=1;i<=n;i++){
      for(j=1;j<=n;j++){
	x=drand48();
	if(x<0.5){
	  a[i][j]=+1;
	}
	else{
	  a[i][j]=-1;
	}
      }
    }

    T=temp[k];
    cout<<T<<"\n";
    
    for(l = 0; l < cycles; l++){

      int  posx[(n+1)*(n+1)] = {0}, posy[(n+1)*(n+1)] = {0};
      count = 0;

      //Wolffy
      
      y1 = drand48(); y2 = drand48(); 
      p = abs(1 + (n-1)*y1); q = abs(1 + (n-1)*y2);
      posx[0] = p; posy[0] = q;
      
      dummy = a[p][q];
      a[p][q] = -a[p][q];
      
      for(r = 0; r <= count; r++){
	
	int forx[4] = {0, 0, 1, -1};
	int fory[4]= {1, -1, 0, 0};
  
	i = posx[r]; j = posy[r];
	
	int  neighbours[4] = { a[i][j+1], a[i][j-1], a[i+1][j], a[i-1][j] };

	//Boundary conditions
	if(i==1){
	  int  neighbours[4] = { a[i][j+1], a[i][j-1], a[i+1][j], a[n][j] };
	  forx[3] = n-1;
	}
	if(j==1){
	  int neighbours[4] = { a[i][j+1], a[i][n], a[i+1][j], a[i-1][j] };
	  fory[1] = n-1;
	}

	if(i==n){
	  int neighbours[4] = { a[i][j+1], a[i][j-1], a[1][j], a[i-1][j] };
	  forx[2] = 1-n;
	}      

	if(j==n){
	  int neighbours[4] = { a[i][1], a[i][j-1], a[i+1][j], a[i-1][j] };
	  fory[0] = 1-n;
	}
	
	for(s = 0; s < 4 ; s++){
	  z = drand48();
	  
	  if( neighbours[s] == dummy && z <= (1-exp(-2/T)) ){
	    
	    count = count + 1; 
	    posx[count] = i + forx[s]  ; posy[count] = j + fory[s];
	    a[posx[count]][posy[count]] = -dummy;
	    
	  }
	}
      } //end of r

      E[k][l] =  M[k][l] = 0.0;  m = up = down =0;

      for(int g = 1; g < n; g++){
	for(int h = 1; h < n; h++){
	  
	  E[k][l]=E[k][l]-a[g][h]*a[g][h+1]-a[g][h]*a[g+1][h];
	  
	  if(a[g][h]==1){
	    up += 1;
	  }
	  else{
	    down += 1;
	  }
	    
	}
      }
      m=abs(up-down);
      M[k][l]=(float)m/float((n-1)*(n-1));
    }//end of l

  }//end of k

  int files=2;
  string filename[files]={"E2 vs time.dat", "M2 vs time.dat"};
 
  ofstream outputFile[files];
  for (int i = 0; i < files; i++){
    outputFile[i].open(filename[i].c_str());
    outputFile[i] << "00" << "\t";
    for (k = 0; k < length; k++){
      outputFile[i]  << setprecision(3) << temp[k] << "\t";
      if(k==length-1){
	outputFile[i]<<"\n";
      }
    }
  }

  for (l=0 ; l<cycles ; l++ )  { //start from exp(eqbmcyc[k])
    outputFile[0] << setprecision(5) << fixed << log(l+1) << "\t";
    outputFile[1] << setprecision(5) << fixed << log(l+1) << "\t";
    
    for(k=0 ; k<length ; k++){     
      outputFile[0] << setprecision(1) << fixed << E[k][l] << "\t";
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
    
      
