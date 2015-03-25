
// The library need LAPACK and FFTW3, compiler option for these libraries are -llapack -lfftw3
#ifndef BRANCHES_ANALYSIS_INCLUDED
#define BRANCHES_ANALYSIS_INCLUDED

#include "polymer_lib.h"
#include <fftw3.h>

vector<double> all_distance_smoothed(vector<double> &distancematrix, double sigma);
double is_minimum(int m,int n,vector<double> &distancematrix);
vector<vector<pair<int,int> > > find_branches(vector<double> &distancematrix);
vector<vector<pair<int,int> > > shorten_branches(vector<vector<pair<int,int> > > & branches,vector<double> &x,vector<double> &y,vector<double> &z);





vector <double> all_distance_smoothed(vector<double> &distancematrix, double sigma)
{

  if(sigma<1e-8)//no smoothing
    return distancematrix;

  int N=(int)(sqrt(distancematrix.size())+0.5);

  double *input=&distancematrix[0];
  double *kernel;
  fftw_complex *input_transformed;
  fftw_complex *kernel_transformed;
  fftw_plan p_input,p_input_inv,p_kernel;
  kernel = (double*) fftw_malloc(sizeof(double) * N*N);
  input_transformed = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N*N);
  kernel_transformed = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N*N);
  p_input=fftw_plan_dft_r2c_2d(N,N, input, input_transformed,FFTW_ESTIMATE);
  p_kernel=fftw_plan_dft_r2c_2d(N,N, kernel, kernel_transformed,FFTW_ESTIMATE);
  p_input_inv= fftw_plan_dft_c2r_2d(N,N, input_transformed, input,FFTW_ESTIMATE);
  //kernel 1D
  vector<double> kernel_1D(N,0);
  double sum=0;
  //new version, no cutoff
  for(int x=0;x<=N/2;x++)
    {
      kernel_1D[x]=exp(-x*x/(2*sigma*sigma));
      sum+=kernel_1D[x];
    }
  for(int x=N/2;x<N;x++)
    {
      kernel_1D[x]=exp(-(N-x)*(N-x)/(2*sigma*sigma));
      sum+=kernel_1D[x];
    }
  for(int x=0;x<N;x++)
    {
      kernel_1D[x]/=sum*N;
    }
  //kernel
  for(int x=0;x<N;x++)
    for(int y=0;y<N;y++)
      kernel[x*N+y]=kernel_1D[x]*kernel_1D[y];

  fftw_execute(p_input); //tranform input
  fftw_execute(p_kernel);//tranform kernel
  

  //multiply (size of complex arrays=N*(N/2+1);
  for (int i = 0; i < N*(N/2+1); ++i)
      input_transformed[i] =  input_transformed[i]*kernel_transformed[i];

  fftw_execute(p_input_inv); //inverse transform input_transformed (=input_transformed*kernel_transformed)


  fftw_destroy_plan(p_input);
  fftw_destroy_plan(p_kernel);
  fftw_destroy_plan(p_input_inv);
  fftw_free(input_transformed);
  fftw_free(kernel);
  fftw_free(kernel_transformed);

  return distancematrix;
}



//return d which is negative if  distance between points x[n],y[n],z[n] and  x[m],y[m],z[m] is a local minima, positive otherwise.
double is_minimum(int m,int n,vector<double> & distancematrix)
{
  int N=(int)(sqrt(distancematrix.size())+0.5);
  m=(m+10*N)%N;
  n=(n+10*N)%N;
  double d0,d1,d2;
  d0=distancematrix[m*N+n];
  d1=distancematrix[((m+1)%N)*N + (n+1)%N];
  d2=distancematrix[((m-1+N)%N)*N + (n-1+N)%N];
  if(d0-d1>d0-d2)
    return d0-d1;
  return d0-d2;
}

//////////////////////////////////////////////////

vector<vector<pair<int,int> > > find_branches(vector<double> & distancematrix)
{
  vector<vector<pair<int,int> > > branches;
  vector<pair<int,int> > branche_tmp;

 ////////////////////////
 //find branches
 ////////////////////////
 int N=(int)(sqrt(distancematrix.size())+0.5);
 int pos1,pos2,pos1tmp,pos2tmp;
 double is_minimum_tmp,tmp;
 vector<int> delta1;
 vector<int> delta2;
 double writhe;
 delta1.push_back(-1); delta2.push_back(0);
 delta1.push_back(-1); delta2.push_back(1);
 delta1.push_back(0); delta2.push_back(1);
 for(int pos=0;pos<N;pos++)
 {
     //check if local min.
     is_minimum_tmp=is_minimum(pos-1,pos+1,distancematrix);
     if(is_minimum_tmp<1e-9)
     {
	   //follow_branche
	   branche_tmp.clear();
	   branche_tmp.push_back(make_pair(pos,pos));
	   pos1=pos-1;
	   pos2=pos+1;
	   while(is_minimum_tmp<1e-9 && abs(pos1-pos2)<N)
	     {
	       branche_tmp.push_back(make_pair(pos1,pos2));
	       is_minimum_tmp=1;
	       for(int i=0;i<delta1.size();i++)
	       {
		       tmp=is_minimum((pos1+delta1[i]+N)%N,(pos2+delta2[i]+N)%N,distancematrix);
		       if(tmp<is_minimum_tmp)
		       {
		         is_minimum_tmp=tmp;
		         pos1tmp=pos1+delta1[i];
		         pos2tmp=pos2+delta2[i];
		       }
	       }
	       pos1=pos1tmp;
	       pos2=pos2tmp;
	     }
	   branches.push_back(branche_tmp);
     }
 }


  return branches;
}


vector<vector<pair<int,int> > > shorten_branches(vector<vector<pair<int,int> > > & branches,vector<double> &x,vector<double> &y,vector<double> &z)
{
  vector<vector<pair<int,int> > > branches_kept=branches;
  //int N=(int)(sqrt(distancematrix.size())+0.5);
 int N=x.size();
 int p1,p2;
double dist;
 vector<double> distances;
 vector<double> distances_smoothed;
 for(int n=0;n<branches_kept.size();n++)
   {
     distances.resize(branches_kept[n].size());
     distances_smoothed.resize(branches_kept[n].size());
     //create distance array
     for(int p=0;p<branches_kept[n].size();p++)
       {
	 p1=branches_kept[n][p].first;
	 p2=branches_kept[n][p].second;
	 //dist=distancematrix[((p1+N)%N)*N+((p2+N)%N)];
	 dist=sqrt((x[(p1+N)%N]-x[(p2+N)%N])*(x[(p1+N)%N]-x[(p2+N)%N])+(y[(p1+N)%N]-y[(p2+N)%N])*(y[(p1+N)%N]-y[(p2+N)%N])+(z[(p1+N)%N]-z[(p2+N)%N])*(z[(p1+N)%N]-z[(p2+N)%N]));
	 distances[p]=dist;
       }	 
     //smooth distance array
     int radius=1;
     for(int p=0;p<branches_kept[n].size();p++)
       {
	 int r=radius;
	 if(r>p)r=p;
	 if(r>branches_kept[n].size()-p-1)r=branches_kept[n].size()-p-1;
	 double tmpsum=0,tmpcount=0;
	 for(int i=-r;i<=r;i++)
	   {
	     tmpsum+=distances[p+i];
	     tmpcount+=1;
	   }
	 distances_smoothed[p]=tmpsum/tmpcount;
       }	 
     //remove points from branches
     for(int p=branches_kept[n].size()-2;p>=0;p--)
       {
	 if(distances_smoothed[p]>distances_smoothed[p+1])
	   {
	     break;
	   }
	 branches_kept[n].pop_back();
       }	 

   }
  return branches_kept;
}


#endif

