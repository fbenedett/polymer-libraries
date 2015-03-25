// POLYMER_LIB.H  Routines made by Fabrizio Benedetti, Julien Dorier, Yannis Burnier.  University of Lausanne.
// 2012 - NOW
// The routine collected are used to analyze several observables of brownian dynamic
// simulations of polymers.  The routine came with NO WARRANTY! Use at your own risk.
// The library need LAPACK, compiler option for lapack is -llapack 

#ifndef POLYMER_LIB_INCLUDED
#define POLYMER_LIB_INCLUDED

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <math.h>
#include <time.h>
#include <functional>
#include <numeric>
#include <vector>
#include <sstream>
#include "Tstatistic.hh"
#include <sys/time.h>
#include <complex.h>


using namespace std;

double pitagora(double a, double b, double c);
double Bound(double val, double min,double max);
double writhe_of_frame(double *chain_x, double *chain_y, double *chain_z, long NbSegments);
double twist_of_frame(double *chain_x, double *chain_y, double *chain_z, double *ph_x, double *ph_y, double *ph_z, long NbSegments);
double twist_of_frame(double *chain_x, double *chain_y, double *chain_z, double *ph_x, double *ph_y, double *ph_z, double *th_x, double *th_y, double *th_z, long NbSegments);
double getx(long current_frame, long num_atoms, long atom_num1, double *alval);
double gety(long current_frame, long num_atoms, long atom_num1, double *alval);
double getz(long current_frame, long num_atoms, long atom_num1, double *alval);
void get_vector(long atom_num, long num_atoms, long current_frame, double *res, double *alval);
double norm_scalarp(long atom_num1, long atom_num2, long num_atoms, long current_frame, double *alval);
double scalar_prod(long atom_num1, long atom_num2, long num_atoms, long current_frame, double *alval);
double scalar_prod(double x, double y, double z, double x1, double y1, double z1, bool norm);
double calc_bending_en(long atom_num1, long atom_num2, long num_atoms, long current_frame, double *alval);
void frames_atoms_size(string filename, long &f_num, long &a_num);
void load_matrix(double *data, long num_atoms, long num_frames);
double error_bar(double *data, long size);
double distance(long current_frame, long num_atoms, long atom_num1, long atom_num2, double *alval);
void center_mass_xyz(long current_frame, long num_atoms, long start, long end, double *res, double *alval);
double ave_frame_distace(long num_frames, long num_atoms, long index_atom_reference, long first_index_interested, long last_index_interested, double *alval, double *results);
double radius_gyration(long current_frame, long num_atoms, double *alval);
double radius_gyration(long current_frame, long num_atoms, double *alval, long start, long end);
double end_to_end_distance(vector <double> &x, vector <double> &y, vector <double> &z, int window, bool circular);
void remap_part(double *chainx, double *chainy, double *chainz, long natoms, double sizebox_x, double sizebox_y, double sizebox_z);
void remap_noboundary(double *data, long num_frames, long num_atoms);
void remap_noboundary(double *data, long num_frames, long num_atoms, double sizebox_x, double sizebox_y, double sizebox_z);
double threeprod(double x1, double x2, double x3, double y1, double y2, double y3, double z1, double z2, double z3);
void crossprod(double x, double y, double z, double x1, double y1, double z1, double *res);
void observable_tstat(vector <vector <double> > &all_data, vector <string> nameobs, long current_f);
void write_CM(long start, long end, long num_frames, long num_atoms, double epsilon, double *data);

void check_periodic_contact(double *chainx, double *chainy, double *chainz, long natoms, long num_frames, double sizebox_x, double sizebox_y, double sizebox_z, int pos1, int pos2, long &real_contact_count, long &periodic_contact_count, double epsilon);
vector<double>  F(vector<double> &x,vector<double> &y,vector<double> &z);

double get_writhe(int pos1,int pos2,vector<double> &x,vector<double> &y,vector<double> &z);
void get_writhe_new(vector<int> half_windows_sizes,const vector<double> &x,const vector<double> &y,const vector<double> &z,vector<vector<double> > & writhe);
double get_twist(int pos1, int pos2, vector <double> &chain_x, vector <double> &chain_y, vector <double> &chain_z, vector <double> &ph_x, vector <double> &ph_y, vector <double> &ph_z, vector <double> &th_x, vector <double> &th_y, vector <double> &th_z);
double get_linking(vector<double> &x,vector<double> &y,vector<double> &z, vector<double> &xs,vector<double> &ys,vector<double> &zs);
void average(vector <double> &avect, double &avev);
double median(vector <double> &a);
void average_sigma(vector <double> &avect, double &average, double&sigma);
void average_bead(vector <double> &matrix_d, int window, double &average, double &sigma);
vector<vector<pair<int,int> > > shorten_branches(vector<vector<pair<int,int> > > & branches,vector<double> &x,vector<double> &y,vector<double> &z);


/* DSYEV prototype */
extern "C" void dsyev_( char* jobz, char* uplo, long int* n, double* a, long int* lda,
                double* w, double* work, long int* lwork, long int* info );

int prol_asph(vector<double> & x,vector<double> & y,vector<double> & z, double &prolactness, double &asphericity, double &sq_radius_gyr);


// these routine assume that we have several frames and that the coordinate
// are stored in a single, one-dimensional, array. So, to correctly get
// the x y z coordinate of a certain atom we have to play with the indexes
// if we are interested to the x coordinate of atom "atom_num1" (where the atom_num1
// is the index) we will have:
// "current_frame*num_atoms*3 + atom_num1*3"



///////////////////////////////////////////////////////////////////////////////////////////////////

// convert a number to a string
template <class T>
string num_to_str(T a)
{
   ostringstream sstream;
   sstream<<a;
    return sstream.str();
}

// pitagora formula
double pitagora(double a, double b, double c)
{
  double v= pow(a*a + b*b +c*c, 0.5);
  return v;
}


double Bound(double val, double min,double max){return (val<min)?min:((val>max)?max:val);}

// cross product between vector (x1,x2,x3) and (y1,y2,y3) then multiplied by (z1,z2,z3)
double threeprod(double x1, double x2, double x3, double y1, double y2, double y3, double z1, double z2, double z3)
{return ((x2*y3-x3*y2)*z1+(x3*y1-x1*y3)*z2+(x1*y2-x2*y1)*z3);}


void crossprod(double x, double y, double z, double x1, double y1, double z1, double *res){
res[0]=y*z1 - z*y1;
res[1]=z*x1 - x*z1;
res[2]=x*y1 - y*x1;
}


/*Using auxiliary beads it calculate the twist of a given chain.
Chain_x _y _z contain the coordinates of the real beads of the polymer.  ph_x _y _z contain the auxiliary beads.
This routine should be applied only when the coordinate frame is overlapped with a real bead, i.e. one virtual bead for every real bead
*/
double twist_of_frame(double *chain_x, double *chain_y, double *chain_z, double *ph_x, double *ph_y, double *ph_z, long NbSegments)
  {//needs the portion of the chain on which we want to calculate the twist in chain_x,y,z 
   //and the corresponding portion of phanom chain ph_x,y,z!
   //sum the dihedral angles along the chain and provide the total twist of the chain.
   //useful together with the writhe to check the total linking number as twist+writhe = DLk
  double twist=0;
  double r1x,r1y,r1z,r2x,r2y,r2z;
  double e1x,e2x,e1y,e2y,e1z,e2z;
  double nx,ny,nz;
  double nn,nr1,nr2,ne1,ne2;
  double see,se1,se2,s,s1,s2;
  double alpha, a1,a2;
  double Pi=M_PI;

  for(long i=0;i<NbSegments;i++)
	{
	//first segment
	r1x=chain_x[(i+1)%NbSegments]-chain_x[i];
	r1y=chain_y[(i+1)%NbSegments]-chain_y[i];
	r1z=chain_z[(i+1)%NbSegments]-chain_z[i];
	//normalize
	nr1=pitagora(r1x,r1y,r1z);
//  if(nr1<1e-6)cerr<<"nr1="<<nr1<<" i="<<i<<endl; 
	r1x=r1x/nr1; r1y=r1y/nr1; r1z=r1z/nr1;
	//second segment
	r2x=chain_x[(i+2)%NbSegments]-chain_x[(i+1)%NbSegments];
	r2y=chain_y[(i+2)%NbSegments]-chain_y[(i+1)%NbSegments];
	r2z=chain_z[(i+2)%NbSegments]-chain_z[(i+1)%NbSegments];
	//normalize
	nr2=pitagora(r2x,r2y,r2z);
//  if(nr2<1e-6)cerr<<"nr2="<<nr2<<" i="<<i<<endl; 
	r2x=r2x/nr2; r2y=r2y/nr2; r2z=r2z/nr2;
	//dihedral vectors
	e1x=-chain_x[i]+ph_x[i];
	e1y=-chain_y[i]+ph_y[i];
	e1z=-chain_z[i]+ph_z[i];
	e2x=-chain_x[(i+1)%NbSegments]+ph_x[(i+1)%NbSegments];
	e2y=-chain_y[(i+1)%NbSegments]+ph_y[(i+1)%NbSegments];
	e2z=-chain_z[(i+1)%NbSegments]+ph_z[(i+1)%NbSegments];
	//normalize
	ne1=pitagora(e1x,e1y,e1z);
	ne2=pitagora(e2x,e2y,e2z);
        //  if(ne1<1e-6)cerr<<"ne1="<<ne1<<" i="<<i<<endl; 
        //  if(ne2<1e-6)cerr<<"ne2="<<ne2<<" i="<<i<<endl; 

	e1x=e1x/ne1; e1y=e1y/ne1; e1z=e1z/ne1;
	e2x=e2x/ne2; e2y=e2y/ne2; e2z=e2z/ne2;
	//normal vector to the plane containing the two segments
	nx=r1y*r2z-r1z*r2y;
	ny=r1z*r2x-r1x*r2z;
	nz=r1x*r2y-r1y*r2x;
	nn=pitagora(nx,ny,nz);
        //cout<<ne1<<" "<<ne2<<" nn "<<nn<<endl;
	if(nn<0.0000000001)
		{//vectors parallel, see=scalar product, s=orientation
		see=e1x*e2x+e1y*e2y+e1z*e2z;
		s=threeprod(e1x,e1y,e1z,e2x,e2y,e2z,r1x,r1y,r1z);
 //   if(fabs(s)<1e-6)cerr<<"s="<<s<<" i="<<i<<endl; 
  //  if(fabs(see)>1-1e-6)cerr<<"see="<<see<<" i="<<i<<endl; 

		alpha=acos(see)*s/fabs(s);//alpha=fabs(fmod(acos(see),Pi))*s/fabs(s);
                cout<<see<<" normal vector was 0 "<<s<<" "<<alpha<<endl;
		}
	else
		{
		nx=nx/nn;ny=ny/nn;nz=nz/nn;
		se1=e1x*nx+e1y*ny+e1z*nz;
		se2=e2x*nx+e2y*ny+e2z*nz;
		s1=threeprod(e1x,e1y,e1z,nx,ny,nz,r1x,r1y,r1z);
    //special case, e1 is parallel to n, then the angle is 0 or Pi
    if(fabs(s1)<1e-10)	{a1=0; cout<<"normal vector parallel to a1-0"<<endl;
			if(e1x*nx+e1y*ny+e1z*nz<0){a1=Pi; //cout<<"normal vector parallel to a1-Pi"<<endl;
             }
		      	}
    //'general' case:
    else a1=acos(se1)*s1/fabs(s1);//a1=fabs(fmod(acos(se1),Pi))*s1/fabs(s1);

		s2=threeprod(nx,ny,nz,e2x,e2y,e2z,r1x,r1y,r1z);//r2
 		if(fabs(s2)<1e-10) {a2=0;
			if(e2x*nx+e2y*ny+e2z*nz<0){a2=Pi; //cout<<"normal vector parallel to a2-Pi"<<endl;
       }
      }
		else a2=acos(se2)*s2/fabs(s2);//a2=fabs(fmod(acos(se2),Pi))*s2/fabs(s2);

       		//    if(fabs(s1)<1e-10)cerr<<"s1="<<s1<<" i="<<i<<" a1="<<a1<<endl; 
 		//   if(fabs(se1)>1-1e-10)cerr<<"se1="<<se1<<" i="<<i<<" a1="<<a1<<endl; 
		//   if(fabs(s2)<1e-10)cerr<<"s2="<<s2<<" i="<<i<<" a2="<<a2<<endl; 
 		//   if(fabs(se2)>1-1e-10)cerr<<"se2="<<se2<<" i="<<i<<" a2="<<a2<<endl; 
		alpha=a1+a2;
		//cout<<se1<<" "<<se2<<" nn was not 0 "<<s1<<" "<<s2<<" "<<alpha<<endl;
		}
	// rescale the angle so that it is between -Pi and Pi
	if(alpha>Pi){alpha=-2*Pi+alpha;}
	if(alpha<-Pi){alpha=2*Pi+alpha;}

	twist=twist+alpha;
	// cout<<a1<<" a12 " <<a2<<" alpha " <<alpha<<" twist "<<twist/2/Pi<<endl;
	}
  return twist/(2.0*Pi);
  }


/*Using auxiliary beads it calculate the twist of a given chain.
Chain_x _y _z contain the coordinates of the real beads of the polymer.  ph_x _y _z contain the middle auxiliary beads.
th_x _y _z cointain the coordinates of the out of axis beads.
This routine should be applied only when the coordinate frame is overlapped with a real bead, i.e. one virtual bead for every real bead
*/
double twist_of_frame(double *chain_x, double *chain_y, double *chain_z, double *ph_x, double *ph_y, double *ph_z, double *th_x, double *th_y, double *th_z, long NbSegments)
  {//needs the portion of the chain on which we want to calculate the twist in chain_x,y,z 
   //and the corresponding portion of phanom chain ph_x,y,z!
   //sum the dihedral angles along the chain and provide the total twist of the chain.
   //useful together with the writhe to check the total linking number as twist+writhe = DLk
  double twist=0;
  double r1x,r1y,r1z,r2x,r2y,r2z, rp1x, rp1y, rp1z, rp2x, rp2y, rp2z;
  double e1x,e2x,e1y,e2y,e1z,e2z;
  double nx,ny,nz;
  double prod;
  double nn,nr1,nr2,ne1,ne2, np1, np2;
  double see,se1,se2,s,s1,s2;
  double alpha, a1,a2;
  double Pi=M_PI;
  
  for(long i=0;i<NbSegments;i++)
	{
	//first segment
	r1x=chain_x[(i+1)%NbSegments]-chain_x[i];
	r1y=chain_y[(i+1)%NbSegments]-chain_y[i];
	r1z=chain_z[(i+1)%NbSegments]-chain_z[i];
	//normalize
	nr1=pitagora(r1x,r1y,r1z);
//  if(nr1<1e-6)cerr<<"nr1="<<nr1<<" i="<<i<<endl; 
	r1x=r1x/nr1; r1y=r1y/nr1; r1z=r1z/nr1;

//  cout<<"r1x r1y r1z"<<r1x<<" "<<r1y<<" "<<r1z<<endl;

	//second segment
	r2x=chain_x[(i+2)%NbSegments]-chain_x[(i+1)%NbSegments];
	r2y=chain_y[(i+2)%NbSegments]-chain_y[(i+1)%NbSegments];
	r2z=chain_z[(i+2)%NbSegments]-chain_z[(i+1)%NbSegments];
	//normalize
	nr2=pitagora(r2x,r2y,r2z);
//  if(nr2<1e-6)cerr<<"nr2="<<nr2<<" i="<<i<<endl; 
	r2x=r2x/nr2; r2y=r2y/nr2; r2z=r2z/nr2;
//  cout<<"r2x r2y r2z"<<r2x<<" "<<r2y<<" "<<r2z<<endl;
	//dihedral vectors
	e1x=-ph_x[i]+th_x[i];
	e1y=-ph_y[i]+th_y[i];
	e1z=-ph_z[i]+th_z[i];
//  cout<<" th_x value  "<<th_x[i]<<endl;
//  cout<<" ph_x value  "<<ph_x[i]<<endl;
//  cout<<" E1X value  "<<e1x<<" E1Y value  "<<e1y<<" E1z value  "<<e1z<<" i  "<<i<<endl;
	e2x=-ph_x[(i+1)%NbSegments]+th_x[(i+1)%NbSegments];
	e2y=-ph_y[(i+1)%NbSegments]+th_y[(i+1)%NbSegments];
	e2z=-ph_z[(i+1)%NbSegments]+th_z[(i+1)%NbSegments];
	//project e1 on the plane normal to r1
	prod=e1x*r1x+e1y*r1y+e1z*r1z;//scalar prod e1*r1
	e1x=e1x-prod*r1x;
	e1y=e1y-prod*r1y;
	e1z=e1z-prod*r1z;
	//project e2 on the plane normal to r2
	prod=e2x*r2x+e2y*r2y+e2z*r2z;//scalar prod e2*r2
	e2x=e2x-prod*r2x;
	e2y=e2y-prod*r2y;
	e2z=e2z-prod*r2z;

	//normalize
	ne1=pitagora(e1x,e1y,e1z);
	ne2=pitagora(e2x,e2y,e2z);
        //  if(ne1<1e-6)cerr<<"ne1="<<ne1<<" i="<<i<<endl; 
        //  if(ne2<1e-6)cerr<<"ne2="<<ne2<<" i="<<i<<endl; 

	e1x=e1x/ne1; e1y=e1y/ne1; e1z=e1z/ne1;
	e2x=e2x/ne2; e2y=e2y/ne2; e2z=e2z/ne2;
	//normal vector to the plane containing the two segments
	nx=r1y*r2z-r1z*r2y;
	ny=r1z*r2x-r1x*r2z;
	nz=r1x*r2y-r1y*r2x;
	nn=pitagora(nx,ny,nz);
//  cout<<"ne1 ne2 nn:"<<ne1<<" "<<ne2<<" "<<nn<<endl;


//////////////////  calculation of dihedral angles
	if(nn<0.0000000001)
		{//vectors parallel, see=scalar product, s=orientation
		see=e1x*e2x+e1y*e2y+e1z*e2z;
		s=threeprod(e1x,e1y,e1z,e2x,e2y,e2z,r1x,r1y,r1z);
 //   if(fabs(s)<1e-6)cerr<<"s="<<s<<" i="<<i<<endl; 
 //   if(fabs(see)>1-1e-6)cerr<<"see="<<see<<" i="<<i<<endl; 

		alpha=acos(see)*s/fabs(s);//alpha=fabs(fmod(acos(see),Pi))*s/fabs(s);
                cout<<see<<" normal vector was 0 "<<s<<" "<<alpha<<endl;
		}
	else
		{
		  nx=nx/nn;ny=ny/nn;nz=nz/nn;
		  se1=e1x*nx+e1y*ny+e1z*nz;
		  se2=e2x*nx+e2y*ny+e2z*nz;
		  s1=threeprod(e1x,e1y,e1z,nx,ny,nz,r1x,r1y,r1z);
		  //special case, e1 is parallel to n, then the angle is 0 or Pi
		  if(fabs(s1)<1e-10)	{a1=0; cout<<"normal vector parallel to a1-0"<<endl;
		    if(e1x*nx+e1y*ny+e1z*nz<0){a1=Pi; //cout<<"normal vector parallel to a1-Pi"<<endl;
		    }
		  }
		  //'general' case:
		  else a1=acos(se1)*s1/fabs(s1);//a1=fabs(fmod(acos(se1),Pi))*s1/fabs(s1);
		  
		  s2=threeprod(nx,ny,nz,e2x,e2y,e2z,r2x,r2y,r2z);//r2
		  if(fabs(s2)<1e-10) {a2=0;
		    if(e2x*nx+e2y*ny+e2z*nz<0){a2=Pi; //cout<<"normal vector parallel to a2-Pi"<<endl;
		    }
		  }
		  else a2=acos(se2)*s2/fabs(s2);//a2=fabs(fmod(acos(se2),Pi))*s2/fabs(s2);
		  
		  //    if(fabs(s1)<1e-10)cerr<<"s1="<<s1<<" i="<<i<<" a1="<<a1<<endl; 
		  // 		if(fabs(se1)>1-1e-10)cerr<<"se1="<<se1<<" i="<<i<<" a1="<<a1<<endl; 
		  //		if(fabs(s2)<1e-10)cerr<<"s2="<<s2<<" i="<<i<<" a2="<<a2<<endl; 
		  // 		if(fabs(se2)>1-1e-10)cerr<<"se2="<<se2<<" i="<<i<<" a2="<<a2<<endl; 
		  alpha=a1+a2;
		  //		cout<<se1<<" "<<se2<<" nn was not 0 "<<s1<<" "<<s2<<" "<<alpha<<endl;
		}
	// rescale the angle so that it is between -Pi and Pi
	if(alpha>Pi){alpha=-2*Pi+alpha;}
	if(alpha<-Pi){alpha=2*Pi+alpha;}
	
   

	twist=twist+alpha;
	// cout<<a1<<" a12 " <<a2<<" alpha " <<alpha<<" twist "<<twist/2/Pi<<endl;
	}

  return twist/(2.0*Pi);
  }




double writhe_of_frame(double *chain_x, double *chain_y, double *chain_z, long NbSegments)
{
double Writhe; //, CrossingNumber;
 Writhe=0;//CrossingNumber=0;
  double ZERO=1e-10;
      double norme;
      double r13x,r13y,r13z;
      double r14x,r14y,r14z;
      double r24x,r24y,r24z;
      double r23x,r23y,r23z;
      double r34x,r34y,r34z;
      double r12x,r12y,r12z;
      double n1x,n1y,n1z;   
      double n2x,n2y,n2z;   
      double n3x,n3y,n3z;   
      double n4x,n4y,n4z;   
      double n5x,n5y,n5z;   
      //The writhe is a double integral on the chain path
  for(long i=1;i<NbSegments;i++)
	{
	  for(long j=0;j<i;j++)
	    {
	      r13x=chain_x[j]-chain_x[i];r13y=chain_y[j]-chain_y[i];r13z=chain_z[j]-chain_z[i];
	      r14x=chain_x[(j+1)%NbSegments]-chain_x[i];r14y=chain_y[(j+1)%NbSegments]-chain_y[i];r14z=chain_z[(j+1)%NbSegments]-chain_z[i];
	      r24x=chain_x[(j+1)%NbSegments]-chain_x[(i+1)%NbSegments];r24y=chain_y[(j+1)%NbSegments]-chain_y[(i+1)%NbSegments];r24z=chain_z[(j+1)%NbSegments]-chain_z[(i+1)%NbSegments];
	      r23x=chain_x[j]-chain_x[(i+1)%NbSegments];r23y=chain_y[j]-chain_y[(i+1)%NbSegments];r23z=chain_z[j]-chain_z[(i+1)%NbSegments];
	      r34x=chain_x[(j+1)%NbSegments]-chain_x[j];r34y=chain_y[(j+1)%NbSegments]-chain_y[j];r34z=chain_z[(j+1)%NbSegments]-chain_z[j];
	      r12x=chain_x[(i+1)%NbSegments]-chain_x[i];r12y=chain_y[(i+1)%NbSegments]-chain_y[i];r12z=chain_z[(i+1)%NbSegments]-chain_z[i];
	      n1x=r13y*r14z-r13z*r14y; n1y=r13z*r14x-r13x*r14z; n1z=r13x*r14y-r13y*r14x;
	      norme=sqrt(n1x*n1x+n1y*n1y+n1z*n1z);
	      if(norme>ZERO){n1x=n1x/norme; n1y=n1y/norme; n1z=n1z/norme;}
	      
	      n2x=r14y*r24z-r14z*r24y; n2y=r14z*r24x-r14x*r24z; n2z=r14x*r24y-r14y*r24x;
	      norme=sqrt(n2x*n2x+n2y*n2y+n2z*n2z);
	      if(norme>ZERO){n2x=n2x/norme; n2y=n2y/norme; n2z=n2z/norme;}
	      
	      n3x=r24y*r23z-r24z*r23y; n3y=r24z*r23x-r24x*r23z; n3z=r24x*r23y-r24y*r23x;
	      norme=sqrt(n3x*n3x+n3y*n3y+n3z*n3z);
	      if(norme>ZERO){n3x=n3x/norme; n3y=n3y/norme; n3z=n3z/norme;}
	      
	      n4x=r23y*r13z-r23z*r13y; n4y=r23z*r13x-r23x*r13z; n4z=r23x*r13y-r23y*r13x;
	      norme=sqrt(n4x*n4x+n4y*n4y+n4z*n4z);
	      if(norme>ZERO){n4x=n4x/norme; n4y=n4y/norme; n4z=n4z/norme;}
	      
	      n5x=r34y*r12z-r34z*r12y; n5y=r34z*r12x-r34x*r12z; n5z=r34x*r12y-r34y*r12x;
	      if(n5x*r13x+n5y*r13y+n5z*r13z > ZERO)
		{
		  Writhe += (asin(Bound(n1x*n2x+n1y*n2y+n1z*n2z,-1,1))+asin(Bound(n2x*n3x+n2y*n3y+n2z*n3z,-1,1))+asin(Bound(n3x*n4x+n3y*n4y+n3z*n4z,-1,1))+asin(Bound(n4x*n1x+n4y*n1y+n4z*n1z,-1,1)))/(2*M_PI);
		 // CrossingNumber += fabs((asin(Bound(n1x*n2x+n1y*n2y+n1z*n2z,-1,1))+asin(Bound(n2x*n3x+n2y*n3y+n2z*n3z,-1,1))+asin(Bound(n3x*n4x+n3y*n4y+n3z*n4z,-1,1))+asin(Bound(n4x*n1x+n4y*n1y+n4z*n1z,-1,1)))/(2*M_PI));
		}
	      else if (n5x*r13x+n5y*r13y+n5z*r13z < -ZERO)
		{
		  Writhe += - (asin(Bound(n1x*n2x+n1y*n2y+n1z*n2z,-1,1))+asin(Bound(n2x*n3x+n2y*n3y+n2z*n3z,-1,1))+asin(Bound(n3x*n4x+n3y*n4y+n3z*n4z,-1,1))+asin(Bound(n4x*n1x+n4y*n1y+n4z*n1z,-1,1)))/(2*M_PI);
		//  CrossingNumber += fabs((asin(Bound(n1x*n2x+n1y*n2y+n1z*n2z,-1,1))+asin(Bound(n2x*n3x+n2y*n3y+n2z*n3z,-1,1))+asin(Bound(n3x*n4x+n3y*n4y+n3z*n4z,-1,1))+asin(Bound(n4x*n1x+n4y*n1y+n4z*n1z,-1,1)))/(2*M_PI));
		}
	    }
	}
return Writhe;
}    

double getx(long current_frame, long num_atoms, long atom_num1, double *alval)
{  return alval[current_frame*num_atoms*3 + atom_num1*3];
}

double gety(long current_frame, long num_atoms, long atom_num1, double *alval)
{  return alval[current_frame*num_atoms*3 + atom_num1*3+1];
}
double getz(long current_frame, long num_atoms, long atom_num1, double *alval)
{  return alval[current_frame*num_atoms*3 + atom_num1*3+2];
}

void get_vector(long atom_num, long num_atoms, long current_frame, double *res, double *alval)
{
double dx, dy, dz;
if(atom_num==0){
  cout<<"Access to this element should not happen."<<endl;
  exit (EXIT_FAILURE);
}
dx=getx(current_frame , num_atoms, atom_num, alval)-getx(current_frame , num_atoms, atom_num-1 , alval);
dy=gety(current_frame , num_atoms, atom_num, alval)-gety(current_frame , num_atoms, atom_num-1 , alval);
dz=getz(current_frame , num_atoms, atom_num, alval)-getz(current_frame , num_atoms, atom_num-1 , alval);
res[0]=dx;
res[1]=dy;
res[2]=dz;
}


double norm_scalarp(long atom_num1, long atom_num2, long num_atoms, long current_frame, double *alval)
{
double v1[3], v2[3], norm1, norm2;
get_vector(atom_num1, num_atoms, current_frame, v1, alval);
get_vector(atom_num2, num_atoms, current_frame, v2, alval);

norm1=pitagora(v1[0],v1[1],v1[2]);
norm2=pitagora(v2[0],v2[1],v2[2]);

return  ( v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2] )/(norm1*norm2);

}


double scalar_prod(long atom_num1, long atom_num2, long num_atoms, long current_frame, double *alval)
{
double v1[3], v2[3];
get_vector(atom_num1, num_atoms, current_frame, v1, alval);
get_vector(atom_num2, num_atoms, current_frame, v2, alval);

return  ( v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2] );

}


double scalar_prod(double x, double y, double z, double x1, double y1, double z1, bool norm){
if (norm==true)
 return ( x*x1 +y*y1 + z*z1 )/(pitagora(x,y,z)*pitagora(x1,y1,z1))  ;
else
 return ( x*x1 +y*y1 + z*z1 );
}


double calc_bending_en(long atom_num1, long atom_num2, long num_atoms, long current_frame, double *alval)
{
double v=0, t=0;
ofstream outfile;

for(long i=atom_num1;i<atom_num2-1;++i)
{
 t=norm_scalarp(i, i+1, num_atoms, current_frame, alval);
 if (t>1.0){
  cout<<"error, scalar product is bigger than 1 for the atoms: "<<i<<" and "<<i+1<<endl;
  return 0;
 }
 t=  acos(norm_scalarp(i, i+1, num_atoms, current_frame, alval));
 v+=t*t;
}

 return v;
}


void frames_atoms_size(string filename, long &f_num, long &a_num)
{
  string headt1="Num_of_atoms", temp1;
  long num_frames=0, num_atoms=0;
  ifstream infile;
  infile.open(filename.c_str());
  string s;
  infile>>temp1; //first part will identify the file

  // old file format were starting with "Num_of_atoms"
  if(temp1.compare(headt1)==0){
  infile>>num_atoms; // this is the number
  cout<<"Num of atoms: "<<num_atoms<<endl;
  // here the two formats that we use differ, one doesn't have "Num_of_frames", the other has it
  infile>>s;
  infile>>num_frames;
  cout<<"Num of frames: "<<num_frames<<endl;
  f_num=num_frames;
  a_num=num_atoms;
  }

  // new file format start with an integer number that give the number of atoms
  else{
  a_num=atof(temp1.c_str());
  long cf=0;
  while(getline(infile,s))
   cf++;
   
  f_num=cf/a_num;
  }
  
  infile.close();  infile.clear();
}


void load_matrix(double *data, string namefile)
{
ifstream infile;
long num_atoms, num_frames;
vector <double> datav;
string s, head1="Num_of_atoms";
infile.open(namefile.c_str());
infile>>s; //in the file this is a string that declare the name of the variable, "num_atoms"
 if(s.compare(head1)==0){
  infile>>num_atoms;
  infile>>s; //in the file this is a string that declare the name of the variable, "num_frames"
  infile>>num_frames;
  
  // we reformat the input file
  long i=0;
  while(infile>>s){
    s.erase(remove(s.begin(), s.end(), '{'), s.end()); 
    s.erase(remove(s.begin(), s.end(), '}'), s.end());
    data[i]=(atof( s.c_str() ));
    i++;
  }
 }  
 else{

  long i=0;
  while(infile>>s && !infile.eof()){
   data[i]=atof(s.c_str());
   i++;
  }
 }

infile.close(); infile.clear();

}




double error_bar(double *data, long size)
{
Tstatistic statx;
for(long i=0;i<size;++i)
  statx.Add(data[i]);

long k=statx.Get_koptimal();


// cout<<"#k BinLength NbBins Average ErrorBars Correlation\n";	  
// cout<<k<<" "<<statx.GetBinLength(k)<<" "<<statx.GetNbBins(k)<<" "<<statx.GetAverage(k)<<" "<<statx.GetMeanSquare(k)<<" "<<statx.GetCorrelation(k)<<"\n";

 return statx.GetMeanSquare(k);

}



//calculate the distance between two beads at a certain frame
double distance(long current_frame, long num_atoms, long atom_num1, long atom_num2, double *alval)
{
  double x_a1, y_a1, z_a1, x_a2, y_a2, z_a2;
  x_a1=getx(current_frame, num_atoms, atom_num1, alval);
  y_a1=gety(current_frame, num_atoms, atom_num1, alval);
  z_a1=getz(current_frame, num_atoms, atom_num1, alval);

  x_a2=getx(current_frame, num_atoms, atom_num2, alval);
  y_a2=gety(current_frame, num_atoms, atom_num2, alval);
  z_a2=getz(current_frame, num_atoms, atom_num2, alval);
  
  return pitagora (x_a1-x_a2, y_a1-y_a2, z_a1-z_a2);
}



void center_mass_xyz(long current_frame, long num_atoms, long start, long end, double *res, double *alval)
{
  double x_a1=0.0, y_a1=0.0, z_a1=0.0;

    for(long i=start; i<end;++i)
    {
    x_a1+=alval[current_frame*num_atoms*3 + i*3];
    y_a1+=alval[current_frame*num_atoms*3 +i*3 +1];
    z_a1+=alval[current_frame*num_atoms*3 +i*3 +2];
    }
    res[0]=x_a1/float(end-start);
    res[1]=y_a1/float(end-start);
    res[2]=z_a1/float(end-start);
}




// calculate the average distance from the bead index_atom_reference to all the others (max: last_index_longerested)
// averaging in all the frames
// the routine allow also to calculate the average distance for a section
// of the chain, starting from "index_atom_reference" and ending with the atom "last_index_longerested"
double ave_frame_distace(long num_frames, long num_atoms, long index_atom_reference, long first_index_interested, long last_index_interested, double *alval, double *results)
{

double xref, yref, zref, xn, yn, zn;
long length_section=fabs(first_index_interested-last_index_interested);
double temporaney[length_section];
for(long g=0; g<length_section;++g)
  temporaney[g]=0.0;

    for(long j=0; j<num_frames; ++j)
    {
      xref=alval[j*num_atoms*3 + index_atom_reference*3 ];
      yref=alval[j*num_atoms*3 + index_atom_reference*3+1];
      zref=alval[j*num_atoms*3 + index_atom_reference*3+2];
      
      for(long i=first_index_interested; i<last_index_interested; ++i)
      {
      xn=alval[j*num_atoms*3 + i*3 ];
      yn=alval[j*num_atoms*3 + i*3+1];
      zn=alval[j*num_atoms*3 + i*3+2];
      
      temporaney[i-first_index_interested]+=pitagora(xref-xn, yref-yn, zref-zn);
      }

    }
    
    for(long j=0; j<length_section;++j)
      results[j]= temporaney[j]/(float(num_frames));

}


// calculate the radius of gyration for a frame  formulas from: http://en.wikipedia.org/wiki/Radius_of_gyration
// the value is root_squared
double radius_gyration(long current_frame, long num_atoms, double *alval){

double cm_pos[3];
center_mass_xyz(current_frame, num_atoms, 0, num_atoms, cm_pos, alval);

double sum=0.0, temp, xn, yn, zn;
  for(long j=0; j<num_atoms; ++j)
  {
  xn=alval[current_frame*num_atoms*3 + j*3 ];
  yn=alval[current_frame*num_atoms*3 + j*3+1];
  zn=alval[current_frame*num_atoms*3 + j*3+2];
  temp=pitagora(xn-cm_pos[0],yn-cm_pos[1],zn-cm_pos[2]);
  sum+= temp*temp;
  }

  sum/=num_atoms;

  return sqrt(sum);
}

// calculate the radius of gyration for a molecule in a frame
// here we use "start" and "end" to compute the Rg of chunks of a molecule or different molecules
// the value is root_squared
double radius_gyration(long current_frame, long num_atoms, double *alval, long start, long end){

double cm_pos[3];
center_mass_xyz(current_frame, num_atoms, start, end, cm_pos, alval);
double sum=0.0, temp, xn, yn, zn;

  for(long j=start; j<end; ++j)
  {
  xn=alval[current_frame*num_atoms*3 + j*3 ];
  yn=alval[current_frame*num_atoms*3 + j*3+1];
  zn=alval[current_frame*num_atoms*3 + j*3+2];
  temp=pitagora(xn-cm_pos[0],yn-cm_pos[1],zn-cm_pos[2]);
  sum+= temp*temp;
  }

  sum/=float(end-start);

  return sqrt(sum);
}


double end_to_end_distance(vector <double> &x, vector <double> &y, vector <double> &z, int window, bool circular)
{


int N=x.size();
if(window<=0 && window>N){
cout<<"Error, windows should be at least bigger than zero and small than the size of the polymer."<<endl;
return sqrt(-1);  //not a number
}


vector <double> temp;
double dval;
if(circular){
  for(int i=0; i<N; ++i){
   dval=pitagora(x[i]-x[(i+window)%N], y[i]-y[(i+window)%N], z[i]-z[(i+window)%N]);
   dval*=dval;
   temp.push_back(dval);
  }
}
else{
  for(int i=0; i<N-window; ++i){
   dval=pitagora(x[i]-x[(i+window)%N], y[i]-y[(i+window)%N], z[i]-z[(i+window)%N]);
   dval*=dval;
   temp.push_back(dval);
  }
}
double avev;
average(temp, avev);

return avev;
}
 



// remap a part of the chain to remove the boundary conditions starting from the first value
void remap_part(double *chainx, double *chainy, double *chainz, long natoms, double sizebox_x, double sizebox_y, double sizebox_z)
{
  double sizemin=min(sizebox_x,sizebox_y);
  sizemin=min(sizemin, sizebox_z);
  long cs=1, nx=0,ny=0,nz=0, min=-1, max=1;
  while (cs< natoms){
     double best=sizebox_x+sizebox_y+sizebox_z;
     long x_best=0,y_best=0,z_best=0;
     for(long x=min; x<=max;++x){
       for(long y=min; y<=max;++y){
         for(long z=min; z<=max;++z){
	   double v= pitagora ( chainx[cs-1]- (chainx[cs]+(x+nx)*sizebox_x),
			        chainy[cs-1]- (chainy[cs]+(y+ny)*sizebox_y),
			        chainz[cs-1]- (chainz[cs]+(z+nz)*sizebox_z) ) ; 
	   if(fabs(v)<best)
	       {
	       best=fabs(v);
	       x_best=x;
	       y_best=y;
	       z_best=z;
	       }
         }
       }
     }
     chainx[cs]=chainx[cs]+(x_best+nx)*sizebox_x;
     chainy[cs]=chainy[cs]+(y_best+ny)*sizebox_y;
     chainz[cs]=chainz[cs]+(z_best+nz)*sizebox_z;
     nx=x_best+nx;
     ny=y_best+ny;
     nz=z_best+nz;
     if(best<sizemin/2,++cs);// in case the distance is several box size
  }

}

void remap_noboundary(double *data, long num_frames, long num_atoms){


// size of the box,  bl "bond length", a temporaney variable
double sizebox_x, sizebox_y, sizebox_z;

string s;
cout<<"Is it a cubic box? (press \"y\" for yes, any other for no )"<<endl;
getline(cin, s);

if(s=="y"){
cout<<"Size of the box?"<<endl;
cin>>sizebox_x;
sizebox_y=sizebox_x;
sizebox_z=sizebox_x;
}

else{
cout<<"Size of the x coordinate of the box?"<<endl;
cin>>sizebox_x;
cout<<"Size of the y coordinate of the box?"<<endl;
cin>>sizebox_y;
cout<<"Size of the z coordinate of the box?"<<endl;
cin>>sizebox_z;
}

remap_noboundary(data, num_frames, num_atoms, sizebox_x, sizebox_y, sizebox_z);
   
}


void remap_noboundary(double *data, long num_frames, long num_atoms, double sizebox_x, double sizebox_y, double sizebox_z){


// size of the box,  bl "bond length", a temporaney variable
//sizebox_x, sizebox_y, sizebox_z
double bl;


vector <double> xcoord;
vector <double> ycoord;
vector <double> zcoord;

for(long j=0;j<num_frames*num_atoms*3;++j){
 if(j%3==0)
   xcoord.push_back(data[j]);
 if(j%3==1)
   ycoord.push_back(data[j]);
 if(j%3==2)
   zcoord.push_back(data[j]);
 }


vector <double> xnew;
vector <double> ynew;
vector <double> znew;



cout<<"Remapping without periodic boundary conditions..."<<endl;
for(long numf=0; numf<num_frames;++numf){
  xnew.push_back(xcoord[numf*num_atoms]);
  ynew.push_back(ycoord[numf*num_atoms]);
  znew.push_back(zcoord[numf*num_atoms]);
  bl=0;
  long cs=1+numf*num_atoms, nx=0,ny=0,nz=0, min=-1, max=1;
  while (cs< (numf+1)*num_atoms  ){
     double best=1000;
     long x_best=0,y_best=0,z_best=0;
     for(long x=min; x<=max;++x){
       for(long y=min; y<=max;++y){
         for(long z=min; z<=max;++z){
	   double v= pitagora ( xnew.back() - (xcoord[cs]+(x+nx)*sizebox_x),
			        ynew.back()- (ycoord[cs]+(y+ny)*sizebox_y),
			        znew.back()- (zcoord[cs]+(z+nz)*sizebox_z) ) ; 
	   if(fabs(v-bl)<best)
	     {
	       best=fabs(v-bl);
	       x_best=x;
	       y_best=y;
	       z_best=z;
	     }
         }
       }
     }
     xnew.push_back(xcoord[cs]+(x_best+nx)*sizebox_x);
     ynew.push_back(ycoord[cs]+(y_best+ny)*sizebox_y);
     znew.push_back(zcoord[cs]+(z_best+nz)*sizebox_z);
     nx=x_best+nx;
     ny=y_best+ny;
     nz=z_best+nz;
    ++cs;
  }

}

  
  for(long i=0;i<num_frames*num_atoms;++i){
    data[i*3]=xnew[i];
    data[i*3+1]=ynew[i];
    data[i*3+2]=znew[i];
  }
  
  xnew.clear();
  ynew.clear();
  znew.clear();
  xcoord.clear();
  ycoord.clear();
  zcoord.clear();
   
}




void observable_tstat(vector <vector <double> > &all_data, vector <string> nameobs, long current_f){
Tstatistic temp;

vector <Tstatistic> all_obs;
ifstream infile;
infile.open("all_stats.txt");
if(infile){  //if the file exist we need to read from the file the status of each observable
 for (long ll=0;ll<all_data.size();++ll){  // reading from the file all the previous states for all the observables
  
  infile>>temp;
  all_obs.push_back(temp);  // all_obs contain a certain number of Tstatistic initialized from file
 }
   
  if(all_data.size()!=all_obs.size())
  { cerr<<"Number of data and number of observable are inconsistent. ERROR."<<endl; exit(1); }  

  // now we need to add the current chunk of data
  for(long i=0;i<all_data.size();++i)
    for(long j=0;j<all_data[i].size();++j)
      all_obs[i].Add(all_data[i][j]);

}

else{  // if there IS NOT a previous state we need to create the Tstatistc from the data, and append it to all_obs
  for(long i=0;i<all_data.size();++i)
  { 

    Tstatistic temp2;
    for(long j=0;j<all_data[i].size();++j)
    {  
      temp2.Add(all_data[i][j]);
    }    
  all_obs.push_back(temp2);
    
  }
}
infile.close(); infile.clear();


// now we need to write everything in a file
ofstream outfile;
outfile.open("all_stats.txt");
for(long i=0;i<all_obs.size();++i)
{
// cout<<"writing stats of data: "<<i<<endl;
 outfile<<all_obs[i]<<" ";
}
outfile.close(); outfile.clear();


//now we are going to write the current average (from the beginning to this chunk) on a file
// here we deal with the header, if the file does NOT exist we create it
infile.open("all_obs_averages.txt", ios::in);
if(!infile)
  {outfile.open("all_obs_averages.txt", ios::out);
  if(outfile)
   {
     outfile<<"Frame# ";
     for(long n=0;n<nameobs.size();++n)
     {
     outfile<<nameobs[n]<<".optimalBin "<<nameobs[n]<<".Average "<<nameobs[n]<<".Errorbar "<<nameobs[n]<<".Autocorr ";     
     }
     outfile<<endl;
    outfile.close(); outfile.clear();
   }
  }
else
  {infile.close(); infile.clear();}

// we need that all the observables have the same range for the binning of the data,
// this mean that all the vectors have to be of the same length, otherwise we are $%&@ed
long firstbin=all_obs[0].GetNbStatistics();
for(long h=0; h< all_obs.size();++h)
  if(all_obs[h].GetNbStatistics()!=firstbin)
  {
   cerr<<"An observable is not coerent in length."<<endl;
   exit(1);
  }


outfile.open("all_obs_averages.txt", ios::app);
outfile<<current_f<<" ";
for(long o=0; o<all_obs.size();++o)
{
long optb=all_obs[o].Get_koptimal();

outfile<<optb<<" "<<all_obs[o].GetAverage(0)<<" "<<all_obs[o].GetMeanSquare(optb)<<" "<<all_obs[o].GetCorrelation(optb)<<" ";
}
outfile<<endl;
outfile.close(); outfile.clear();


// in this file we write all the values for all the bin considering the entire trajectory
outfile.open("all_bins_all_obs.txt", ios::out);


for(long n=0;n<nameobs.size();++n){
for(long h=0;h<all_obs[0].GetNbStatistics();++h){
  outfile<<nameobs[n]<<".Frame# "<<nameobs[n]<<".Bin "<<nameobs[n]<<".Average "<<nameobs[n]<<".Errorbar "<<nameobs[n]<<".Autocorr ";
  }
}
outfile<<endl;



for(long l=0;l<all_obs.size();++l){
for(long h=0;h<all_obs[0].GetNbStatistics();++h){
   outfile<<current_f<<" "<<h<<" "<<all_obs[l].GetAverage(h)<<" "<<all_obs[l].GetMeanSquare(h)<<" "<<all_obs[l].GetCorrelation(h)<<" ";
  }
}

outfile.close(); outfile.clear();


}





void write_CM(long start, long end, long num_frames, long num_atoms, double epsilon, double *data){
//contact map calculation
 ofstream outfile;
 ifstream infile;
// initialize a single block of memory to store the contact map matrix
long sizem=end-start;

vector <vector <double> > map;
for(long j=0;j<sizem;++j){
    vector <double> temp;
    for(long k=0;k<sizem; ++k){
        temp.push_back(0);
    }
    map.push_back(temp);
    temp.clear();
}

timeval t1, t2;
double elapsedTime;

//start timer
//gettimeofday(&t1, NULL);
//cout<<"Calculating the contact map..."<<endl;
 
  string eps_str = num_to_str(epsilon);
  string namef="contactmap_eps_"+eps_str+".txt";
  infile.open(namef.c_str(), ios::in);
    if(infile){
    cout<<"previous contact map exist"<<endl; 
      for(long j=0;j<sizem;++j)
       for(long i=0;i<sizem;++i)
       {
        double mv;
        infile>>mv;
        map[j][i]=mv;
       }

     }
  infile.close(); infile.clear(); 
     #pragma omp parallel for shared(map)
     for(long i=start; i<end;++i){
      for(long j=0;j<num_frames;++j){
      long k=start;
       while(k<end){
        double dis=distance( j, num_atoms, i, k, data)-epsilon;
        if(dis<=0 ){
          map[i-start][k-start]+=1.0;
          k++;     
        }
        else{ k+= ceil(dis*0.7); }
       }
      }
     }
 
  outfile.open(namef.c_str());
  outfile.precision(12);
  for(long i=0; i<sizem;++i)
  {
    for(long k=0; k<sizem;++k)
     outfile<<map[i][k]<<" ";
     outfile<<endl;
  }
  outfile.close();  outfile.clear();

//// stop timer
//gettimeofday(&t2, NULL);

// compute and print the elapsed time in millisec
//elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;
//cout<<"Time to calculate the map in millisec: "<<elapsedTime<<endl;


vector <double> avevals;
vector <double> stdvals_mean;
vector <double> medians;
for(long j=0;j<sizem;++j){
 vector <double> temp;
 long i=0, k=j;
 while(i<(sizem-j) && k<sizem){
  temp.push_back(map[i][k]/map[0][0]);  // dividing by the first position give us a normalized contact map
  i+=1;
  k+=1;
 }
  Tstatistic a;
 for(int vali=0;vali<temp.size();++vali)
 a.Add(temp[vali]);
 
 avevals.push_back(a.GetAverage());
 stdvals_mean.push_back(a.GetMeanSquare(a.Get_koptimal()));
 temp.clear();
 }

for(long i=0;i<map.size();++i)
    map[i].clear();
map.clear();

namef="1D_contactmap_eps"+eps_str+".txt";
outfile.open(namef.c_str());
outfile.precision(12);
outfile<<"Distance Average Std_dev_av Median"<<endl;
for(long i=0;i<avevals.size();++i)
 outfile<<i<<" "<<avevals[i]<<" "<<stdvals_mean[i]<<endl;
outfile.close(); outfile.clear();

avevals.clear();
stdvals_mean.clear();


}



//////////////////////////////////////////////////////////////////
void calc_dihedral(double *data, long num_atoms, long num_frames, long num_real_beads, long start, long stop)
{
  double twist=0;
  double r1x,r1y,r1z,r2x,r2y,r2z;
  double e1x,e2x,e1y,e2y,e1z,e2z;
  double nx,ny,nz;
  double nn,nr1,nr2,ne1,ne2;
  double see,se1,se2,s,s1,s2;
  double alpha, a1,a2;
  double Pi=M_PI;

ofstream outfile;
outfile.open("dihedral.dat");

for(int f=0;f<num_frames;++f){
//  cout<<"frame n:"<<f<<" on : "<<num_frames<<endl;
  for(int a=start; a<stop-2;++a){
  r1x=getx(f, num_atoms, a, data);
  r1y=gety(f, num_atoms, a, data);
  r1z=getz(f, num_atoms, a, data);
  r1x-=getx(f, num_atoms, a+1, data);
  r1y-=gety(f, num_atoms, a+1, data);
  r1z-=getz(f, num_atoms, a+1, data);
  nr1=pitagora(r1x,r1y,r1z);
  r1x/=nr1; r1y/=nr1; r1z/=nr1;

  r2x=getx(f, num_atoms, a+1, data);
  r2y=gety(f, num_atoms, a+1, data);
  r2z=getz(f, num_atoms, a+1, data);
  r2x-=getx(f, num_atoms, a+2, data);
  r2y-=gety(f, num_atoms, a+2, data);
  r2z-=getz(f, num_atoms, a+2, data);
  nr2=pitagora(r2x,r2y,r2z);
  r2x/=nr2; r2y/=nr2; r2z/=nr2;

  e1x=getx(f, num_atoms, a+2*num_real_beads, data);
  e1y=gety(f, num_atoms, a+2*num_real_beads, data);
  e1z=getz(f, num_atoms, a+2*num_real_beads, data);
  e1x-=getx(f, num_atoms, a+num_real_beads, data);
  e1y-=gety(f, num_atoms, a+num_real_beads, data);
  e1z-=getz(f, num_atoms, a+num_real_beads, data);
  ne1=pitagora(e1x,e1y,e1z);
  e1x/=ne1; e1y/=ne1; e1z/=ne1;

  e2x=getx(f, num_atoms, a+2*num_real_beads+1, data);
  e2y=gety(f, num_atoms, a+2*num_real_beads+1, data);
  e2z=getz(f, num_atoms, a+2*num_real_beads+1, data);
  e2x-=getx(f, num_atoms, a+num_real_beads+1, data);
  e2y-=gety(f, num_atoms, a+num_real_beads+1, data);
  e2z-=getz(f, num_atoms, a+num_real_beads+1, data);
  ne2=pitagora(e2x,e2y,e2z);
  e2x/=ne2; e2y/=ne2; e2z/=ne2;


	//normal vector to the plane containing the two segments
	nx=r1y*r2z-r1z*r2y;
	ny=r1z*r2x-r1x*r2z;
	nz=r1x*r2y-r1y*r2x;
	nn=pitagora(nx,ny,nz);
        //cout<<ne1<<" "<<ne2<<" nn "<<nn<<endl;

	if(nn<0.0000000001)
		{//vectors parallel, see=scalar product, s=orientation
		see=e1x*e2x+e1y*e2y+e1z*e2z;
		s=threeprod(e1x,e1y,e1z,e2x,e2y,e2z,r1x,r1y,r1z);
 //   if(fabs(s)<1e-6)cerr<<"s="<<s<<" i="<<i<<endl; 
 //   if(fabs(see)>1-1e-6)cerr<<"see="<<see<<" i="<<i<<endl; 

		alpha=acos(see)*s/fabs(s);//alpha=fabs(fmod(acos(see),Pi))*s/fabs(s);
                cout<<see<<" normal vector was 0 "<<s<<" "<<alpha<<endl;
		}
	else
		{
		nx=nx/nn;ny=ny/nn;nz=nz/nn;
		se1=e1x*nx+e1y*ny+e1z*nz;
		se2=e2x*nx+e2y*ny+e2z*nz;
		s1=threeprod(e1x,e1y,e1z,nx,ny,nz,r1x,r1y,r1z);
    //special case, e1 is parallel to n, then the angle is 0 or Pi
    if(fabs(s1)<1e-10)	{a1=0; cout<<"normal vector parallel to a1-0"<<endl;
			if(e1x*nx+e1y*ny+e1z*nz<0){a1=Pi; //cout<<"normal vector parallel to a1-Pi"<<endl;
      }
		}
    //'general' case:
    else a1=acos(se1)*s1/fabs(s1);//a1=fabs(fmod(acos(se1),Pi))*s1/fabs(s1);

		s2=threeprod(nx,ny,nz,e2x,e2y,e2z,r1x,r1y,r1z);//r2
 		if(fabs(s2)<1e-10) {a2=0;
			if(e2x*nx+e2y*ny+e2z*nz<0){a2=Pi; //cout<<"normal vector parallel to a2-Pi"<<endl;
       }
      }
		else a2=acos(se2)*s2/fabs(s2);//a2=fabs(fmod(acos(se2),Pi))*s2/fabs(s2);

//    if(fabs(s1)<1e-10)cerr<<"s1="<<s1<<" i="<<i<<" a1="<<a1<<endl; 
// 		if(fabs(se1)>1-1e-10)cerr<<"se1="<<se1<<" i="<<i<<" a1="<<a1<<endl; 
//		if(fabs(s2)<1e-10)cerr<<"s2="<<s2<<" i="<<i<<" a2="<<a2<<endl; 
// 		if(fabs(se2)>1-1e-10)cerr<<"se2="<<se2<<" i="<<i<<" a2="<<a2<<endl; 
		alpha=a1+a2;
//		cout<<se1<<" "<<se2<<" nn was not 0 "<<s1<<" "<<s2<<" "<<alpha<<endl;
		}
	// rescale the angle so that it is between -Pi and Pi
	if(alpha>Pi){alpha=-2*Pi+alpha;}
	if(alpha<-Pi){alpha=2*Pi+alpha;}

	twist=twist+alpha;
	// cout<<a1<<" a12 " <<a2<<" alpha " <<alpha<<" twist "<<twist/2/Pi<<endl;
  outfile<<alpha<<" ";
	}
  outfile<<endl; 

//return twist/(2.0*Pi);
  }
  outfile.close(); 
}


// check if two beads A B that belong to a chain, in a periodic box, get in contact with the periodic copies B A 
//chain x,y,z contain the coordinates of the (real) beads of the chain for all the frames
//pos1 and pos2 are the index of the special sites
//natoms is the size of the ring
void check_periodic_contact(double *chainx, double *chainy, double *chainz, long natoms, long num_frames, double sizebox_x, double sizebox_y, double sizebox_z, int pos1, int pos2, long &real_contact_count, long &periodic_contact_count, double epsilon)
{
  for(int nf=0;nf<num_frames;++nf){
    //the following value is the first coordinate of the first bead in the frame nf
    int fc=nf*natoms;
    double pos1_x=chainx[fc+pos1], pos2_x=chainx[fc+pos2],pos1_y=chainy[fc+pos1], pos2_y=chainy[fc+pos2],pos1_z=chainz[fc+pos1], pos2_z=chainz[fc+pos2];
    //measure the distance between the two special bead: pos1 and pos2
    double maximum_dist=max(pos1_x-pos2_x,pos1_y-pos2_y);
    maximum_dist=max(maximum_dist,pos1_z-pos2_z);
    //get the minimum size of the box (even if it's usually a cube)
    double sizemin=min(sizebox_x,sizebox_y);
    sizemin=min(sizemin, sizebox_z);
    int max_periodic=ceil(maximum_dist/sizemin);
    for(long x=-max_periodic; x<=max_periodic;++x){
      int control_rcc=0, control_pcc=0;
      for(long y=-max_periodic; y<=max_periodic;++y){
        for(long z=-max_periodic; z<=max_periodic;++z){

          // create the coordinate for the periodic copy of pos2
          double pos2_x_pc=pos2_x + x*sizebox_x , pos2_y_pc=pos2_y + y*sizebox_y , pos2_z_pc=pos2_z + z*sizebox_z ;

          //measure the distance between the pos1 and the copy of pos2
          double copy_dist=pitagora(pos1_x-pos2_x_pc,pos1_y-pos2_y_pc,pos1_z-pos2_z_pc);
          if(copy_dist<=epsilon){
            if(x==0 && y==0 && z==0){
             control_rcc++;
             real_contact_count++;
            }
            else{
             periodic_contact_count++;
             control_pcc++;
             }
          }
        }
      }
    
    }
  }
}


/* Generate the distance matrix for a given chain x,y,z.  Putting in every cell between bead_i and bead_j
the value of their distance */
vector<double> get_distance_matrix(vector<double> & x,vector<double> & y,vector<double> & z)
{
  int N=x.size();
  vector<double> matrix(N*N);
  #pragma omp parallel for shared(matrix)
  for(int m=0;m<N;m++)
    for(int n=0;n<N;n++)
      {
	matrix[m*N+n]=sqrt((x[n]-x[m])*(x[n]-x[m])+(y[n]-y[m])*(y[n]-y[m])+(z[n]-z[m])*(z[n]-z[m]));
	matrix[n*N+m]=matrix[m*N+n];
      }
  return matrix;
}

//////////////////////////////////////////////////







//////////////////////////////////////////////////




//////////////////////////////////////////////////////
//full circle obtained with pos1=0, pos2=N
double get_writhe(int pos1,int pos2,vector<double> &x,vector<double> &y,vector<double> &z)//if pos1>pos2 => take subchain from pos1 to 0 and from 0 to pos2
{
  double ZERO=1e-10;
  int N=x.size();
  if(pos1<0||pos2<0)
    {
      pos1+=N;
      pos2+=N;
    }
  if(pos1>pos2)
    pos2=pos2+N;
  
  double Writhe=0;

  #pragma omp parallel for reduction(+:Writhe)
  for(int i1=pos1+1;i1<pos2;i1++)
  {
  int j1,j2;
    double norme;
    double r13x,r13y,r13z;
    double r14x,r14y,r14z;
    double r24x,r24y,r24z;
    double r23x,r23y,r23z;
    double r34x,r34y,r34z;
    double r12x,r12y,r12z;
    double n1x,n1y,n1z;   
    double n2x,n2y,n2z;   
    double n3x,n3y,n3z;   
    double n4x,n4y,n4z;   
    double n5x,n5y,n5z;  

  
      j1=(i1+1)%N;
      for(int i2=pos1;i2<i1;i2++)
	{
	  j2=(i2+1)%N;
	  r13x=x[i2%N]  -x[i1%N];
	  r13y=y[i2%N]  -y[i1%N];
	  r13z=z[i2%N]  -z[i1%N];
	  r14x=x[j2]-x[i1%N];
	  r14y=y[j2]-y[i1%N];
	  r14z=z[j2]-z[i1%N];
	  r24x=x[j2]-x[j1];
	  r24y=y[j2]-y[j1];
	  r24z=z[j2]-z[j1];
	  r23x=x[i2%N]  -x[j1];
	  r23y=y[i2%N]  -y[j1];
	  r23z=z[i2%N]  -z[j1];
	  r34x=x[j2]-x[i2%N];
	  r34y=y[j2]-y[i2%N];
	  r34z=z[j2]-z[i2%N];
	  r12x=x[j1]-x[i1%N];
	  r12y=y[j1]-y[i1%N];
	  r12z=z[j1]-z[i1%N];
	  n1x=r13y*r14z-r13z*r14y; n1y=r13z*r14x-r13x*r14z; n1z=r13x*r14y-r13y*r14x;
	  norme=sqrt(n1x*n1x+n1y*n1y+n1z*n1z);
	  if(norme>ZERO){n1x=n1x/norme; n1y=n1y/norme; n1z=n1z/norme;}
		    
	      n2x=r14y*r24z-r14z*r24y; n2y=r14z*r24x-r14x*r24z; n2z=r14x*r24y-r14y*r24x;
	      norme=sqrt(n2x*n2x+n2y*n2y+n2z*n2z);
	      if(norme>ZERO){n2x=n2x/norme; n2y=n2y/norme; n2z=n2z/norme;}
		  
	      n3x=r24y*r23z-r24z*r23y; n3y=r24z*r23x-r24x*r23z; n3z=r24x*r23y-r24y*r23x;
	      norme=sqrt(n3x*n3x+n3y*n3y+n3z*n3z);
	      if(norme>ZERO){n3x=n3x/norme; n3y=n3y/norme; n3z=n3z/norme;}
		  
	      n4x=r23y*r13z-r23z*r13y; n4y=r23z*r13x-r23x*r13z; n4z=r23x*r13y-r23y*r13x;
	      norme=sqrt(n4x*n4x+n4y*n4y+n4z*n4z);
	      if(norme>ZERO){n4x=n4x/norme; n4y=n4y/norme; n4z=n4z/norme;}
		  
	      n5x=r34y*r12z-r34z*r12y; n5y=r34z*r12x-r34x*r12z; n5z=r34x*r12y-r34y*r12x;
	      if(n5x*r13x+n5y*r13y+n5z*r13z > ZERO)
		{
		    Writhe += -(asin(Bound(n1x*n2x+n1y*n2y+n1z*n2z,-1,1))+asin(Bound(n2x*n3x+n2y*n3y+n2z*n3z,-1,1))+asin(Bound(n3x*n4x+n3y*n4y+n3z*n4z,-1,1))+asin(Bound(n4x*n1x+n4y*n1y+n4z*n1z,-1,1)))/(2*M_PI);
		}
	      else if (n5x*r13x+n5y*r13y+n5z*r13z < -ZERO)
		{
		    Writhe += (asin(Bound(n1x*n2x+n1y*n2y+n1z*n2z,-1,1))+asin(Bound(n2x*n3x+n2y*n3y+n2z*n3z,-1,1))+asin(Bound(n3x*n4x+n3y*n4y+n3z*n4z,-1,1))+asin(Bound(n4x*n1x+n4y*n1y+n4z*n1z,-1,1)))/(2*M_PI);
		}		    
	    }
	}
  return -Writhe;
}

/*
After completion, writhe[pos][w] contains the writhe for postion pos=0,..,Nbeads-1 and window of size 2*half_windows_sizes[w]

*/
void get_writhe_new(vector<int> half_windows_sizes,const vector<double> &x,const vector<double> &y,const vector<double> &z,vector<vector<double> > & writhe)//if pos1>pos2 => take subchain from pos1 to 0 and from 0 to pos2
{
  double ZERO=1e-10;
  int N=x.size();
  int maxwindow=2*(*max_element(half_windows_sizes.begin(),half_windows_sizes.end()));
//  cerr<<"maxwindow="<<maxwindow<<endl;
  writhe.resize(N);
  for(int i=0;i<writhe.size();i++){
    writhe[i].resize(half_windows_sizes.size());
  }


  /////precompute writhe_local
//  cerr<<"precomputing"<<endl;
  double *writhe_local=new double[N*N];
  for(int i=0;i<N*N;++i)
   writhe_local[i]=0;
  
  int pos1,pos2;
  for(int i1=0;i1<N;i1++)
  {
      int j1,j2;
      double norme;
      double r13x,r13y,r13z;
      double r14x,r14y,r14z;
      double r24x,r24y,r24z;
      double r23x,r23y,r23z;
      double r34x,r34y,r34z;
      double r12x,r12y,r12z;
      double n1x,n1y,n1z;   
      double n2x,n2y,n2z;   
      double n3x,n3y,n3z;   
      double n4x,n4y,n4z;   
      double n5x,n5y,n5z;  
      j1=(i1+1)%N;
      for(int i2=i1-maxwindow+N;i2<i1+N;i2++)
//      for(int i2=0;i2<N;i2++)
	    {
	      j2=(i2+1)%N;
	      r13x=x[i2%N]  -x[i1%N];
	      r13y=y[i2%N]  -y[i1%N];
	      r13z=z[i2%N]  -z[i1%N];
	      r14x=x[j2]-x[i1%N];
	      r14y=y[j2]-y[i1%N];
	      r14z=z[j2]-z[i1%N];
	      r24x=x[j2]-x[j1];
	      r24y=y[j2]-y[j1];
	      r24z=z[j2]-z[j1];
	      r23x=x[i2%N]  -x[j1];
	      r23y=y[i2%N]  -y[j1];
	      r23z=z[i2%N]  -z[j1];
	      r34x=x[j2]-x[i2%N];
	      r34y=y[j2]-y[i2%N];
	      r34z=z[j2]-z[i2%N];
	      r12x=x[j1]-x[i1%N];
	      r12y=y[j1]-y[i1%N];
	      r12z=z[j1]-z[i1%N];
	      n1x=r13y*r14z-r13z*r14y; n1y=r13z*r14x-r13x*r14z; n1z=r13x*r14y-r13y*r14x;
	      norme=sqrt(n1x*n1x+n1y*n1y+n1z*n1z);
	      if(norme>ZERO){n1x=n1x/norme; n1y=n1y/norme; n1z=n1z/norme;}
	      
	      n2x=r14y*r24z-r14z*r24y; n2y=r14z*r24x-r14x*r24z; n2z=r14x*r24y-r14y*r24x;
	      norme=sqrt(n2x*n2x+n2y*n2y+n2z*n2z);
	      if(norme>ZERO){n2x=n2x/norme; n2y=n2y/norme; n2z=n2z/norme;}
	      
	      n3x=r24y*r23z-r24z*r23y; n3y=r24z*r23x-r24x*r23z; n3z=r24x*r23y-r24y*r23x;
	      norme=sqrt(n3x*n3x+n3y*n3y+n3z*n3z);
	      if(norme>ZERO){n3x=n3x/norme; n3y=n3y/norme; n3z=n3z/norme;}
	      
	      n4x=r23y*r13z-r23z*r13y; n4y=r23z*r13x-r23x*r13z; n4z=r23x*r13y-r23y*r13x;
	      norme=sqrt(n4x*n4x+n4y*n4y+n4z*n4z);
	      if(norme>ZERO){n4x=n4x/norme; n4y=n4y/norme; n4z=n4z/norme;}
	      
	      n5x=r34y*r12z-r34z*r12y; n5y=r34z*r12x-r34x*r12z; n5z=r34x*r12y-r34y*r12x;
	      if(n5x*r13x+n5y*r13y+n5z*r13z > ZERO)
	      {
	          writhe_local[i1*N+i2%N] = (asin(Bound(n1x*n2x+n1y*n2y+n1z*n2z,-1,1))+asin(Bound(n2x*n3x+n2y*n3y+n2z*n3z,-1,1))+asin(Bound(n3x*n4x+n3y*n4y+n3z*n4z,-1,1))+asin(Bound(n4x*n1x+n4y*n1y+n4z*n1z,-1,1)))/(2*M_PI);
	      }
	      else if (n5x*r13x+n5y*r13y+n5z*r13z < -ZERO)
	      {
	          writhe_local[i1*N+i2%N] = -(asin(Bound(n1x*n2x+n1y*n2y+n1z*n2z,-1,1))+asin(Bound(n2x*n3x+n2y*n3y+n2z*n3z,-1,1))+asin(Bound(n3x*n4x+n3y*n4y+n3z*n4z,-1,1))+asin(Bound(n4x*n1x+n4y*n1y+n4z*n1z,-1,1)))/(2*M_PI);
	      }
	    }
  }

  //evaluate writhe
//  cerr<<"evaluating writhe"<<endl;

  for(int w=0; w<half_windows_sizes.size();++w)
  {
   for(int pos=0; pos<N;++pos)
	 {
	  int half_wind=half_windows_sizes[w];
	  pos1=pos-half_wind;
	  pos2=pos+half_wind;
	  
	  if(pos1<0||pos2<0)
	  {
	      pos1+=N;
	      pos2+=N;
	  }
	  if(pos1>pos2)
	    pos2=pos2+N;
	  //compute only first
	  if(pos==0)
	  {
	    double Writhe=0;
	    for(int i1=pos1+1;i1<pos2;i1++)
		  {
		    for(int i2=pos1;i2<i1;i2++)
		    {
		      Writhe +=writhe_local[(i1%N)*N+i2%N];
		    }
		  }
	    writhe[pos][w]=Writhe;
	  }
	  else
	  {
	    double Writhe=writhe[pos-1][w];
	      
	    int i2=pos1-1+N;
	    for(int i1=pos1;i1<pos2-1;i1++)
		   Writhe -=writhe_local[(i1%N)*N+i2%N];
	      
	    int i1=pos2-1+N;
	    for(int i2=pos1;i2<pos2-1;i2++)
		   Writhe +=writhe_local[(i1%N)*N+i2%N];
	    writhe[pos][w]=Writhe;
	  }
	 }
  }
  delete[] writhe_local;
}








//full circle obtained with pos1=0, pos2=N+1
double get_twist(int pos1, int pos2, vector <double> &chain_x, vector <double> &chain_y, vector <double> &chain_z, vector <double> &ph_x, vector <double> &ph_y, vector <double> &ph_z, vector <double> &th_x, vector <double> &th_y, vector <double> &th_z)
{

  int NbSegments=chain_x.size();
  int N=chain_x.size();
  if(pos1<0||pos2<0)
    {
      pos1+=N;
      pos2+=N;
    }
  if(pos1>pos2)
    pos2=pos2+N;


  double twist=0;
  double r1x,r1y,r1z,r2x,r2y,r2z, rp1x, rp1y, rp1z, rp2x, rp2y, rp2z;
  double e1x,e2x,e1y,e2y,e1z,e2z;
  double nx,ny,nz;
  double prod;
  double nn,nr1,nr2,ne1,ne2, np1, np2;
  double see,se1,se2,s,s1,s2;
  double alpha, a1,a2;
  double Pi=M_PI;
  
   for(long i=pos1;i<pos2-1;i++)
	{
	//first segment
	r1x=chain_x[(i+1)%NbSegments]-chain_x[i%NbSegments];
	r1y=chain_y[(i+1)%NbSegments]-chain_y[i%NbSegments];
	r1z=chain_z[(i+1)%NbSegments]-chain_z[i%NbSegments];
	//normalize
	nr1=pitagora(r1x,r1y,r1z);
//  if(nr1<1e-6)cerr<<"nr1="<<nr1<<" i="<<i<<endl; 
	r1x=r1x/nr1; r1y=r1y/nr1; r1z=r1z/nr1;
	//second segment
	r2x=chain_x[(i+2)%NbSegments]-chain_x[(i+1)%NbSegments];
	r2y=chain_y[(i+2)%NbSegments]-chain_y[(i+1)%NbSegments];
	r2z=chain_z[(i+2)%NbSegments]-chain_z[(i+1)%NbSegments];
	//normalize
	nr2=pitagora(r2x,r2y,r2z);
//  if(nr2<1e-6)cerr<<"nr2="<<nr2<<" i="<<i<<endl; 
	r2x=r2x/nr2; r2y=r2y/nr2; r2z=r2z/nr2;
	//dihedral vectors
	e1x=-ph_x[i%NbSegments]+th_x[i%NbSegments];
	e1y=-ph_y[i%NbSegments]+th_y[i%NbSegments];
	e1z=-ph_z[i%NbSegments]+th_z[i%NbSegments];
	
//  cout<<" th_x value  "<<th_x[i%NbSegments]<<endl;
//  cout<<" ph_x value  "<<ph_x[i%NbSegments]<<endl;
//  cout<<" E1X value  "<<e1x<<" E1Y value  "<<e1y<<" E1z value  "<<e1z<<" i  "<<i<<endl;
	e2x=-ph_x[(i+1)%NbSegments]+th_x[(i+1)%NbSegments];
	e2y=-ph_y[(i+1)%NbSegments]+th_y[(i+1)%NbSegments];
	e2z=-ph_z[(i+1)%NbSegments]+th_z[(i+1)%NbSegments];
	//project e1 on the plane normal to r1
	prod=e1x*r1x+e1y*r1y+e1z*r1z;//scalar prod e1*r1
	e1x=e1x-prod*r1x;
	e1y=e1y-prod*r1y;
	e1z=e1z-prod*r1z;
	//project e2 on the plane normal to r2
	prod=e2x*r2x+e2y*r2y+e2z*r2z;//scalar prod e2*r2
	e2x=e2x-prod*r2x;
	e2y=e2y-prod*r2y;
	e2z=e2z-prod*r2z;

	//normalize
	ne1=pitagora(e1x,e1y,e1z);
	ne2=pitagora(e2x,e2y,e2z);
        //  if(ne1<1e-6)cerr<<"ne1="<<ne1<<" i="<<i<<endl; 
        //  if(ne2<1e-6)cerr<<"ne2="<<ne2<<" i="<<i<<endl; 

	e1x=e1x/ne1; e1y=e1y/ne1; e1z=e1z/ne1;
	e2x=e2x/ne2; e2y=e2y/ne2; e2z=e2z/ne2;
	//normal vector to the plane containing the two segments
	nx=r1y*r2z-r1z*r2y;
	ny=r1z*r2x-r1x*r2z;
	nz=r1x*r2y-r1y*r2x;
	nn=pitagora(nx,ny,nz);
//  cout<<"ne1 "<<ne1<<" "<<"ne2 "<<ne2<<" nn "<<nn<<endl;


//////////////////  calculation of dihedral angles
	if(nn<0.0000000001)
		{//vectors parallel, see=scalar product, s=orientation
		see=e1x*e2x+e1y*e2y+e1z*e2z;
		s=threeprod(e1x,e1y,e1z,e2x,e2y,e2z,r1x,r1y,r1z);
 //   if(fabs(s)<1e-6)cerr<<"s="<<s<<" i="<<i<<endl; 
 //   if(fabs(see)>1-1e-6)cerr<<"see="<<see<<" i="<<i<<endl; 

		alpha=acos(see)*s/fabs(s);//alpha=fabs(fmod(acos(see),Pi))*s/fabs(s);
                cout<<see<<"normal vector was 0 "<<s<<" "<<alpha<<endl;
		}
	else
		{
		  nx=nx/nn;ny=ny/nn;nz=nz/nn;
		  se1=e1x*nx+e1y*ny+e1z*nz;
		  se2=e2x*nx+e2y*ny+e2z*nz;
		  s1=threeprod(e1x,e1y,e1z,nx,ny,nz,r1x,r1y,r1z);
		  //special case, e1 is parallel to n, then the angle is 0 or Pi
		  if(fabs(s1)<1e-10)	{a1=0; cout<<"normal vector parallel to a1-0"<<endl;
		    if(e1x*nx+e1y*ny+e1z*nz<0){a1=Pi; //cout<<"normal vector parallel to a1-Pi"<<endl;
		    }
		  }
		  //'general' case:
		  else a1=acos(se1)*s1/fabs(s1);//a1=fabs(fmod(acos(se1),Pi))*s1/fabs(s1);
		  
		  s2=threeprod(nx,ny,nz,e2x,e2y,e2z,r2x,r2y,r2z);//r2
		  if(fabs(s2)<1e-10) {a2=0;
		    if(e2x*nx+e2y*ny+e2z*nz<0){a2=Pi; //cout<<"normal vector parallel to a2-Pi"<<endl;
		    }
		  }
		  else a2=acos(se2)*s2/fabs(s2);//a2=fabs(fmod(acos(se2),Pi))*s2/fabs(s2);
		  
		  //    if(fabs(s1)<1e-10)cerr<<"s1="<<s1<<" i="<<i<<" a1="<<a1<<endl; 
		  // 		if(fabs(se1)>1-1e-10)cerr<<"se1="<<se1<<" i="<<i<<" a1="<<a1<<endl; 
		  //		if(fabs(s2)<1e-10)cerr<<"s2="<<s2<<" i="<<i<<" a2="<<a2<<endl; 
		  // 		if(fabs(se2)>1-1e-10)cerr<<"se2="<<se2<<" i="<<i<<" a2="<<a2<<endl; 
		  alpha=a1+a2;
//		  cout<<se1<<" "<<se2<<" nn was not 0 "<<s1<<" "<<s2<<" "<<alpha<<endl;


		}
	// rescale the angle so that it is between -Pi and Pi
	if(alpha>Pi){alpha=-2*Pi+alpha;}
	if(alpha<-Pi){alpha=2*Pi+alpha;}
	
   
//  cout<<"alpha calc:"<<alpha<<endl<<endl;
	twist=twist+alpha;
	// cout<<a1<<" a12 " <<a2<<" alpha " <<alpha<<" twist "<<twist/2/Pi<<endl;
	}

//  cout<<"library twist:"<<twist<<endl;
  return twist/(2.0*Pi);
}






int prol_asph(vector<double> & x,vector<double> & y,vector<double> & z, double &prolactness, double &asphericity, double &sq_radius_gyr){


  // let's calculate prolactness and asphericity
  // first we need to rearrange the vectors
  vector <vector <double> > data;
  data.push_back(x); data.push_back(y); data.push_back(z);
  int coord=3, N=x.size();
  double inertia[coord][coord];
  //set zero the inertia
  for(int i=0;i<coord;++i)
   for(int j=0;j<coord;++j)
    inertia[i][j]=0.0;

  
  /*   T_ij= (1/2N^2) SUM(n=1 to N) SUM(m=1 to N) (x_n^i - x_m^i)(x_n^j - x_m^j)  */

  //we can skip some calculation with c2=c because the matrix is simmetric
  //and because "dsyev" use only the lower or upper triangle of the matrix
  double CM[coord]; for(int i=0;i<coord;++i) CM[i]=0;

  // NEW ROUTINE 
  for(int i=0;i<N;i++)
   for(int c=0;c<coord;++c)
    CM[c]+=data[c][i];

  for(int c=0;c<coord;++c)
   CM[c]/=N;

  for(int c=0; c<coord;++c){
   for(int c1=c;c1<coord;++c1){
    double tc1=0;
    #pragma omp parallel for reduction(+:tc1)
    for(int i=0; i<N; ++i){
     tc1+=(CM[c]-data[c][i])*(CM[c1]-data[c1][i])/N;
    }
    inertia[c][c1]+=tc1;
   }
  }
  

  /* lapack accept lower or upper triangular matrices as 1D array.
  However, it seems that the standart it's reverse. Left/Lower triangular matrix
  are for LAPACK Upper, Right/Upper triangular matrix are Lower */

  //converting inertia to 1D inertia
  double inertia_1D[9];
  for(int i=0;i<coord;++i){
   for(int j=0;j<coord;++j){
    inertia_1D[i*coord+j]=inertia[i][j];
//   cout<<inertia[i][j]<<"  ";
   }
  }

  /* Locals */
  long int n = coord, lda = coord, info, lwork;
  double wkopt;
  double* work;
  /* Local arrays */
  double w[coord];

  /* Query and allocate the optimal workspace */
  lwork = -1;
  string nns="N";
  char *nn=new char[nns.length()+1];
  strcpy(nn,nns.c_str());

  string lls="L";
  char *ll=new char[lls.length()+1];
  strcpy(ll,lls.c_str());
  
  
  dsyev_( nn, ll, &n, inertia_1D, &lda, w, &wkopt, &lwork, &info );
  lwork = (long int)wkopt;
  work = (double*)malloc( lwork*sizeof(double) );


  /* Solve eigenproblem */
  dsyev_( nn, ll, &n, inertia_1D, &lda, w, work, &lwork, &info );

//  cout<<w[0]<<"  "<<w[1]<<"  "<<w[2]<<endl;
  /* Check for convergence */
  // something strange happen... sometimes even if the results are correct the info is > 0...
  // I will ignore this
//  if( info > 0 ) {
//   printf( "The algorithm failed to compute eigenvalues.\n" );
//   return 0;
//  }


  double a=pow(3*w[0],0.5), b=pow(3*w[1],0.5), c=pow(3*w[2],0.5);
  asphericity= (  (a-b)*(a-b)  + (a-c)*(a-c) + (b-c)*(b-c)  )/(2.0*(a+b+c)*(a+b+c));
  prolactness= (  (2.0*a-b-c)*(2.0*b-a-c)*(2.0*c-a-b)  )/( 2.0* pow( a*a + b*b + c*c -a*b -a*c -b*c, 1.5) );
  sq_radius_gyr= (a*a +b*b + c*c)/3.0;

  return 1;


}


void average_bead(vector <double> &matrix_d, int window, double &average, double &sigma){

int nbeads= pow(matrix_d.size(),0.5);

vector <double> storev; storev.resize(nbeads);

#pragma omp parallel for shared(storev)
for(int i=0; i<nbeads; ++i){
 vector <double> tempv;

// for(int j=i+window/2; j<i+nbeads-window/2; ++j){
//    int jres=j%nbeads;
//    tempv.push_back(matrix_d[i*nbeads+jres]);
// }
 
 for(int j=0; j<nbeads; ++j){
//  cout<<"i:"<<i<<"  j:"<<j<<endl;
  if(i<(window/2) || i>=(nbeads-window/2) ){
     if( j>=(i+window/2)%nbeads && j <=(i+nbeads-window/2)%nbeads  )
      tempv.push_back(matrix_d[i*nbeads+j]);
  }
  else{
     if(  j>(i+window/2) || j<(i-window/2)  ){
//      cout<<"for i:"<<i<<"  j:"<<j<<"  we are on the else"<<endl;
      tempv.push_back(matrix_d[i*nbeads+j]);
      }
  }
 }

// cout<<"size of tempv: "<<tempv.size()<<endl;
// cout<<"min el:"<<*min_element(tempv.begin(), tempv.end())<<endl<<endl<<endl;
 storev[i]=(*min_element(tempv.begin(), tempv.end()));
 tempv.clear();
}

 average_sigma(storev, average, sigma);

}

void average(vector <double> &avect, double &avev){
double sum = std::accumulate(avect.begin(), avect.end(), 0.0);
avev = sum / avect.size();
}


double median(vector <double> &a)
{
int size=a.size();
sort(a.begin(),a.begin()+size);
return size % 2 ? a[size/2] : (a[size/2 - 1] + a[size/2])/2.0;
}

void average_sigma(vector <double> &avect, double &avev, double&sigma){
average(avect,avev);

double sq_sum = inner_product(avect.begin(), avect.end(), avect.begin(), 0.0);
sigma = sqrt(sq_sum / avect.size() - avev * avev);
}


void orientation(double *data, int currentf, int total_last_frame, int num_atoms, int num_atoms_per_ring, int real_atoms_per_ring, int num_rings){
  
  int coord=3;
  double *Sab=new double[coord*coord]; for(int i=0; i< coord*coord;++i) Sab[i]=0;
  double inertia[coord][coord];
  string nns="V";
  char *nn=new char[nns.length()+1];
  strcpy(nn,nns.c_str());

  string lls="L";
  char *ll=new char[lls.length()+1];
  strcpy(ll,lls.c_str());

  /* Locals */
  long int n = coord, lda = coord, info, lwork;
  double wkopt;
  double* work;
  /* Local arrays */
  double w[coord];

  
  ifstream infile;
  ofstream outfile;
  infile.open("eigenvalues_rings.txt");
  if(!infile){
    infile.close(); infile.clear();
    outfile.open("eigenvalues_rings.txt");
    outfile<<"frame "<<endl;
    for(int i=0; i<num_rings; ++i)
     for(int k=0; k<coord; ++k)
      outfile<<"ring_"+num_to_str(i)+"_eigenvalue_"+num_to_str(k)+" ";
    outfile<<endl;
    outfile.close(); outfile.clear();
  }
  
  else infile.close(); infile.clear();
  
  outfile.open("eigenvalues_rings.txt",ios::app);
  outfile<<currentf+total_last_frame<<" ";

//  ofstream outfile4;
//  outfile4.open("largest_eigenv.txt", ios::app);
//  outfile4<<endl<<endl;
  
  
  for(int nr=0; nr<num_rings;++nr){
   vector <double> x, y, z;
   for(int pos=currentf*num_atoms*coord + nr*num_atoms_per_ring*coord; pos<currentf*num_atoms*coord + nr*num_atoms_per_ring*coord + real_atoms_per_ring*coord;pos+=coord){
   x.push_back(data[pos]);
   y.push_back(data[pos+1]);
   z.push_back(data[pos+2]);
   }
  
//  ofstream outfile2;
//  string rname="ring"+num_to_str(nr)+".txt";
//  outfile2.open(rname.c_str());
//  for(int jj=0; jj<x.size();++jj)
//   outfile2<<x[jj]<<" "<<y[jj]<<" "<<z[jj]<<endl;
//  outfile2.close();
  
  // first we need to rearrange the vectors
  vector <vector <double> > tdata;
  tdata.push_back(x); tdata.push_back(y); tdata.push_back(z);
  x.clear(); y.clear(); z.clear();
  

  //set zero the inertia
  for(int i=0;i<coord;++i)
   for(int j=0;j<coord;++j)
    inertia[i][j]=0.0;

  
  /*   T_ij= (1/2N^2) SUM(n=1 to N) SUM(m=1 to N) (x_n^i - x_m^i)(x_n^j - x_m^j)  */

  //we can skip some calculation with c2=c because the matrix is simmetric
  //and because "dsyev" use only the lower or upper triangle of the matrix
  double CM[coord]; for(int i=0;i<coord;++i) CM[i]=0;


  for(int i=0;i<real_atoms_per_ring;i++)
   for(int c=0;c<coord;++c)
    CM[c]+=tdata[c][i];

  for(int c=0;c<coord;++c)
   CM[c]/=real_atoms_per_ring;

  for(int c=0; c<coord;++c){
   for(int c1=c;c1<coord;++c1){
    double tc1=0;
    #pragma omp parallel for reduction(+:tc1)
    for(int i=0; i<real_atoms_per_ring; ++i){
     tc1+=(CM[c]-tdata[c][i])*(CM[c1]-tdata[c1][i])/real_atoms_per_ring;
    }
    inertia[c][c1]+=tc1;
   }
  }
  

  /* lapack accept lower or upper triangular matrices as 1D array.
  However, it seems that the standart it's reverse. Left/Lower triangular matrix
  are for LAPACK Upper, Right/Upper triangular matrix are Lower */

  //converting inertia to 1D inertia
  double inertia_1D[9];
  for(int i=0;i<coord;++i){
   for(int j=0;j<coord;++j){
    inertia_1D[i*coord+j]=inertia[i][j];
//   cout<<inertia[i][j]<<"  ";
   }
  }

  /* Query and allocate the optimal workspace */
  lwork = -1;
  dsyev_(nn, ll, &n, inertia_1D, &lda, w, &wkopt, &lwork, &info);
  lwork = (long int)wkopt;
  work = (double*)malloc( lwork*sizeof(double) );


//  ofstream outfile3;
//  string rname="all_details.txt";
//  outfile3.open(rname.c_str(), ios::app);
//  outfile3<<"inertia matrix"<<endl;
//  for(int i=0;i<coord;++i){
//   for(int j=0;j<coord;++j){
//    outfile3<<inertia_1D[i*coord +j]<<" ";
//    }
//   outfile3<<endl;
//  }


  /* Solve eigenproblem */
  dsyev_(nn, ll, &n, inertia_1D, &lda, w, work, &lwork, &info);
  // the maximum eigenvalue corresponding to the maximum eigenvector is stored in position [2][...]

  outfile<<w[0]<<" "<<w[1]<<" "<<w[2]<<" ";

//  outfile3<<"eigenvalues: ";
//  outfile3<<w[0]<<" "<<w[1]<<" "<<w[2]<<endl;
//  outfile3<<"eigenvectors: "<<endl;  
//  for(int i=0;i<coord;++i){
//   for(int j=0;j<coord;++j){
//    outfile3<<"i:"<<i<<"  j:"<<j<<"  "<<inertia_1D[i*coord +j]<<" ";
//    }
//   outfile3<<endl;
//  }
//  outfile3<<endl;
//  outfile3<<"largest eigenvector:"<<endl;

  vector <double> biggest_eigenv;
  for(int i=0;i<coord;++i){
   biggest_eigenv.push_back(inertia_1D[2*coord + i]);
//  outfile3<<inertia_1D[2*coord + i]<<" ";
//  outfile4<<inertia_1D[2*coord + i]<<" ";
  }
//  outfile3<<endl;
//  outfile4<<endl;
//  outfile3.close();
  
  double delta_k=0;
  for(int i=0;i<coord;++i){
   for(int k=0;k<coord;++k){
    if(i==k)
     delta_k=1.0/3.0;
    else
     delta_k=0.0;
    Sab[i*coord+k]+=(biggest_eigenv[i]*biggest_eigenv[k]-delta_k)/num_rings;
   }
  }
  
 }
// outfile4.close();
 
 outfile<<endl;
 outfile.close(); outfile.clear();
 
//   /* Query and allocate the optimal workspace */
// lwork = -1;
// dsyev_(nn, ll, &n, inertia_1D, &lda, w, &wkopt, &lwork, &info);
// lwork = (long int)wkopt;
// work = (double*)malloc( lwork*sizeof(double) );

// ofstream outfile3;
// string rname="Sab_matrix.txt";
// outfile3.open(rname.c_str(), ios::app);
// outfile3<<"order matrix"<<endl;
// for(int i=0;i<coord;++i){
//  for(int k=0;k<coord;++k){
//   outfile3<<Sab[i*coord+k]<<" ";
//  }
//  outfile3<<endl;
// }
//  outfile3.close();
 
 //we find again the diagonal value of the matrix Sab
 dsyev_(nn, ll, &n, Sab, &lda, w, work, &lwork, &info);

 infile.open("order_matrix.txt");
 if(!infile){
   infile.close(); infile.clear();
   outfile.open("order_matrix.txt");
   outfile<<"frame ";
    for(int k=0; k<coord; ++k)
     outfile<<"eigenvalue_"+num_to_str(k)+" ";
   outfile<<endl;
   outfile.close(); outfile.clear();
 }
  
 else infile.close(); infile.clear();
 
 outfile.open("order_matrix.txt",ios::app);
 outfile<<currentf+total_last_frame<<" ";
 for(int i=0;i<coord;++i)
  outfile<<w[i]<<" ";
 outfile<<endl;
 
 outfile.close(); outfile.clear();

 delete []work;
 delete []nn; delete []ll;

}






double get_linking(vector<double> &x,vector<double> &y,vector<double> &z, vector<double> &xs,vector<double> &ys,vector<double> &zs)
{
  double ZERO=1e-10;
  int N=x.size();
  int Ns=xs.size();

  double Linking=0;

  #pragma omp parallel for reduction(+:Linking)
  for(int i1=0;i1<N;i1++)
  {
    int j1,j2;
    double norme;
    double r13x,r13y,r13z;
    double r14x,r14y,r14z;
    double r24x,r24y,r24z;
    double r23x,r23y,r23z;
    double r34x,r34y,r34z;
    double r12x,r12y,r12z;
    double n1x,n1y,n1z;   
    double n2x,n2y,n2z;   
    double n3x,n3y,n3z;   
    double n4x,n4y,n4z;   
    double n5x,n5y,n5z;  

    j1=(i1+1)%N;

	  r12x=x[j1]-x[i1];
	  r12y=y[j1]-y[i1];
	  r12z=z[j1]-z[i1];

    for(int i2=0;i2<Ns;i2++)
	  {
	  j2=(i2+1)%Ns;
	  r13x=xs[i2]  -x[i1];
	  r13y=ys[i2]  -y[i1];
	  r13z=zs[i2]  -z[i1];
	  r14x=xs[j2]-x[i1];
	  r14y=ys[j2]-y[i1];
	  r14z=zs[j2]-z[i1];
	  r24x=xs[j2]-x[j1];
	  r24y=ys[j2]-y[j1];
	  r24z=zs[j2]-z[j1];
	  r23x=xs[i2]  -x[j1];
	  r23y=ys[i2]  -y[j1];
	  r23z=zs[i2]  -z[j1];

	  r34x=xs[j2]-xs[i2];
	  r34y=ys[j2]-ys[i2];
	  r34z=zs[j2]-zs[i2];


	  n1x=r13y*r14z-r13z*r14y; n1y=r13z*r14x-r13x*r14z; n1z=r13x*r14y-r13y*r14x;
	  norme=sqrt(n1x*n1x+n1y*n1y+n1z*n1z);
	  if(norme>ZERO){n1x=n1x/norme; n1y=n1y/norme; n1z=n1z/norme;}
		    
	      n2x=r14y*r24z-r14z*r24y; n2y=r14z*r24x-r14x*r24z; n2z=r14x*r24y-r14y*r24x;
	      norme=sqrt(n2x*n2x+n2y*n2y+n2z*n2z);
	      if(norme>ZERO){n2x=n2x/norme; n2y=n2y/norme; n2z=n2z/norme;}
		  
	      n3x=r24y*r23z-r24z*r23y; n3y=r24z*r23x-r24x*r23z; n3z=r24x*r23y-r24y*r23x;
	      norme=sqrt(n3x*n3x+n3y*n3y+n3z*n3z);
	      if(norme>ZERO){n3x=n3x/norme; n3y=n3y/norme; n3z=n3z/norme;}
		  
	      n4x=r23y*r13z-r23z*r13y; n4y=r23z*r13x-r23x*r13z; n4z=r23x*r13y-r23y*r13x;
	      norme=sqrt(n4x*n4x+n4y*n4y+n4z*n4z);
	      if(norme>ZERO){n4x=n4x/norme; n4y=n4y/norme; n4z=n4z/norme;}
		  
	      n5x=r34y*r12z-r34z*r12y; n5y=r34z*r12x-r34x*r12z; n5z=r34x*r12y-r34y*r12x;
	      if(n5x*r13x+n5y*r13y+n5z*r13z > ZERO)
		{
		    Linking += -(asin(Bound(n1x*n2x+n1y*n2y+n1z*n2z,-1,1))+asin(Bound(n2x*n3x+n2y*n3y+n2z*n3z,-1,1))+asin(Bound(n3x*n4x+n3y*n4y+n3z*n4z,-1,1))+asin(Bound(n4x*n1x+n4y*n1y+n4z*n1z,-1,1)))/(2*M_PI);
		}
	      else if (n5x*r13x+n5y*r13y+n5z*r13z < -ZERO)
		{
		    Linking += (asin(Bound(n1x*n2x+n1y*n2y+n1z*n2z,-1,1))+asin(Bound(n2x*n3x+n2y*n3y+n2z*n3z,-1,1))+asin(Bound(n3x*n4x+n3y*n4y+n3z*n4z,-1,1))+asin(Bound(n4x*n1x+n4y*n1y+n4z*n1z,-1,1)))/(2*M_PI);
		}		    
	    }
	}
  return -Linking/2;
}




#endif
