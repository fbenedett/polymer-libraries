#include "polymer_lib.h"
/*
Algorithm to measure the Accessible surface area

*/


void temp_MAP(long start, long end, long num_frames, long num_atoms, double epsilon, double *data, string namef);
int measure_sasa(vector <double> &x, vector <double> &y, vector <double> &z, double radius, double &sasa );


void temp_MAP(long start, long end, long num_frames, long num_atoms, double epsilon, double *data, string namef){
//contact map calculation
 ofstream outfile;

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



for(long i=0;i<map.size();++i)
    map[i].clear();
map.clear();

}


int measure_sasa(vector <double> &x, vector <double> &y, vector <double> &z, double radius, double &sasa ){

  int i;
  int npts = 500;
  double maxrad = -1;

#if 0
  // Query min/max atom radii for the entire molecule
  mol->get_radii_minmax(minrad, maxrad);
#endif

  //all the atoms have the same radius here we look at the distance "radius"
  maxrad= radius;
  
  // to generate the pair list we create a contact map with the optimized alghoritm, however we need to merge the vectors
  // in a single 1D array
  vector <double> data; data.resize(x.size()*3);
  for(int i=0;i<x.size();i+=3){
    data[i]=x[i];
    data[i+1]=y[i];
    data[i+2]=z[i];
  }
  

  string cmapn="map_sasa_"+num_to_str(radius)+".txt";
  vector < vector <int> > pairlist;
  temp_MAP(0, x.size(), 1, x.size(), radius, &data[0], cmapn);
  data.clear();
  

  // now read the file and create the pair list.
  // the pair list is a simple 2D array where we store the INDEX of the contacts (excluded a given bead with itself)
  ifstream infile;
  infile.open(cmapn.c_str());
  int val=0, counter=0, nump=0;
  vector <int> tempcont;
  while(infile>>val){
    if(nump!=counter && val!=0)
      tempcont.push_back(counter);
    
    counter++;
    
    if(counter==x.size()){
     counter=0;   
     nump++;
     pairlist.push_back(tempcont);
     tempcont.clear();
    }
  }
  infile.close();  infile.clear();  
  
  // we need to remove this contact map 
  

  // Seed the random number generator before each calculation.  This gives
  // reproducible results and still allows a more accurate answer to be
  // obtained by increasing the samples size.  I don't know if this is a
  // "good" seed value or not, I just picked something random-looking.
  srand(38572111);
  static const double RAND_MAX_INV = 1.0f/RAND_MAX;
  double PI=3.14159265358979323846264338;
  // All the spheres use the same random points.  
  double *spherepts = new double[3*npts];
  for (i=0; i<npts; i++) {
    double u1 = (double) rand();
    double u2 = (double) rand();
    double z = 2.0f*u1*RAND_MAX_INV -1.0f;
    double phi = (double) (2.0f*PI*u2*RAND_MAX_INV);
    double R = sqrtf(1.0f-z*z);
    spherepts[3*i  ] = R*cosf(phi);
    spherepts[3*i+1] = R*sinf(phi);
    spherepts[3*i+2] = z;
  }

  const double prefac = (double) (4 * PI / npts);
  double totarea = 0.0f;

  // compute area for each atom based on its pairlist
  #pragma omp parallel for reduction(+: totarea)
  for (int i=0; i<x.size(); ++i) {

      double loc[3];
      loc[0]= x[i]; loc[1]=y[i]; loc[2]=z[i];

      double rad = maxrad;
      double surfpos[3];
      int surfpts = npts;
      vector <int> nbrs = pairlist[i];
      for (int j=0; j<npts; j++) {
        surfpos[0] = loc[0] + rad*spherepts[3*j  ];
        surfpos[1] = loc[1] + rad*spherepts[3*j+1];
        surfpos[2] = loc[2] + rad*spherepts[3*j+2];
        int on = 1;
        for (int k=0; k<nbrs.size(); k++) {
          int ind = nbrs[k];
          
          double nbrloc[3];
          nbrloc[0]= x[ind]; nbrloc[1]= y[ind]; nbrloc[2]= z[ind];
          
          double radsq = maxrad; radsq *= radsq;
          double dx = surfpos[0]-nbrloc[0];
          double dy = surfpos[1]-nbrloc[1];
          double dz = surfpos[2]-nbrloc[2];
          if (dx*dx + dy*dy + dz*dz <= radsq){
            on = 0;
            break;
          }
        }
        if (on==0){
          surfpts--;
        }
      }
     double atomarea = prefac * rad * rad * surfpts;
     totarea += atomarea;

   }
 

   sasa = totarea;
   
   pairlist.clear();

   return 1;
}
