//************************************************************
//****************** Tstatistic **************************
//************************************************************

//---------------------------------------------------------------------------

#include "Tstatistic.hh"

Tstatistic::Tstatistic()
{ 
#ifdef  DEBUG_FLAG
  cout<<"constructeur Tstatistic"<<endl;
#endif
  NbMesuresX=0;

}

//---------------------------------------------------------------------------
//Tstatistic::Tstatistic(const Tstatistic& stat)
//{
//  cout<<"constructeur copie Tstatistic"<<endl;
//
//  if(ArraySize!=stat.ArraySize)
//    {
//      ArraySize=stat.ArraySize;
//      mean.clear();
//      mean_tmp.clear();
//      M2.clear();
//      Correlations_mean1.clear();
//      Correlations_mean2.clear();
//      Correlations_M2.clear();
//      Correlations_lastX.clear();
//      BinLength.clear();
//      NbBins.clear();
//
//      try{mean.resize(ArraySize);}catch(...){cerr<<"Tstatistic::Tstatistic(const Tstatistic&): can't create mean. Aborting..."<<endl;exit(1);}
//      try{mean_tmp.resize(ArraySize);}catch(...){cerr<<"Tstatistic::Tstatistic(const Tstatistic&): can't create mean_tmp. Aborting..."<<endl;exit(1);}
//      try{M2.resize(ArraySize);}catch(...){cerr<<"Tstatistic::Tstatistic(const Tstatistic&): can't create M2. Aborting..."<<endl;exit(1);}
//      try{Correlations_mean1.resize(ArraySize);}catch(...){cerr<<"Tstatistic::Tstatistic(const Tstatistic&): can't create Correlations_mean1. Aborting..."<<endl;exit(1);}
//      try{Correlations_mean2.resize(ArraySize);}catch(...){cerr<<"Tstatistic::Tstatistic(const Tstatistic&): can't create Correlations_mean2. Aborting..."<<endl;exit(1);}
//      try{Correlations_M2.resize(ArraySize);}catch(...){cerr<<"Tstatistic::Tstatistic(const Tstatistic&): can't create Correlations_M2. Aborting..."<<endl;exit(1);}
//      try{Correlations_lastX.resize(ArraySize);}catch(...){cerr<<"Tstatistic::Tstatistic(const Tstatistic&): can't create Correlations_lastX. Aborting..."<<endl;exit(1);}
//      try{BinLength.resize(ArraySize);}catch(...){cerr<<"Tstatistic::Tstatistic(const Tstatistic&): can't create BinLength. Aborting..."<<endl;exit(1);}
//      try{NbBins.resize(ArraySize);}catch(...){cerr<<"Tstatistic::Tstatistic(const Tstatistic&): can't create NbBins. Aborting..."<<endl;exit(1);}
//    }  
//
//  NbStatistics=stat.NbStatistics;
//  NbMesuresX=stat.NbMesuresX;
//  for(int k=0;k<ArraySize;k++)
//    {
//      mean[k]=stat.mean[k];
//      mean_tmp[k]=stat.mean_tmp[k];
//      M2[k]=stat.M2[k];
//      Correlations_mean1[k]=stat.Correlations_mean1[k];
//      Correlations_mean2[k]=stat.Correlations_mean2[k];
//      Correlations_M2[k]=stat.Correlations_M2[k];
//      Correlations_lastX[k]=stat.Correlations_lastX[k];
//      BinLength[k]=stat.BinLength[k];
//      NbBins[k]=stat.NbBins[k];
//    }
//}

//---------------------------------------------------------------------------

//Tstatistic::~Tstatistic()
//{
//#ifdef  DEBUG_FLAG
//  cout<<"destructeur Tstatistic"<<endl;
//#endif
//      mean.clear();
//      mean_tmp.clear();
//      M2.clear();
//      Correlations_mean1.clear();
//      Correlations_mean2.clear();
//      Correlations_M2.clear();
//      Correlations_lastX.clear();
//      BinLength.clear();
//      NbBins.clear();
//
//}

//---------------------------------------------------------------------------

//Tstatistic& Tstatistic::operator=(const Tstatistic& stat)//event: return *this
//{
//  NbMesuresX=stat.NbMesuresX;
//      mean=stat.mean;
//      mean_tmp=stat.mean_tmp;
//      M2=stat.M2;
//      Correlations_mean1=stat.Correlations_mean1;
//      Correlations_mean2=stat.Correlations_mean2;
//      Correlations_M2=stat.Correlations_M2;
//      Correlations_lastX=stat.Correlations_lastX;
//      BinLength=stat.BinLength;
//      NbBins=stat.NbBins;
//
//  return *this;
//}

//---------------------------------------------------------------------------

void Tstatistic::Add(double X)
{
 
  NbMesuresX++;

  //////////////////check if NbStatistics or ArraySize should increase ////////////////
  if(mean.size()==0)
    {
      long long int lbin=1;
      BinLength.push_back(lbin);
      NbBins.push_back(0);
      mean.push_back(0);
      mean_tmp.push_back(0);
      M2.push_back(0);
      Correlations_mean1.push_back(0);
      Correlations_mean2.push_back(0);
      Correlations_M2.push_back(0);
      Correlations_lastX.push_back(0);
      lbin=2*lbin;
      
    }
  else if(NbMesuresX%BinLength.back()==0)//need to increment NbStatistics
    {
      mean.push_back(0);
      mean_tmp.push_back((mean[mean.size()-2]+mean_tmp[mean.size()-2])/2.0);//1/2.0=BinLength[NbStatistics-1]/BinLength[NbStatistics];
      M2.push_back(0);
      Correlations_lastX.push_back(0);
      Correlations_mean1.push_back(0);
      Correlations_mean2.push_back(0);
      Correlations_M2.push_back(0);
      BinLength.push_back(2*BinLength[mean.size()-2]);
      NbBins.push_back(0);
    }
  //version with correlation=correlation between bins (mean_tmp instead of X)
  ////////////////////////// Add X to the stat /////////////////
  for(int k=0;k<mean.size();k++)
    {
      mean_tmp[k]+=X/BinLength[k];
      if(NbMesuresX%BinLength[k]==0)
	{
	  double delta;
	  NbBins[k]++;
	  if(NbBins[k]>1)
	    {
	      Correlations_M2[k]=Correlations_M2[k]+(Correlations_lastX[k]-Correlations_mean1[k])*(mean_tmp[k]-Correlations_mean2[k]);
	      
	      delta=mean_tmp[k]-Correlations_mean2[k];
	      Correlations_mean2[k]=Correlations_mean2[k]+delta/(double)(NbBins[k]-1);
	    }
	  delta=mean_tmp[k]-Correlations_mean1[k];
	  Correlations_mean1[k]=Correlations_mean1[k]+delta/(double)NbBins[k];
	  
	  delta=mean_tmp[k]-mean[k];
	  mean[k]=mean[k]+delta/(double)NbBins[k];
	  M2[k]=M2[k]+delta*(mean_tmp[k]-mean[k]);
	  
	  Correlations_lastX[k]=mean_tmp[k];
	  mean_tmp[k]=0;
	}
    }
  
}

//---------------------------------------------------------------------------

void Tstatistic::Reset()
{
  //ArraySize=1;//ATTENTION: we need ArraySize>=NbStatistics+1
  NbMesuresX=0;

      mean.clear();
      mean_tmp.clear();
      M2.clear();
      Correlations_mean1.clear();
      Correlations_mean2.clear();
      Correlations_M2.clear();
      Correlations_lastX.clear();
      BinLength.clear();
      NbBins.clear();

}

//---------------------------------------------------------------------------

int Tstatistic::Get_koptimal()
{
//koptimal= k which minimize: correlation + Relative error on the average errorbars 
//= correlation + 1/sqrt(2*(NbBins-1))
//= Correlations_M2[k]*NbBins[0]/(double)((NbBins[k]-1)*M2[0]) + 1/sqrt(2*(NbBins[k]-1))
  double min=Correlations_M2[0]*NbBins[0]/(double)((NbBins[0]-1)*M2[0]) + 1/sqrt(2*(NbBins[0]-1)),val;
  int kopt=0;
  for(int k=1;(k<mean.size())&&(NbBins[k]>1);k++)
    {
      val=Correlations_M2[k]*NbBins[0]/(double)((NbBins[k]-1)*M2[0]) + 1/sqrt(2*(NbBins[k]-1));
      if(min>val)
	{
	  min=val;
	  kopt=k;
	}
    }
  return kopt;
}
//---------------------------------------------------------------------------

ostream& operator << (ostream& os, const Tstatistic& obj)
{
  os<<obj.mean.size()<<" ";

  for(int k=0;k<obj.mean.size();k++)
    os<<obj.mean[k]<<" ";
  for(int k=0;k<obj.mean.size();k++)
    os<<obj.mean_tmp[k]<<" ";
  for(int k=0;k<obj.mean.size();k++)
    os<<obj.M2[k]<<" ";
  for(int k=0;k<obj.mean.size();k++)
    os<<obj.Correlations_mean1[k]<<" ";
  for(int k=0;k<obj.mean.size();k++)
    os<<obj.Correlations_mean2[k]<<" ";
  for(int k=0;k<obj.mean.size();k++)
    os<<obj.Correlations_M2[k]<<" ";
  for(int k=0;k<obj.mean.size();k++)
    os<<obj.Correlations_lastX[k]<<" ";
  for(int k=0;k<obj.mean.size();k++)
    os<<obj.BinLength[k]<<" ";
  for(int k=0;k<obj.mean.size();k++)
    os<<obj.NbBins[k]<<" ";

  os<<obj.NbMesuresX;

  return os;
}
 
//---------------------------------------------------------------------------

istream& operator >> (istream& os,  Tstatistic& obj)
{
  string strtemp;
  int size;
  if(!os.good()){cerr<<"Tstatistic::Tstatistic(filename): Error with file 1"<<endl;exit(1);}
  os>>size;
  if(!os.good()){cerr<<"Tstatistic::Tstatistic(filename): Error with file 2"<<endl;exit(1);}
  obj.mean.resize(size);
  obj.mean_tmp.resize(size);
  obj.M2.resize(size);
  obj.Correlations_mean1.resize(size);
  obj.Correlations_mean2.resize(size);
  obj.Correlations_M2.resize(size);
  obj.Correlations_lastX.resize(size);
  obj.BinLength.resize(size);
  obj.NbBins.resize(size);

  for(int k=0;k<size;k++)
    {
      os>>obj.mean[k];
      if(!os.good()){cerr<<"Tstatistic::Tstatistic(filename): nan or inf in obj.mean[k] k="<<k<<" size="<<size<<endl;os.clear();os>>strtemp;obj.mean[k]=strtod(strtemp.c_str(),NULL);}//in case of nan or inf;
      if(!os.good()){cerr<<"Tstatistic::Tstatistic(filename): Error with file 5 "<<k<<" size="<<size<<endl;exit(1);}
    }
  for(int k=0;k<size;k++)
    {
      os>>obj.mean_tmp[k];
      if(!os.good()){cerr<<"Tstatistic::Tstatistic(filename): nan or inf in obj.mean_tmp[k] k="<<k<<" size="<<size<<endl;os.clear();os>>strtemp;obj.mean_tmp[k]=strtod(strtemp.c_str(),NULL);}//in case of nan or inf;
      if(!os.good()){cerr<<"Tstatistic::Tstatistic(filename): Error with file 6 "<<k<<" size="<<size<<endl;exit(1);}
    }
  for(int k=0;k<size;k++)
    {
      os>>obj.M2[k];
      if(!os.good()){cerr<<"Tstatistic::Tstatistic(filename): nan or inf in obj.M2[k] k="<<k<<" size="<<size<<endl;os.clear();os>>strtemp;obj.M2[k]=strtod(strtemp.c_str(),NULL);}//in case of nan or inf;
      if(!os.good()){cerr<<"Tstatistic::Tstatistic(filename): Error with file 7 "<<k<<" size="<<size<<endl;exit(1);}
    }
  for(int k=0;k<size;k++)
    {
      os>>obj.Correlations_mean1[k];
      if(!os.good()){cerr<<"Tstatistic::Tstatistic(filename): nan or inf in obj.Correlations_mean1[k] k="<<k<<" size="<<size<<endl;os.clear();os>>strtemp;obj.Correlations_mean1[k]=strtod(strtemp.c_str(),NULL);}//in case of nan or inf;
      if(!os.good()){cerr<<"Tstatistic::Tstatistic(filename): Error with file 8 "<<k<<" size="<<size<<endl;exit(1);}
    }
  for(int k=0;k<size;k++)
    {
      os>>obj.Correlations_mean2[k];
      if(!os.good()){cerr<<"Tstatistic::Tstatistic(filename): nan or inf in obj.Correlations_mean2[k] k="<<k<<" size="<<size<<endl;os.clear();os>>strtemp;obj.Correlations_mean2[k]=strtod(strtemp.c_str(),NULL);}//in case of nan or inf;
      if(!os.good()){cerr<<"Tstatistic::Tstatistic(filename): Error with file 9 "<<k<<" size="<<size<<endl;exit(1);}
    }
  for(int k=0;k<size;k++)
    {
      os>>obj.Correlations_M2[k];
      if(!os.good()){cerr<<"Tstatistic::Tstatistic(filename): nan or inf in obj.Correlations_M2[k] k="<<k<<" size="<<size<<endl;os.clear();os>>strtemp;obj.Correlations_M2[k]=strtod(strtemp.c_str(),NULL);}//in case of nan or inf;
      if(!os.good()){cerr<<"Tstatistic::Tstatistic(filename): Error with file 10 "<<k<<" size="<<size<<endl;exit(1);}
    }
  for(int k=0;k<size;k++)
    {
      os>>obj.Correlations_lastX[k];
      if(!os.good()){cerr<<"Tstatistic::Tstatistic(filename): nan or inf in obj.Correlations_lastX[k] k="<<k<<" size="<<size<<endl;os.clear();os>>strtemp;obj.Correlations_lastX[k]=strtod(strtemp.c_str(),NULL);}//in case of nan or inf;
      if(!os.good()){cerr<<"Tstatistic::Tstatistic(filename): Error with file 11 "<<k<<" size="<<size<<endl;exit(1);}
    }
  for(int k=0;k<size;k++)
    {
      os>>obj.BinLength[k];
      if(!os.good()){cerr<<"Tstatistic::Tstatistic(filename): Error with file 12 "<<k<<" size="<<size<<endl;exit(1);}
    }
  for(int k=0;k<size;k++)
    {
      os>>obj.NbBins[k];
     if(!os.good()){cerr<<"Tstatistic::Tstatistic(filename): Error with file 13 "<<k<<" size="<<size<<endl;exit(1);}
    }
  
  os>>obj.NbMesuresX;
  if(!os.good()){cerr<<"Tstatistic::Tstatistic(filename): Error with file 4"<<endl;exit(1);}

  return os;
}
 
//************************************************************
//****************** Thistogram **************************
//************************************************************

Thistogram::Thistogram()
{ 
  Histo_NbInterval=1;
  Histo_min=0;
  Histo_max=1;
  xmin=1e200;
  xmax=-1e200;
  BufferSize=1;
  samplecount=0;
  try{HistoTemp=new double[Histo_NbInterval+2];}catch(...){cerr<<"Thistogram::Thistogram: can't create HistoTemp. Aborting..."<<endl;exit(1);}
  try{stat_Histogram=new Tstatistic[Histo_NbInterval+2];}catch(...){cerr<<"Thistogram::Thistogram: can't create Thistogram. Aborting..."<<endl;exit(1);}
 
  for(int n=0;n<Histo_NbInterval+2;n++)
    {
      HistoTemp[n]=0;
    }
}
//---------------------------------------------------------------------------

Thistogram::Thistogram(int histo_nbinterval,double histo_min,double histo_max,int buffersize)
{ 
  Histo_NbInterval=histo_nbinterval;
  Histo_min=histo_min;
  Histo_max=histo_max;
  xmin=1e200;
  xmax=-1e200;
  BufferSize=buffersize;
  samplecount=0;
  try{HistoTemp=new double[Histo_NbInterval+2];}catch(...){cerr<<"Thistogram::Thistogram: can't create HistoTemp. Aborting..."<<endl;exit(1);}
  try{stat_Histogram=new Tstatistic[Histo_NbInterval+2];}catch(...){cerr<<"Thistogram::Thistogram: can't create Thistogram. Aborting..."<<endl;exit(1);}

  for(int n=0;n<Histo_NbInterval+2;n++)
    {
      HistoTemp[n]=0;
    }
}

//---------------------------------------------------------------------------
Thistogram::Thistogram(const Thistogram& stat)
{
  cout<<"constructeur copie Thistogram not implemented."<<endl;
  exit(1);
}

//---------------------------------------------------------------------------

Thistogram::~Thistogram()
{
  delete[] HistoTemp;
  delete[] stat_Histogram;
}

//---------------------------------------------------------------------------

Thistogram& Thistogram::operator=(const Thistogram& obj)//event: return *this
{

  Histo_min=obj.Histo_min;
  Histo_max=obj.Histo_max;
  xmin=obj.xmin;
  xmax=obj.xmax;
  BufferSize=obj.BufferSize;
  samplecount=obj.samplecount;
  if(Histo_NbInterval!=obj.Histo_NbInterval)
    {
      Histo_NbInterval=obj.Histo_NbInterval;
      delete[] HistoTemp;
      delete[] stat_Histogram;
      try{HistoTemp=new double[Histo_NbInterval+2];}catch(...){cerr<<"Thistogram::Thistogram: can't create HistoTemp. Aborting..."<<endl;exit(1);}
      try{stat_Histogram=new Tstatistic[Histo_NbInterval+2];}catch(...){cerr<<"Thistogram::Thistogram: can't create Thistogram. Aborting..."<<endl;exit(1);}
    }
  for(int i=0;i<Histo_NbInterval+2;i++)
    {
      HistoTemp[i]=obj.HistoTemp[i];
      stat_Histogram[i]=obj.stat_Histogram[i];
    }
  return *this;
}

//---------------------------------------------------------------------------

void Thistogram::Add(double X)
{
  int n;
  if(X<xmin)xmin=X;
  if(X>xmax)xmax=X;
  if(BufferSize<2)
    {
      n=lround(floor((X-Histo_min)/(Histo_max-Histo_min)*Histo_NbInterval));
      if(n>=Histo_NbInterval)n=Histo_NbInterval+1;
      if(n<0)n=Histo_NbInterval;
      for(int i=0;i<Histo_NbInterval+2;i++)
	if(n==i)
	  stat_Histogram[i].Add(1);
	else
	  stat_Histogram[i].Add(0);
    }
  else
    {
      samplecount++; 
      n=lround(floor((X-Histo_min)/(Histo_max-Histo_min)*Histo_NbInterval));
      if(n>=Histo_NbInterval)n=Histo_NbInterval+1;
      if(n<0)n=Histo_NbInterval;
      HistoTemp[n]++;
      if(samplecount==BufferSize)
	{
	  for(int n=0;n<Histo_NbInterval+2;n++)
	    {
	      stat_Histogram[n].Add(HistoTemp[n]/(double)BufferSize);
	      HistoTemp[n]=0;
	    }
	  samplecount=0;
	}
    }
}


//---------------------------------------------------------------------------

void Thistogram::Reset()
{
  xmin=1e200;
  xmax=-1e200;
  samplecount=0;
  for(int i=0;i<Histo_NbInterval+2;i++)
    {
      HistoTemp[i]=0;
      stat_Histogram[i].Reset();
    }
  
}

//---------------------------------------------------------------------------
double Thistogram::GetArea()
{
     double area=0,dx;
     dx=Histo_min-xmin;
     if(dx>0)
       area+=stat_Histogram[Histo_NbInterval].GetAverage();
     for(int l=0;l<Histo_NbInterval;l++)
       {
	 dx=(Histo_max-Histo_min)/Histo_NbInterval;
	 area+=stat_Histogram[l].GetAverage();
       }
     dx=xmax-Histo_max;
     if(dx>0)
       area+=stat_Histogram[Histo_NbInterval+1].GetAverage();
     return area;
}
//---------------------------------------------------------------------------

ostream& operator << (ostream& os, const Thistogram& obj)
{
  os<<obj.Histo_NbInterval<<" ";
  os<<obj.Histo_min<<" ";
  os<<obj.Histo_max<<" ";
  os<<obj.xmin<<" ";
  os<<obj.xmax<<" ";
  for(int i=0;i<obj.Histo_NbInterval+2;i++)
    os<<obj.HistoTemp[i]<<" ";
  for(int i=0;i<obj.Histo_NbInterval+2;i++)
    os<<obj.stat_Histogram[i]<<" ";
  os<<obj.samplecount<<" ";
  os<<obj.BufferSize;
  
  return os;
}
 
//---------------------------------------------------------------------------

istream& operator >> (istream& os,  Thistogram& obj)
{
  string strtemp;
  int  Histo_NbInterval_tmp;
  if(!os.good()){cerr<<"Thistogram::Thistogram(filename): Error with file 1"<<endl;exit(1);}
  os>>Histo_NbInterval_tmp;
  if(!os.good()){cerr<<"Thistogram::Thistogram(filename): Error with file 2"<<endl;exit(1);}

  os>>obj.Histo_min;
  if(!os.good()){cerr<<"Thistogram::Thistogram(filename): nan or inf in obj.Histo_min"<<endl;os.clear();os>>strtemp;obj.Histo_min=strtod(strtemp.c_str(),NULL);}//in case of nan or inf;
  if(!os.good()){cerr<<"Thistogram::Thistogram(filename): Error with file 3 "<<endl;exit(1);}

  os>>obj.Histo_max;
  if(!os.good()){cerr<<"Thistogram::Thistogram(filename): nan or inf in obj.Histo_max"<<endl;os.clear();os>>strtemp;obj.Histo_max=strtod(strtemp.c_str(),NULL);}//in case of nan or inf;
  if(!os.good()){cerr<<"Thistogram::Thistogram(filename): Error with file 4 "<<endl;exit(1);}

  os>>obj.xmin;
  if(!os.good()){cerr<<"Thistogram::Thistogram(filename): nan or inf in obj.xmin"<<endl;os.clear();os>>strtemp;obj.xmin=strtod(strtemp.c_str(),NULL);}//in case of nan or inf;
  if(!os.good()){cerr<<"Thistogram::Thistogram(filename): Error with file 5"<<endl;exit(1);}

  os>>obj.xmax;
  if(!os.good()){cerr<<"Thistogram::Thistogram(filename): nan or inf in obj.xmax"<<endl;os.clear();os>>strtemp;obj.xmax=strtod(strtemp.c_str(),NULL);}//in case of nan or inf;
  if(!os.good()){cerr<<"Thistogram::Thistogram(filename): Error with file 6"<<endl;exit(1);}

  if(Histo_NbInterval_tmp!=obj.Histo_NbInterval)
    {
      obj.Histo_NbInterval=Histo_NbInterval_tmp;
      delete[] obj.HistoTemp;
      delete[] obj.stat_Histogram;
      try{obj.HistoTemp=new double[obj.Histo_NbInterval+2];}catch(...){cerr<<"Thistogram::Thistogram: can't create HistoTemp. Aborting..."<<endl;exit(1);}
      try{obj.stat_Histogram=new Tstatistic[obj.Histo_NbInterval+2];}catch(...){cerr<<"Thistogram::Thistogram: can't create Thistogram. Aborting..."<<endl;exit(1);}
    }

  for(int i=0;i<obj.Histo_NbInterval+2;i++)
    {
      os>>obj.HistoTemp[i];
      if(!os.good()){cerr<<"Thistogram::Thistogram(filename): nan or inf in obj.HistoTemp[i] i="<<i<<endl;os.clear();os>>strtemp;obj.HistoTemp[i]=strtod(strtemp.c_str(),NULL);}//in case of nan or inf;
      if(!os.good()){cerr<<"Thistogram::Thistogram(filename): Error with file 7"<<endl;exit(1);}
    }
  for(int i=0;i<obj.Histo_NbInterval+2;i++)
    {
      if(!os.good()){cerr<<"Thistogram::Thistogram(filename): Error with file 8 i="<<i<<endl;exit(1);}
      os>>obj.stat_Histogram[i];
    }
  if(!os.good()){cerr<<"Thistogram::Thistogram(filename): Error with file 9"<<endl;exit(1);}
  os>>obj.samplecount;
  if(!os.good()){cerr<<"Thistogram::Thistogram(filename): Error with file 10"<<endl;exit(1);}
  os>>obj.BufferSize;
  if(!os.good()){cerr<<"Thistogram::Thistogram(filename): Error with file 11"<<endl;exit(1);}

  return os;
}


