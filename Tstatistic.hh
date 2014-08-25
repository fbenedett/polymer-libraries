#ifndef STATISTIC_HH
#define STATISTIC_HH
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <sstream>
#include <vector>

using namespace std;
//************************************************************
//****************** Tstatistic **************************
//************************************************************

class Tstatistic//evalue X/Y
{
public:
  Tstatistic();
  //Tstatistic(const Tstatistic&);
  //~Tstatistic();
  //Tstatistic& operator =(const Tstatistic& stat);
  void Add(double X);
  void Reset();
  int GetNbStatistics(){return mean.size();}
  long long int GetBinLength(int k){return BinLength[k];}
  long long int GetNbBins(int k){return NbBins[k];}
  double GetAverage(int k=0){return mean[k];}
  double GetMeanSquare(int k){return  sqrt(M2[k])/(double)(NbBins[k]);}//=sqrt(Variance(Average)). relative error bars 1/sqrt(2(NbBins[k]-1))
  double GetCorrelation(int k){return  Correlations_M2[k]*NbBins[0]/(double)((NbBins[k]-1)*M2[0]);}//normalise par variance
  double Get_Distribution_MeanSquare(){return  sqrt(M2[0]/(double)(NbBins[0]));}//=sqrt(Variance(Samples))
  double GetIntegratedCorrelationTime(int k){return  M2[k]*((NbBins[0])*(NbBins[0]))/(double)(2.0*M2[0]*(NbBins[k])*(NbBins[k]));}//=sqrt(Variance(Average)). relative error bars sqrt(2/(NbBins[k]-1))
  int Get_koptimal();//koptimal= k which minimize: correlation + Relative error on the average errorbars = correlation + 1/sqrt(2*(NbBins-1))
  string Get_Average_MeanSquare_formated(){double x=GetAverage(),dx=GetMeanSquare(Get_koptimal()); stringstream str;str.precision(9);if(dx>1e-11){double prec=pow(10.0,round(log10(dx)))/100.0;str<<round(x/prec)*prec<<" "<<round(dx/prec)*prec;}else str<<x<<" "<<dx;return str.str();}

private:
  friend ostream& operator << (ostream& os, const Tstatistic& obj);
  friend istream& operator >> (istream& os,  Tstatistic& obj);

  long long int NbMesuresX;//nbre de  mesures effectuees jusqu'ici (NbMesures=NbBins*BinLength)
  vector <double> mean;
  vector <double> mean_tmp;
  vector <double> M2;
  vector <double> Correlations_mean1;
  vector <double> Correlations_mean2;
  vector <double> Correlations_M2;
  vector <double> Correlations_lastX;
  vector <long long int> BinLength;
  vector <long long int> NbBins;//nbre de bins mesures jusqu'ici

};

//************************************************************
//********************** Thistogram **************************
//************************************************************

class Thistogram//evalue X/Y
{
public:
  Thistogram();
  Thistogram(int Histo_NbInterval,double Histo_min,double Histo_max,int BufferSize=0);
  Thistogram(const Thistogram&);
  ~Thistogram();
  Thistogram& operator =(const Thistogram& stat);
  void Add(double X);
  int GetNbInterval(){return Histo_NbInterval;}
  double GetHisto_min(){return Histo_min;}
  double GetHisto_max(){return Histo_max;}
  void Reset();

  double GetArea();
  double Getx1(int n){return Histo_min+n*(Histo_max-Histo_min)/Histo_NbInterval;}
  double Getx2(int n){return Histo_min+(n+1)*(Histo_max-Histo_min)/Histo_NbInterval;}
  double Getx(int n){return Histo_min+(n+0.5)*(Histo_max-Histo_min)/Histo_NbInterval;}
  double GetHistogramAverage(int n){return stat_Histogram[n].GetAverage()*Histo_NbInterval/(double)(Histo_max-Histo_min);}
  double GetHistogramMeanSquare(int n){return stat_Histogram[n].GetMeanSquare(stat_Histogram[n].Get_koptimal())*Histo_NbInterval/(double)(Histo_max-Histo_min);}
  
  //underflow x<Histo_min
  double Getx1_underflow(){if(xmin<Histo_min)return xmin;else return Histo_min;}
  double Getx2_underflow(){return Histo_min;}
  double Getx_underflow(){if(xmin<Histo_min)return (xmin+Histo_min)/2.0;else return Histo_min;}
  double GetHistogramAverage_underflow(){if(xmin<Histo_min)return stat_Histogram[Histo_NbInterval].GetAverage()/(Histo_min-xmin);else return 0;}
  double GetHistogramMeanSquare_underflow(){if(xmin<Histo_min)return stat_Histogram[Histo_NbInterval].GetMeanSquare(stat_Histogram[Histo_NbInterval].Get_koptimal())/(Histo_min-xmin);else return 0;}

  //overflow x>Histo_max
  double Getx1_overflow(){return Histo_max;}
  double Getx2_overflow(){if(xmax>Histo_max)return xmax;else return Histo_max;}
  double Getx_overflow(){if(xmax>Histo_max)return (xmax+Histo_max)/2.0;else return Histo_max;}
  double GetHistogramAverage_overflow(){if(xmax>Histo_max)return stat_Histogram[Histo_NbInterval+1].GetAverage()/(xmax-Histo_max);else return 0;}
  double GetHistogramMeanSquare_overflow(){if(xmax>Histo_max)return stat_Histogram[Histo_NbInterval+1].GetMeanSquare(stat_Histogram[Histo_NbInterval+1].Get_koptimal())/(xmax-Histo_max);else return 0;}

private:
  friend ostream& operator << (ostream& os, const Thistogram& obj);
  friend istream& operator >> (istream& os,  Thistogram& obj);
  
  int Histo_NbInterval;
  double Histo_min;
  double Histo_max;
  double xmin;
  double xmax;
  int BufferSize;//Nb samples collected before sending the histogram to Tstatistics.
  int samplecount;//Nb samples collected.
  double* HistoTemp;
  Tstatistic* stat_Histogram;
};

#endif

