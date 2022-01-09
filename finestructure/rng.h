#ifndef __RNG_H__
#define __RNG_H__

#include <ctime>
#include<iostream>
#include <cstring>
#include <fstream>
#include <vector>
#include <random> 

#define ABS(x) ((x) >= 0 ? (x) : -(x))
#define mymin(a,b) ((a) <= (b) ? (a) : (b))
#define mymax(a,b) ((a) >= (b) ? (a) : (b))
#define OVERFLO 1e100
#define UNDERFLO 1e-100


namespace fines
{

  class my_rng
  {
  public:
    std::mt19937 rng;
    unsigned long seed;
  };

  my_rng * seedrng(unsigned long seed=0,int verbose=1);
  int saverng(my_rng * rng,std::string fname);
  int loadrng(my_rng * rng,std::string fname);
  void freerng(my_rng * rng);
  
  double RandomReal(my_rng * rng,double low, double high);
  int RandomInteger(my_rng * rng,int low, int high);
  double rnd(my_rng * rng);
  double RGamma(my_rng * rng,double n,double lambda);
  void RDirichlet(my_rng * rng,const std::vector<double> * a, std::vector<double>  * b);
  long RPoisson(my_rng * rng,double mu);
  double RExpon(my_rng * rng,double av);
  double RNormal(my_rng * rng,double mu,double sd) ;
  double sexpo(my_rng * rng);
  double snorm(my_rng * rng);
  double genexp(my_rng * rng,double av);   
  long ignpoi(my_rng * rng,double mean);  
  long ignuin(my_rng * rng,int low, int high);   
  double genunf(my_rng * rng,double low, double high);   
  long   Binomial(my_rng * rng,int n, double p);
  long   Binomial1(my_rng * rng,int n, double p);
  void LogRDirichlet (my_rng * rng,const double *a, const int k, double *b,double *c);

  double BinoProb(int n, double p,int i);
  double fsign(double num, double sign );

} // end namespace fines
#endif // __RNG_H__
