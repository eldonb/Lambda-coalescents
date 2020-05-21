#include <iostream>
#include <vector>
#include <random>
#include <functional>
#include <memory>
#include <utility>
#include <algorithm>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <list>
#include <string>
#include <fstream>
#include <forward_list>
#include <assert.h>
#include <math.h>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/beta_distribution.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <boost/random/lagged_fibonacci.hpp>
#include <boost/random/uniform_01.hpp>
#include <Rcpp.h>
using namespace Rcpp;
using namespace std;
using namespace boost;
using namespace boost::random;
using namespace math;


 /* obtain a seed out of thin air for the random number engine */
  std::random_device randomseed;
  /* Standard mersenne twister  random number engine seeded with rng() */
  std::mt19937_64 rng(randomseed());
 /* construct a uniform on the unit interval; call with uni01( rng ) */
std::uniform_real_distribution<double> uni01(0.0, 1.0);
/* construct a beta distribution */
   boost::random::beta_distribution<double> randombeta(1.0, 1.0);
/* constuct a binomial */
std::binomial_distribution<unsigned> randombinom(1, .5) ;


static double xreject( const double a, const double u, const double M)
{
/* need boost library */
 double x;

 randombeta =  boost::random::beta_distribution<double>(2.-a, a);
  x = randombeta(rng); 
  /* x =  boost::math::ibeta_inv(2.-a, a, u); */

 while( (x > M) || ( x <= 0.) ){
    x = randombeta(rng); }
 /* x =   boost::math::ibeta_inv(2.-a, a, u) ;} */


  assert( x > 0. );
  assert( x <= M);
  return( x);
}


static double itime( const double n, const double a,   const double M,  const double b,   double* x, double* Sjp)
{
  /* $n$ is  number of  blocks; $a$ is $\alpha$; $M=M$, $b$ is exponential growth parameter */
double timi = 0. ;
double ct = n*(n - 1.)/2. ;

/*
  timi = timi + ( -log(1. - uni01(rng) )/(ct) ); 
 x[0] = xreject( a, uni01(rng), M); 
*/

if ( a < 2. ){ 
do{
  timi = timi + ( -log(1. - uni01(rng) )/(ct) ); 
  x[0] = xreject( a, uni01(rng), M); }
  while ( uni01(rng) >  (((1. - pow(1.-x[0], n) - (n*x[0]*pow(1.-x[0], n-1.)))/(x[0]*x[0]))/(ct) ) );
}
else{
   timi =  ( -log(1. - uni01(rng) )/( ct ) );} 

  double Sj = (a < 2 ? 0. : (b > 0. ? log( exp(b*Sjp[0]) - (b*log(uni01(rng)) / ( ct )))/b : Sjp[0] + timi ) ); 

  double svar =   (a < 2. ? timi : Sj - Sjp[0]) ;
  Sjp[0] = Sj ;
  return( svar ); 
}

static unsigned int samplek( const double n,  const double x )
{

   randombinom = std::binomial_distribution<unsigned>(n-2., x) ;
   unsigned  k = 2 +  randombinom(rng)  ;
   while( uni01(rng) > 2./( (double)(k*(k - 1)) ) ){
      k = 2 +  randombinom( rng ); }
 
 assert( k > 1);
 assert( k <= (unsigned)n);
 
  return(k);
}


// [[Rcpp::export]]

void lambdaijk( const double n, const double a, const double K, const double b, const double theta, const unsigned Imin, const unsigned Imax, const unsigned Jmin, const unsigned Jmax,  const unsigned Kmin, const unsigned Kmax,  const unsigned trials, const unsigned numer )
{
  /*   NumericMatrix& res */
 /* std::cout << n << ' ' << a << ' ' << K << ' ' << b << '\n' ; */
  /* set dimenions of matrix `res` 
  res.attr("dim") = Dimension(trials, 3) ;
  */  

  std::string resskra  = "resskra";
  resskra = resskra + std::to_string(numer);
  
  std::vector<unsigned> v;
  std::vector<double> ebi;
  std::list<unsigned> tre;
  std::list<unsigned>::iterator it;

  unsigned r, size, kblocks, index ;
  double   timi, z;
  /*
   const double   a = atof(argv[2]);
   const double K = atof(argv[3]);
   const double b = atof(argv[4]); 
   const double theta = atof(argv[5]);
   */
   const double M = (a < 2. ? (K > 0. ? K/(K + 1. + (pow(2,1.-a)/(a-1.))) : 1.) : 1.);
   double *x = (double *)calloc(1, sizeof(double)); 
   double *Sjp = (double *)calloc(1, sizeof(double)); 

   /* generating a poisson random number generator */
    std::poisson_distribution<unsigned> rpois(1.0) ;
    
   /* generating a beta random number generator using boost */

  /* generating an exponential dist */

  /* initialise spectrum */
    /* ebi.assign( (unsigned)n  + 1, 0. );  */
    ebi.assign( 4, 0. ); 
 r = 0; 

       
 while( r < trials ){
   r = r + 1;
   
    assert( n > 1. );
   tre.assign( (unsigned)n, 1.);
   ebi.assign( 4, 0. ); 
   Sjp[0] = 0. ;
   
  while( (unsigned)tre.size() > 1){ 
    /*
       for( auto& x: tre){ std::cout << ' ' << x; }
        std::cout << '\n';
        */
   /* sample time and $x$ */
 
    timi = itime( (double)tre.size(), a, M, b, x, Sjp);

  /* given sampled time update spectrum */
   if (theta > 0. ){
     rpois = std::poisson_distribution<unsigned>(theta * timi);}
      for( unsigned& y: tre){
      z = (theta > 0. ? (double)rpois(rng) : timi) ;
      ebi[(y < Kmin ? 3 :(y > Kmax ? 3 : 0))] = ebi[(y < Kmin ? 3 :(y > Kmax ? 3 : 0))] + z ;
      ebi[ (y < Imin ? 3 : (y > Imax ? 3 : 1))] = ebi[ (y < Imin ? 3 : (y > Imax ? 3 : 1))] + z ;
      ebi[ (y < Jmin ? 3 : (y > Jmax ? 3 : 2))] = ebi[ (y < Jmin ? 3 : (y > Jmax ? 3 : 2))] + z ; }

   if( (unsigned)tre.size() > 2){
    /* first sample  blocks to merge */
    /* v.resize( (unsigned)tre.size(), 0); */
    kblocks = samplek( (double)tre.size(),  x[0]);

    v.clear();
    v.assign( tre.begin(), tre.end()); 
    std::shuffle( v.begin(), v.end(), rng );
    tre.assign( v.begin(), v.end());
    
    size = 0 ; 
    for( it = tre.begin(); it != std::next( tre.begin(), kblocks); ++it){
    assert( *it > 0 );
    size = size + (*it);
    /* label the block pointed to by $it$ to be removed */
    *it = 0; }
      tre.push_back(size);
    tre.remove(0); }  
          else{
             tre.resize(1); }
    /* generated one tree */
  }  
/* print out average of spectrum */
       /* for(  double& x: ebi){ std::cout << x/ebi[0] << '\n' ; } */
       /*  std::cout << ebi[1]/ebi[0] << ' ' << ebi[2]/ebi[0] <<  '\n' ; */
  /*  
   FILE * f;
  f = fopen(resskra, "a");
  fprintf(f,  "%g %g\n", ebi[1]/ebi[0], ebi[2]/ebi[0]  );
  fclose(f); */
  ofstream f;
  f.open(resskra, ios::app );
  f <<  (ebi[0] > 0. ? ebi[1]/ebi[0] : 0.) << ' ' << (ebi[0] > 0. ? ebi[2]/ebi[0] : 0.)  << ' ' << ebi[0] << '\n' ;
  f.close();
  /*  if working with matrix `res` : 
res(r-1, 0) =  ebi[1]/ebi[0] ;
  res(r-1, 1) =  ebi[2]/ebi[0];
  res(r-1, 2) =  ebi[0];
  */ 
  /* close while loop over $trials$ */
  } 
 free(x); 
 free(Sjp);
 
/* close espectrum */
    }
