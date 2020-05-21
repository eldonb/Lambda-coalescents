\pdfoutput=1
\documentclass[a4paper,10pt]{cweb}
\usepackage[T1]{fontenc}
\usepackage[utf8]{inputenc}%
\RequirePackage{fix-cm}%
\usepackage{amsfonts, amsmath, amssymb}
\usepackage{fullpage}%
\usepackage{marvosym}%
\usepackage{bm}
\usepackage{upgreek}
\usepackage[round,numbers,super]{natbib}
\usepackage[all]{xy}%
\usepackage{xcolor}%
%%\usepackage{sectsty}
\usepackage{rotating}%
\usepackage{a4wide,fullpage}%
\usepackage{setspace}%
\usepackage{enumerate}%
\usepackage{textcomp}
\usepackage{eufrak}%
\usepackage{hyperref}%
\usepackage{showkeys}
\usepackage{environ}
\usepackage{dsfont}%
\usepackage[right]{lineno}%
\usepackage{verbatim}
\usepackage{tabto}
\usepackage{environ}%
\usepackage{lipsum}%
\setstretch{1.5}%
\newcommand{\one}[1]{\ensuremath{\mathds{1}_{\left\{ #1 \right\}}}}
\newcommand{\g}{\,\boldsymbol{|}\,}%
\newcommand{\EE}[1]{\mathds{E}\left[ #1 \right]}%
\newcommand{\im}{\ensuremath{\imath} }%
\newcommand{\jm}{\ensuremath{\jmath} }%
\newcommand{\be}{\begin{equation}}%
\newcommand{\ee}{\end{equation}}%
\newcommand{\prb}[1]{\ensuremath{\mathds{P}\left( #1 \right) } }%
\newcommand{\h}[1]{\ensuremath{\uptheta_{ #1 } } }%
\newcommand{\VV}[1]{\ensuremath{ \mathbb{V}\left( #1 \right)}}%
\newcommand{\hp}{\ensuremath{\theta_1}}%
\newcommand{\hs}{\ensuremath{\theta_2}}%
\newcommand{\D}{\ensuremath{\mathbb{D}}}%
\newcommand{\F}{\ensuremath{\mathbb{F}} }%
\newcommand{\G}{\ensuremath{\mathbb{G}} }%
\newcommand{\bt}[1]{\textcolor{blue}{\tt #1}}%
\newcommand{\bi}[1]{\textcolor{blue}{\it #1}}%
\newcommand{\bb}[1]{\textcolor{blue}{\bf #1}}%
\newcommand{\Cpp}{{\texttt{C++ }}}
\NewEnviron{esplit}[1]{%
\begin{equation}
\label{#1}
\begin{split}
  \BODY
\end{split}\end{equation}
}
%%
\makeatletter
\renewcommand{\maketitle}{\bgroup\setlength{\parindent}{0pt}
\begin{flushright}
  \texttt{\textcolor{blue}{\Large \textit \@@title}}

 \textit{ \@@author}

  \@@date
\end{flushright}\egroup
}
\makeatother
%%
\title{ sampling $\Lambda$-Beta }%%
%\author{ Bjarki Eldon\\ leibniz institute at mfn berlin \\ 10115 Berlin,Germany\\ \href{mailto:eldon@@mfn-berlin.de}{\texttt{\Letter: my.name(\emph{you know the symbol})mfn-berlin.de}} }%
\author{ Bjarki Eldon \\  Leibniz Institute at MfN Berlin \\   \href{mailto:bjarki.eldon@@mfn-berlin.de}{\texttt{\Letter: my.name(\emph{you know the symbol})mfn-berlin.de}} }%
\date{\today }%
\begin{document}
\maketitle


\rule{\textwidth}{.8pt}


\begin{abstract}%
 The {\tt C++}  code embedded in this document samples branch lengths from the (incomplete) Beta$(2-\alpha,\alpha)$-coalescent.  This  CWEB
       \citep{knuth1994cweb} document (the {\tt .w} file) can be compiled
       with {\tt cweave} to generate a {\tt .tex} file, and with {\tt
       ctangle} to generate a {\tt .c} \citep{kernighan1988c} file, which can be compiled using a \Cpp compiler.
\end{abstract}%


\tableofcontents



@* \bb{Compilation}. 

You will need the Rcpp library (available from CRAN)  for calling C++ code for seamless integration of C++ and R, 
and the Boost C++ library

%% to be verified: should call it something else like wlambdacoal.w
first use {\tt ctangle} to compile a .cpp file: 
{\tt ctangle +s wclambdabeta.w wclambdabeta.cpp }

then from within R use:

{\tt 
> library(Rcpp) \\
> Rcpp::sourceCpp("wclambdabeta.cpp") \\
}


@* {Code}. 




@*1 {Includes}. 

the Include statements

@<includes>=@#
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


@*1 {Random number generators}. 

Use the standard random number generators for generating constructors

@<generators@>=@#
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


@*1 {Sampling a random Beta$(z;2-\alpha,\alpha)$ variate}. 

Sample a random variate from a Beta$(z; 2-\alpha, \alpha)$ distribution with 
density  
 \be
      f(x) =  x^{1 - \alpha}(1-x)^{\alpha - 1}\one{0<x \le z}dx
 \ee

We use rejection sampling

@<reject@>=@#
static double xreject( const double a, const double u, const double M)
{
/* need boost library */
 double x;

 randombeta =  boost::random::beta_distribution<double>(2.-a, a);
  x = randombeta(rng); 
  /* x =  boost::math::ibeta_inv(2.-a, a, u); */

Need to ensure $0 < x \le M$ where


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
