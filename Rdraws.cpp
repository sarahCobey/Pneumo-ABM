/*
* 
* Rdraws.cpp
* Pneumo-ABM - S. Cobey
*
*/

#include <cstdlib>
using namespace std;
#include <cmath>
#include <iostream>
#include <cassert>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/exponential_distribution.hpp>
#include <boost/random/uniform_real.hpp>

#include "Rdraws.h"

// Random double from [0,1)
double r01(  boost::mt19937& rng ) {
  boost::uniform_real<> uni_dist(0,1);
  boost::variate_generator<boost::mt19937&,boost::uniform_real<> > uni( rng, uni_dist );
  return uni();
}

// This function updates array of # individuals making each transition (for >2 transitions, i.e., multinomial)
void rmultinom( double p_trans[], int numTrials, int numTrans, int numEachTrans[], boost::mt19937& rng )
{
  double sumTrans = 0;
  for ( int i = 0; i < numTrans; i++) {
    sumTrans += p_trans[ i ];
  }
  int numLeaving = numTrials;

  // Implement multinomial through binom iterations
  numTrans -= 1;
  for ( int k = 0; k < numTrans; k++) { 
    double tmp;
    if ( p_trans[ k ] > sumTrans) {
      sumTrans = p_trans[ k ];
    }
    tmp = ( ( numLeaving > 0 ) && ( sumTrans > 0)) ? rbinom( numLeaving, p_trans[ k ] / sumTrans, rng ) : 0;
    numEachTrans[ k ] = (int)tmp;
    numLeaving -= numEachTrans[ k ];
    sumTrans -= p_trans[ k ];
  }
  numEachTrans[ numTrans ] = numLeaving;
}

// This function updates array of # individuals making each transition (for >2 transitions, i.e., multinomial)
void rmultinom( const double p_trans[], const int numTrials, int numTrans, int numEachTrans[], boost::mt19937& rng )
{
  double sumTrans = 0;
  for ( int i = 0; i < numTrans; i++) {
    sumTrans += p_trans[ i ];
  }
  int numLeaving = numTrials;

  // Implement multinomial through binom iterations
  numTrans -= 1;
  for ( int k = 0; k < numTrans; k++) {
    double tmp;
    if ( p_trans[ k ] > sumTrans) {
      sumTrans = p_trans[ k ];
    }
    tmp = ( ( numLeaving > 0 ) && ( sumTrans > 0)) ? rbinom( numLeaving, p_trans[ k ] / sumTrans, rng ) : 0;
    numEachTrans[ k ] = (int)tmp;
    numLeaving -= numEachTrans[ k ];
    sumTrans -= p_trans[ k ];
  }
  numEachTrans[ numTrans ] = numLeaving;
}

int rmultinom( const double p_trans[], int numTrans, boost::mt19937& rng )
{
  double sumTrans = 0;
  int ageAtDeath = numTrans-1;
  for ( int i = 0; i < numTrans; i++) {
    sumTrans += p_trans[ i ];
  }
  int numLeaving = 1;

  // Implement multinomial through binom iterations
  numTrans -= 1;
  for ( int k = 0; k < numTrans; k++) {
    if ( numLeaving == 1 ) {
      double tmp;
      if ( p_trans[ k ] > sumTrans) {
	sumTrans = p_trans[ k ];
      }
      tmp = ( ( numLeaving > 0 ) && ( sumTrans > 0)) ? rbinom( numLeaving, p_trans[ k ] / sumTrans, rng ) : 0;
      numLeaving -= (int)tmp;
      if ( numLeaving == 0 ) { 
	ageAtDeath = k; 
      }
      sumTrans -= p_trans[ k ];
    } // end while numLeaving == 1
  } // end for k in numTrans
  return ageAtDeath;
}

int rbinom( int nTrials, double pLeaving, boost::mt19937& rng )
{
  // THIS PRNG ONLY WORKS FOR SMALL TRIALS. If large trials needed, investigate TR1 or GSL PRNGS
  int successes = 0;
  for ( int n = 0; n < nTrials; n++ ) {
    if ( r01(rng) < pLeaving ) {
      successes++;
    }
  }
  return successes;
}

double rnorm( double mean, double sd, boost::mt19937& rng ) {
  boost::normal_distribution<> nd( mean, sd );
  boost::variate_generator<boost::mt19937&,boost::normal_distribution<> > var_nor( rng, nd );
  double val = -1.0;
  while ( val < 0 ) {
    val = var_nor(); 
  }
  return val;
}

double rgamma( double mean, double var, boost::mt19937& rng ) {
  // Note that boost's only gamma distribution requires scale = 1
  // See http://www.boost.org/doc/libs/1_43_0/boost/math/distributions/gamma.hpp
  // Convert by multiplying by scale after: see http://www.johndcook.com/cpp_TR1_random.html#gamma
  double shape = ( mean*mean )/var;
  double scale = var/mean;
  assert( shape > 0 );
  assert( scale > 0 );
  boost::gamma_distribution<> gd( shape );
  boost::variate_generator<boost::mt19937&,boost::gamma_distribution<> > var_gamma( rng, gd );
  return scale*var_gamma();
}

double rexp( double mean, boost::mt19937& rng ) {
  double lambda = 1.0/mean;
  boost::exponential_distribution<> ed( lambda );
  boost::variate_generator<boost::mt19937&,boost::exponential_distribution<> > var_exp( rng, ed );
  double val = var_exp();
  return val;
}
