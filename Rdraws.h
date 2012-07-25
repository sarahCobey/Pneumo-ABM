/* 
*
* Rdraws.h
* Pneumo-ABM - S. Cobey
*
*/

#ifndef RDRAWS
#define RDRAWS

#include <boost/random.hpp>

double r01( boost::mt19937& );
void rmultinom( double[], int, int, int[], boost::mt19937& );
void rmultinom( const double[], const int, int, int[], boost::mt19937& );
int rbinom( int, double, boost::mt19937& );
int rmultinom( const double[], int, boost::mt19937& );
double rnorm( double, double, boost::mt19937& );
double rgamma( double, double, boost::mt19937& );
double rexp( double, boost::mt19937& );

#endif
