/*
*
* Host.cpp
* Pneumo-ABM - S. Cobey
*
*/

#include <cstdlib>
using namespace std;
#include <iostream>
#include <queue>
#include <vector>
#include <boost/unordered_map.hpp>
#include <boost/random.hpp>
#include <cmath>

#include "Event.h"
#include "Host.h"
#include "Parameters.h"
#include "Rdraws.h"
#include "Infection.h"
#include "SimPars.h"

#define HFLU_INDEX (INIT_NUM_STYPES-1)

// HOST CONSTRUCTOR & DESTRUCTOR
Host::Host( double t, double dob, int i, int h, int n, EventPQ & ce, SimPars * spPtr, boost::mt19937& rng ) {
  simParsPtr = spPtr;
  DOB = dob;
  id = i;
  household = h;
  neighborhood = n;
  calcLifeHist( t, ce, t - dob, rng );
  partner = 0;
  fledge = false;
  inf = false;
  for ( int s = 0; s < INIT_NUM_STYPES; s++ ){
    carriageSummary[ s ] = 0;
    immune[ s ] = 0;
    susc[ s ] = 1.0;
  } 
} 

Host::~Host() {
}

// SET FUNCTION DEFINITIONS
void Host::setHousehold( int h ) {
  household = h;
}

void Host::setNeighborhood( int n ) {
  neighborhood = n;
}

void Host::setGroup( int g ) {
  group = g;
}

void Host::setDOB( double d ) {
  DOB = d;
}

void Host::incrementAge( void ) {
  age++;
}

void Host::setPartner( int b ) {
  partner = b;
}
 
void Host::setFledge( bool f ) {
  fledge = f;
}

void Host::setInf( bool b ) {
  inf = b;
}

void Host::getVaccinated( void ) {
  for ( int s = 0; s < NUM_VACCINE_SEROTYPES; s++ ) {
    int thisSerotype = VACCINE_SEROTYPES[ s ];
    susc[thisSerotype] = min( susc[thisSerotype], 1.0-VACCINE_EFFICACY);
  }
}

// GET FUNCTION DEFINITIONS
int Host::getAgeInY( void ) const {
  return age;
}

int Host::getID( void ) const {
  return id;
}

int Host::getGroup( void ) const {
  return group;
}

int Host::getHousehold( void ) const {
  return household;
}

int Host::getNeighborhood( void ) const {
  return neighborhood;
}

bool Host::isPaired( void ) const {
  return ( partner > 0 );
}

bool Host::isAdult( void ) const {
  return ( (double)age >= MATURITY_AGE );
}

bool Host::isEligible( void ) const {
  return ( (double)age >= MATURITY_AGE && partner==0 ); 
}

bool Host::hasFledged( void ) const {
  return fledge;
}

int Host::isInfectedZ( int z ) const {
  return carriageSummary[ z ];
}

bool Host::isInfected( void ) const {
  return( isInfectedPneumo() || isInfectedHflu() );
}

bool Host::isInfectedPneumo( void ) const {
  return ( carriage.size() - carriageSummary[ HFLU_INDEX ] > 0.5 );
}

bool Host::isInfectedHflu( void ) const {
  return ( carriageSummary[ HFLU_INDEX ] > 0 );
}

int Host::totStrains( void ) const {
  int m = 0;
  for ( int s = 0; s < INIT_NUM_STYPES; s++ ) {
    m += (carriageSummary[s] > 0);
  }
  return carriage.size();
}

double Host::getSusc( int z ) const {
  return susc[ z ];
}

double Host::getDOD() const {
  return DOD;
}

double Host::getDOB() const {
  return DOB;
}

int Host::getSummedTheta() const {
  int numPastColonizations = 0;
  for ( int s = 0; s < HFLU_INDEX; s++ ) {
    if ( immune[s] > 0 ) {
      numPastColonizations++;	
    }
  }
  return numPastColonizations;
}


// MEMBER FUNCTION DEFINITIONS

void Host::calcLifeHist( double t, EventPQ &cePtr, double initAge, boost::mt19937 & rng ) {
  // Schedule death
  double ageAtDeath = 0.0;
  while ( ageAtDeath <= initAge ) {
    ageAtDeath = calcDeath( rng );
  }
  addEvent( t + ageAtDeath - initAge, DEATH_EVENT, id, cePtr );
  DOD = ageAtDeath - initAge + t;

  // Set current age (in years)
  age = floor( initAge/365.0); 

  // Make birthdays
  double thisT = DOB + 365.0;
  while ( thisT < t ) {
    thisT += 365.0;
  }
  while ( thisT < DOD ) {
    addEvent( thisT, BIRTHDAY, id, cePtr );
    thisT += 365.0;
  }

#ifdef SIM_PCV
  // Schedule vaccination -- implemented as change to susceptibility starting at VACCINATION_START
  if ( VACCINE_AGE < ageAtDeath && VACCINE_AGE > initAge && t >= ( DEM_SIM_LENGTH + VACCINATION_START - VACCINE_AGE ) ) { // means current n-month-olds will get vaccinated starting then
    addEvent( t + VACCINE_AGE - initAge, VACCINATION, id, cePtr );	
  }
#endif  

  // Leave home (if not already dead)
  double ageAtFledge = calcFledge( rng );
  if ( ageAtFledge < ageAtDeath && ageAtFledge > initAge ) {
    addEvent( t + ageAtFledge - initAge, FLEDGE_EVENT, id, cePtr );
  }
    
  // Potentially reproduce - no check for min time b/w births
  int numKids = calcNumBirths( rng );
  for ( int k = 0; k < numKids; k++ ) {
    double birthAge = calcBirthAge( rng );
    if ( birthAge < ageAtDeath && birthAge > initAge ) {
      addEvent( t + birthAge - initAge, BIRTH_EVENT, id, cePtr );
    }
  } 
    
  // Pair (unless dead) -- corrected for pairing with other initiators
  double perCapProb = ( 2.0 - sqrt( 4.0 - 4.0*PROB_PAIR ) )/2.0;
  if ( r01(rng) < perCapProb ) { 
    double pairAge = calcPairAge( rng );
    if ( pairAge < ageAtDeath && pairAge > initAge ) {
      addEvent( t + pairAge - initAge, PAIR_EVENT, id, cePtr );
    }
  }
}

double Host::calcDeath( boost::mt19937 & rng ) { 
  int age_death = rmultinom( simParsPtr->get_demPMF_row( LSPAN_INDEX ), INIT_NUM_AGE_CATS, rng );
  return 365.0*( (double)age_death + r01(rng) );
}

void Host::addEvent( double et, int eid, int hid, EventPQ & ceptr ) {
  Event thisEvent( et, eid, hid, 0 ); 
  ceptr.insert( thisEvent );
}

void Host::addEvent( double et, int eid, int hid, int s, EventPQ & ceptr ) {
  Event thisEvent( et, eid, hid, s );
  ceptr.insert( thisEvent );
}

double Host::calcFledge( boost::mt19937 & rng ) {
  int age_fledge = rmultinom( simParsPtr->get_demPMF_row( FLEDGE_INDEX ), INIT_NUM_AGE_CATS, rng );
  return 365.0*( (double)age_fledge + r01(rng) );
}

double Host::calcPairAge( boost::mt19937 & rng ) {
  int age_pair = rmultinom( simParsPtr->get_demPMF_row( PAIR_INDEX ), INIT_NUM_AGE_CATS, rng );
  return 365.0*( (double)age_pair + r01(rng) ); 
}

int Host::calcNumBirths( boost::mt19937 & rng ) {
  return rmultinom( simParsPtr->get_parity(), PARITY_BUFFER, rng );
}

double Host::calcBirthAge( boost::mt19937 & rng ) {
  int birth_age = rmultinom( simParsPtr->get_demPMF_row( BIRTH_AGE_INDEX ), INIT_NUM_AGE_CATS, rng );
  return 365.0*( (double)birth_age + r01(rng) );
}

void Host::becomeInfected( int s, double currentTime, EventPQ & ce, boost::mt19937& rng) {
    Infection thisInf( currentTime, 0.0 );
    carriage.insert( InfectionMap::value_type( s, thisInf ) ); 
    carriageSummary[ s ]++;
    if ( carriageSummary[ s ] == 1 ) { // imperfect (calls calcSusc() for unnecessary pneumo cases) but potentially most efficient
      calcSusc( currentTime );
    }
    double newRecTime, oldRecTime, infTime;
    for ( InfectionMap::iterator iter = carriage.begin(); iter != carriage.end(); iter++ ) {
      oldRecTime = (iter->second).recT;
      if ( oldRecTime == 0.0 ) {	
	if ( s < HFLU_INDEX && isInfectedHflu() && r01(rng) < simParsPtr->get_Hflu_prob(s) ) { // challenged with pneumo, has Hflu, and this pneumo strain can't handle it
	  newRecTime = currentTime + r01(rng)*EPSILON;
	} else { // challenging strain colonizes normally
	  infTime = (iter->second).infT;
	  newRecTime = currentTime + calcRecovery( s, currentTime, infTime, rng );
	} // end for pneumo and hflu colonization options
	while ( ce.count( Event(newRecTime) ) == 1 ) { 
	  cout << "\tHost id " << id << " is adjusting new recovery time to strain " << s << " to avoid collision with pre-existing events." << endl;
	  newRecTime += pow(10,APPROX_NOW);
	}
	(iter->second).recT = newRecTime;
	if ( newRecTime < DOD ) {
	  addEvent( newRecTime, RECOVERY_EVENT, id, s, ce );
	}
      } else if ( s == HFLU_INDEX && isInfectedPneumo() ) { // end if oldRecTime == 0.0 (i.e., dealing with recovery to challenging strain); now looking only for pneumo residents to Hflu challenge
	int thisS = iter->first; // figure out strain currently carried; do not want to update hflu strains
	if ( thisS < HFLU_INDEX && r01(rng) < simParsPtr->get_Hflu_prob(thisS) ) {
	  newRecTime = currentTime + r01(rng)*EPSILON;
	  oldRecTime = (iter->second).recT;
	  if ( oldRecTime != 0.0 && oldRecTime < DOD ) { 
	    EventPQ::iterator epqItr = ce.find( Event(oldRecTime) );
	    ce.erase( epqItr );
	  }
	  if ( newRecTime < DOD ) {
	    while ( ce.count( Event(newRecTime) ) == 1 ) { 
	      newRecTime += pow(10,APPROX_NOW);  
	    }
	    (iter->second).recT = newRecTime;
	    addEvent( newRecTime, RECOVERY_EVENT, id, thisS, ce );
	  } else { // ... else just update infection struct
	    (iter->second).recT = newRecTime;
	  }
	} // end if strain requires updating
      } // end if getting colonized with Hflu and already colonized with pneumo
    } // end for each strain carried
}

void Host::recover( int s, double recoverTime, EventPQ & ce ) {
  immune[ s ]++;
  carriageSummary[ s ]--;
  pair<InfectionMap::iterator,InfectionMap::iterator> ret = carriage.equal_range( s );
  InfectionMap::iterator it;
  for ( it = ret.first; it != ret.second; it++ ) {
    if ( ((*it).second).recT == recoverTime ) {
      carriage.erase( it );
      break;
    }
  }
  // Update susceptibilities if this is a new recovery or have lost current carriage
  if ( immune[ s ] == 1 || carriageSummary[ s ] == 0 ) { 
    calcSusc( recoverTime );
  }
}

double Host::calcRecovery( int s, double currentTime, double infectionTime, boost::mt19937& rng ) {
    int pastInfections = 0; // contribution of preexisting immunity (total number of conspecific carriage events cleared)
    double effectiveMean = 0;
    if ( s != HFLU_INDEX ) {
      for ( int z = 0; z < HFLU_INDEX; z++ ) {
	pastInfections += immune[z];
      }
    } else {
      pastInfections = immune[s];
    }
    effectiveMean = BASE_DURATION + ( simParsPtr->get_serotypePar_ij( MEAN_DURATION_INDEX, s ) - BASE_DURATION ) * exp( -REC_EPS * pastInfections );
    double rt = rexp( effectiveMean, rng );     // pulls duration from exponential distribution given effective mean
    return max( EPSILON*r01(rng), rt - ( currentTime - infectionTime ) );     // returns the amount of time left in infection
}

double Host::calcRecovery( int s, double currentTime, Infection & thisInf, boost::mt19937& rng ) {
  double infectionTime = thisInf.infT;
  return ( calcRecovery( s, currentTime, infectionTime, rng ) );
}

void Host::calcSusc( double t ) {

#ifdef SIM_PCV
  // For each pneumo serotype
  double runSum;
  double maxReduction;
  for ( int s = 0; s < HFLU_INDEX; s++ ) {
    runSum = 0.0;
    maxReduction = 0.0;
    for ( int z = 0; z < HFLU_INDEX; z++ ) {
      runSum += simParsPtr->get_XI_ij(s,z)*(double)( immune[z]>0 ); // contribution of past infections
      if ( isInfectedZ(z) > 0 ) {
	maxReduction = max( maxReduction, simParsPtr->get_reductions(z)); // contribution of current carriage
      }
    }
    if ( (t - DOB >= VACCINE_AGE) && (t >= DEM_SIM_LENGTH + VACCINATION_START) && ( IN_VACCINE[s] == 1 ) ) { // host vaccinated
      susc[ s ] = min( 1.0-VACCINE_EFFICACY, 1.0-min(1.0,runSum) );
    } else { // host not yet vaccinated
      susc[ s ] = 1.0 - min(1.0,runSum);
    }
    susc[ s ] *= ( 1.0 - maxReduction );
  } // end for each serotype
#else 
  // For each pneumo serotype
  double runSum;
  double maxReduction;
  for ( int s = 0; s < HFLU_INDEX; s++ ) {
    runSum = 0.0;
    maxReduction = 0.0;
    for ( int z = 0; z < HFLU_INDEX; z++ ) {
      runSum += simParsPtr->get_XI_ij(s,z)*(double)( immune[z]>0 ); // contribution of past infections
      if ( isInfectedZ(z) > 0 ) {
	maxReduction = max( maxReduction, simParsPtr->get_reductions(z));
      }
    }
    susc[ s ] = 1.0 - min(1.0,runSum);
    susc[ s ] *= ( 1.0 - maxReduction );
  } // end for each serotype
#endif

  // For Hflu
  susc[ HFLU_INDEX ] = 1.0 - min(1.0,simParsPtr->get_XI_ij(HFLU_INDEX,HFLU_INDEX)*(double)(immune[ HFLU_INDEX ]>0));
  if ( isInfectedHflu() ) {
    susc[ HFLU_INDEX ] *= ( 1.0 - RSCC_HFLU );
  }

}
