/*
*
* Simulation.cpp
* Pneumo-ABM - S. Cobey
*
*/

#include <cstdlib>
using namespace std;

#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <cassert>
#include <sstream>
#include <boost/random.hpp>

#include "Simulation.h"
#include "Parameters.h"
#include "Host.h"
#include "Event.h"
#include "Rdraws.h"
#include "Containers.h"

#define HFLU_INDEX (INIT_NUM_STYPES-1)

Simulation::Simulation( int trt, int sid, SimPars * spPtr ) {
  simID = sid;
  treatment = trt;
  simParsPtr = spPtr;
  t = 0;
  demOutputStrobe = 0;
  eventCtr = 0;
  idCtr = 1;
  hholdCtr = 1;
  numInfecteds[ INIT_NUM_AGE_CATS ][ INIT_NUM_STYPES ][ NUM_NEIGHBORHOODS ];
  demComplete = 0.0;
  rng.seed( (unsigned int)sid ); 
  for ( int a = 0; a < INIT_NUM_AGE_CATS; a++ ) {
    for ( int s = 0; s < INIT_NUM_STYPES; s++ ) {
      for ( int n = 0; n < NUM_NEIGHBORHOODS; n++ ) {
	numInfecteds[a][s][n] = 0;
      }
    }
  }

  // Set up multinomial to determine # of households in each size category
  // and # of people in each age class
  int init_age_dist[ INIT_NUM_AGE_CATS ];
  for ( int a = 0; a < INIT_NUM_AGE_CATS; a++ ) {
    init_age_dist[ a ] = 0;
  } 
  rmultinom( simParsPtr->get_demPMF_row( INIT_AGE_INDEX ), N0, INIT_NUM_AGE_CATS, init_age_dist, rng );
      
  // Create hosts
  double thisAge;
  double DOB;
  int randNeighborhood;
  for ( int a = 0; a < INIT_NUM_AGE_CATS; a++ ) {
    while ( init_age_dist[ a ] > 0 ) {
      thisAge = 365.0*( (double)a + r01(rng) );
      DOB = t - thisAge;
      randNeighborhood = (int)(floor( (double)NUM_NEIGHBORHOODS*r01(rng) ));
      Host * newHostPtr = new Host( t, DOB, idCtr, hholdCtr, randNeighborhood, currentEvents, simParsPtr, rng );
      allHosts.insert( boost::shared_ptr<Host>(newHostPtr) );
      allHouseholds.insert( hholdCtr );
      idCtr++;
      hholdCtr++;
      init_age_dist[ a ]--;
    } // end while (= this age category exhausted)
  } // all age categories done

  // Assign to households
  HostsByAge::iterator it = allHosts.get<age>().end();
  while ( it != allHosts.get<age>().begin() ) {
    it--;
    if ( (*it)->isEligible() ) { // If host is an adult and not yet paired
      double pairs = r01(rng);
      if ( pairs < PROB_PAIR/2.0 ) {
	int thisID = (*it)->getID();
	pairHost( thisID );
      }
    } else if ( (*it)->isAdult() == false ) { // If host is a kid, pick hosts randomly until find adult and join household
      bool a = false;
      int rID = 0;
      HostsByID::iterator ita;
      while ( a == false ) {
	rID = ceil( (double)N0*r01(rng) );
	ita = allHosts.find( rID );
	a = (*ita)->isAdult();
      }
      int newHH = (*ita)->getHousehold();
      int newNH = (*ita)->getNeighborhood();
      HostsByID::iterator itt = allHosts.find( (*it)->getID() );
      allHouseholds.erase( (*itt)->getHousehold() );
      allHosts.modify( itt, updateHousehold( newHH ) );
      allHosts.modify( itt, updateNeighborhood( newNH ) );
    }
  } // end while
  initOutput();
}

Simulation::~Simulation() {
  closeOutput();
}


// Stream functions
void Simulation::closeOutput() {
  ageDistStream.close();
  hhDistStream.close();
  demTimesStream.close();

  coinfectionHistStream.close();
  coinfectionHFHistStream.close();
  epidTimesStream.close();
  totCarriageStream.close();
}

void Simulation::initOutput() {
  for ( int n = 0; n < NUM_NEIGHBORHOODS; n++ ) {
    ageDistFile = makeBigName( "age_dist_neighborhood", n );
    ageDistStream.open( ageDistFile.c_str(),ios::out );
    ageDistStream.close();
  }
  hhDistFile = makeName( "hh_dist" ); 
  demTimesFile = makeName( "dem_times" );
  epidTimesFile = makeName( "epid_times" );
  coinfectionHistFile = makeName( "coinfection_dist_pneumo" );
  coinfectionHFHistFile = makeName( "coinfection_dist_hflu-pneumo" );
  totCarriageFile = makeName( "totCarriage" );

  hhDistStream.open( hhDistFile.c_str(),ios::out );
  demTimesStream.open( demTimesFile.c_str(),ios::out );
  epidTimesStream.open( epidTimesFile.c_str(), ios::out );
  coinfectionHistStream.open(coinfectionHistFile.c_str(),ios::out );
  coinfectionHFHistStream.open(coinfectionHFHistFile.c_str(),ios::out );
  totCarriageStream.open( totCarriageFile.c_str(),ios::out );

  for ( int s = 0; s < INIT_NUM_STYPES; s++ ) {
    for ( int n = 0; n < NUM_NEIGHBORHOODS; n++ ) {
      infectionStream.open( makeBiggerName( "infections", s, "neighborhood", n ).c_str(),ios::out);
      infectionStream.close();
      infectedStream.open( makeBiggerName( "infecteds", s, "neighborhood", n ).c_str(),ios::out );
      infectedStream.close();
    }
  }
}

void Simulation::writeDemOutput() {
  
  // I. Output age distribution to file
  HostsByAge::iterator it = allHosts.get<age>().end();
  it--;
  int maxAge = INIT_NUM_AGE_CATS-1;

  // For 0 to max age, count ages 
  HostsByAN& sorted_index = allHosts.get<an>();
  for ( int n = 0; n < NUM_NEIGHBORHOODS; n++ ) {
    ageDistFile = makeBigName( "age_dist_neighborhood", n );
    ageDistStream.open( ageDistFile.c_str(),ios::app);
    for ( int a = 0; a < maxAge+1; a++) {
      ageDistStream << allHosts.get<an>().count( boost::make_tuple(n, a )) << "\t";
    } 
    ageDistStream << endl;
    ageDistStream.close();
  }
  
  // Update dem times
  demTimesStream << demOutputStrobe << endl;

  // II. Output household distribution to file
  HostsByHH& sorted_index2 = allHosts.get<household>();
  int thisHH = 0;
  int hhSize = 0;
  int maxSize = 1;
  int households[ HHOLD_SIZE_BUFFER ]; // holds counts for sizes 1...HHOLD_SIZE_BUFFER+1
  for ( int s = 0; s < HHOLD_SIZE_BUFFER; s++ ) {
    households[ s ] = 0;
  }
  for ( HHSet::iterator hhit = allHouseholds.begin(); hhit != allHouseholds.end(); hhit++ ) {
    thisHH = *hhit;
    hhSize = sorted_index2.count( thisHH );
    if ( hhSize > maxSize ) {
      maxSize = hhSize;
    }
    if ( hhSize > floor(0.8*(double)HHOLD_SIZE_BUFFER ) ) {
      cout << "Set MAX_HHOLD_SIZE to be larger. Encountered household with " << hhSize << " members." << endl;
      assert( hhSize < HHOLD_SIZE_BUFFER );
    }
    households[ hhSize - 1 ]++;
  }
  for ( int s = 0; s < HHOLD_SIZE_BUFFER; s++ ) {
    hhDistStream << households[ s ] << "\t";
  }
  hhDistStream << endl;
}

void Simulation::writeTheta( void ) {
  thetaFile = makeName( "theta" );
  thetaStream.open( thetaFile.c_str(),ios::out);
  for ( HostsByAge::iterator it = allHosts.get<age>().begin(); it != allHosts.get<age>().end(); it++ ) { // For each host...
    if ( (*it)->getAgeInY() < 6 ) {
      thetaStream << t - (*it)->getDOB() << "\t" << (*it)->getSummedTheta() << endl;
    }
  }
  thetaStream.close();
}

void Simulation::writeEpidOutput( void ) {
  // I. Output infections by age to file
  for ( int n = 0; n < NUM_NEIGHBORHOODS; n++ ) { // Fix when adapt to more than one neighborhood
    for ( int s = 0; s < INIT_NUM_STYPES; s++ ) {
      infectionStream.open( makeBiggerName( "infections", s, "neighborhood", n ).c_str(),ios::app);
      for ( int a = 0; a < INIT_NUM_AGE_CATS; a++ ) {
	infectionStream << numInfecteds[ a ][ s ][ n ] << "\t";
      }
      infectionStream << endl;
      infectionStream.close();
    }
  }

  // II. Write time file
  epidTimesStream << epidOutputStrobe << endl;

  // III. Write infecteds by age to file for each serotype
  int actualInfecteds[ INIT_NUM_AGE_CATS ][ INIT_NUM_STYPES ][ NUM_NEIGHBORHOODS ];
  for ( int n = 0; n < INIT_NUM_STYPES; n++ ) {
    for ( int a = 0; a < INIT_NUM_AGE_CATS; a++ ) {
      for ( int s = 0; s < INIT_NUM_STYPES; s++ ) {
	actualInfecteds[ a ][ s ][ n ] = 0;
      }
    }
  }

  // IV. Write coinfections & total carriage to files
  int numCarryingPneumo = 0;
  int coinf[ HFLU_INDEX ];
  int coinfHflu[ HFLU_INDEX ];
  for ( int c = 0; c < HFLU_INDEX; c++ ) {
    coinf[ c ] = 0;
    coinfHflu[ c ] = 0;
  }
  int id, hhold, thisC;
  std::pair< HostsByInf::const_iterator, HostsByInf::const_iterator > pit = allHosts.get<inf>().equal_range(true);
  for ( HostsByInf::const_iterator fit = pit.first; fit != pit.second; fit++ ) {
    if ( (*fit)->isInfectedPneumo() ) {
      numCarryingPneumo++;
    }
    if ( (*fit)->isInfectedPneumo() && (*fit)->getAgeInY() < COCOL_AGE_LIMIT ) {
      bool hasHflu =  (*fit)->isInfectedHflu();
      thisC = -1;
      for ( int z = 0; z < HFLU_INDEX; z++ ) {
	if ( (*fit)->isInfectedZ(z) ) {
	  thisC++;
	  if ( hasHflu ) {
	    coinfHflu[ z ]++;
	  }
	}
      }
      coinf[ thisC ]++;
    }

    for ( int s = 0; s < INIT_NUM_STYPES; s++ ) {
      if ( (*fit)->isInfectedZ( s ) ) {
	actualInfecteds[(*fit)->getAgeInY()][s][(*fit)->getNeighborhood()]++;
      }
    } // end for each serotype
  } // end for all infected hosts
  for ( int c = 0; c < HFLU_INDEX; c++ ) {
    coinfectionHistStream << coinf[ c ] << "\t";
    coinfectionHFHistStream << coinfHflu[ c ] << "\t";
  }
  coinfectionHistStream << endl;
  coinfectionHFHistStream << endl;
  totCarriageStream << numCarryingPneumo << endl;

  for ( int n = 0; n < NUM_NEIGHBORHOODS; n++ ) {
    for ( int s = 0; s < INIT_NUM_STYPES; s++ ) { 
      infectedStream.open( makeBiggerName( "infecteds", s, "neighborhood", n ).c_str(),ios::app);
      for ( int a = 0; a < INIT_NUM_AGE_CATS; a++ ) {
	infectedStream << actualInfecteds[ a ][ s ][ n ] << "\t"; 
      }
      infectedStream << endl;
      infectedStream.close();
    }
  }
}


// Member function definitions

void Simulation::runDemSim( void ) {
  double percentDone = 0.0;
  cout << "  Starting demographic component with seed=" << simID << "." << endl;
  writeDemOutput();
  demOutputStrobe += STROBE_DEM;
  EventPQ::iterator eventIter = currentEvents.begin();
  while ( ( *eventIter ).time < DEM_SIM_LENGTH  && currentEvents.size() != 0 && allHosts.size() > 0 )   {
    Event thisEvent = *eventIter;
    t = thisEvent.time;
    while ( demOutputStrobe < t ) { 
      writeDemOutput();
      demOutputStrobe += STROBE_DEM;
    }
    while ( percentDone/100.0 <= t/DEM_SIM_LENGTH  ) {
      cout << "\t" << percentDone << "% of this component complete." << endl;
      percentDone += PROGRESS_INTERVAL;
    }
    executeEvent( thisEvent );
    eventCtr++;
    currentEvents.erase( eventIter );
    eventIter = currentEvents.begin();
  } // while time < DEM_SIM_LENGTH

  // Set time to end of simulation and reset strobe to calibrate with epid strobing
  t = DEM_SIM_LENGTH;
  demComplete = t;
  cout << "\t100% of this component complete." << endl;
}
 
double Simulation::runTestEpidSim( void ) {
  if ( allHosts.size() == 0 ) {  
    cerr << "No hosts remaining for epidemiological simulation. Cancelling." << endl;
    assert(false);
  }
  cout << "  Entering test simulation at t="  << demComplete << "." << endl;
  double percentDone = 0.0;

  // Initialize host population with infections
  demOutputStrobe = t;
  epidOutputStrobe = t;
  seedInfections();
  EventPQ::iterator eventIter = currentEvents.begin();
  double nextTimeStep = t + EPID_DELTA_T;
  double prevalences[ NUM_TEST_SAMPLES ]; // holds prevalences at strobing periods
  for ( int p = 0; p < NUM_TEST_SAMPLES; p++ ) {
    prevalences[ p ] = 0.0;
  }
  double prevYear = DEM_SIM_LENGTH + TEST_EPID_SIM_LENGTH - ( NUM_TEST_SAMPLES * 365.0 ); // first simulation time to start sampling
  int prevSamples = 0;

  while ( t < TEST_EPID_SIM_LENGTH + demComplete )
    {
      // Calculate new infections for every host and add events to stack
      calcSI();
      eventIter = currentEvents.begin();  
      while ( ( *eventIter ).time < nextTimeStep ) { 
	while ( demOutputStrobe < t ) { 
	  writeDemOutput();
	  demOutputStrobe += STROBE_DEM;
	}
	while ( epidOutputStrobe < t ) { 
	  writeEpidOutput();
	  epidOutputStrobe += STROBE_EPID;
	}
	if ( prevYear < t ) {
	  prevalences[ prevSamples ] = calcPrev();
	  cout << "\tOutputting prevalence sample #" << prevSamples+1 << "; prevalence of pneumo under 5 is " << prevalences[ prevSamples ] << endl;
	  prevSamples++;
	  prevYear += 365.0;
	}
      	while ( percentDone/100.0 <= ( t - demComplete )/TEST_EPID_SIM_LENGTH  ) {
	  cout << "\t" << percentDone << "% of this test component complete." << endl;
	  percentDone += PROGRESS_INTERVAL; 
	} 

	Event thisEvent = *eventIter;
	t = thisEvent.time;
	executeEvent( thisEvent );
	eventCtr++;
	currentEvents.erase( eventIter );
	eventIter = currentEvents.begin();
      }
      t = nextTimeStep;
      nextTimeStep += EPID_DELTA_T;
    }
  cout << "\t100% of this test component complete." << endl;
  double meanPrev = 0.0;
  double sumPrev = 0.0;
  int totSamples = 0;
  for (  int p = 0; p < NUM_TEST_SAMPLES; p++ ) {
    sumPrev += prevalences[ p ];
    totSamples++;
  }
  cout << "Counted " << totSamples << " total prevalence samples." << endl;
  meanPrev = sumPrev/(double)totSamples;
  return( meanPrev );
}

void Simulation::runEpidSim( void ) {
  if ( allHosts.size() == 0 ) {  
    cerr << "No hosts remaining for epidemiological simulation. Cancelling." << endl;
    assert(false);
  }
  cout << "  Entering epidemiological simulation at t="  << demComplete << "." << endl;
  double percentDone = 0.0;

  // Initialize host population with infections
  demOutputStrobe = t;
  epidOutputStrobe = t;
  seedInfections();
  EventPQ::iterator eventIter = currentEvents.begin();
  double nextTimeStep = t + EPID_DELTA_T;

  while ( t < EPID_SIM_LENGTH + demComplete )
    {
      // Calculate new infections for every host and add events to stack
      calcSI();
      eventIter = currentEvents.begin();      
      while ( ( *eventIter ).time < nextTimeStep ) { 
	while ( demOutputStrobe < t ) { 
	  writeDemOutput();
	  demOutputStrobe += STROBE_DEM;
	}
	while ( epidOutputStrobe < t ) { 
	  writeEpidOutput();
	  epidOutputStrobe += STROBE_EPID;
	}
      	while ( percentDone/100.0 < ( t - demComplete )/EPID_SIM_LENGTH  ) {
	  cout << "\t" << percentDone << "% of this component complete." << endl;
	  percentDone += PROGRESS_INTERVAL; 
	  } 
	Event thisEvent = *eventIter;
	t = thisEvent.time;
	executeEvent( thisEvent );
	eventCtr++;
	currentEvents.erase( eventIter );
	eventIter = currentEvents.begin();
      }
      t = nextTimeStep;
      nextTimeStep += EPID_DELTA_T;
    }
  cout << "\t100% of this component complete." << endl;
  writeTheta();
  cout << "  At end of simulation, " << allHosts.size() << " hosts and " << allHouseholds.size() << " households; cumulatively, " << idCtr-1 << " hosts and " << hholdCtr-1 << " households. " << eventCtr << " total events." << endl;
}


// PRIVATE SIMULATION FUNCTIONS
void Simulation::executeEvent( Event & te ) {
  switch ( te.eventID ) {
  case DEATH_EVENT :
    killHost( te.hostID );
    break;
  case FLEDGE_EVENT :
    fledgeHost( te.hostID );
    break; 
  case PAIR_EVENT :
    pairHost( te.hostID );
    break;
  case BIRTH_EVENT :
    birthHost( te.hostID );
    break; 
  case BIRTHDAY :
    ageHost( te.hostID );
    break;
  case INFECTION_EVENT :
    infectHost( te.hostID, te.s );
    break; 
  case RECOVERY_EVENT :
    recoverHost( te.hostID, te.s );
    break;
  case VACCINATION :
    vaccinateHost( te.hostID );
    break;
  default :
    cerr << "Event ID " << te.eventID << " attempted and failed." << endl;
    assert(false);
  } 
}

void Simulation::seedInfections( void ) {
  // Currently assuming: 
  // .... only single colonizations occur for each serotype
  // .... not considering co-colonizations when calculating duration of infection at this point
  HostsByAge::iterator it = allHosts.get<age>().begin();
  while ( it != allHosts.get<age>().end() ) {
    for ( int s = 0; s < INIT_NUM_STYPES; s++ ) {
      double thisRNG = r01(rng);
      if ( thisRNG < simParsPtr->get_serotypePar_ij( INIT_INFECTEDS_INDEX, s ) ) {
	int id = (*it)->getID();
	infectHost( id, s );
      }
    }
    it++;
  }
}

void Simulation::killHost( int id ) {
  HostsByID::iterator it = allHosts.find( id );
  // If only member of household, remove household from allHouseholds
  int hid = (*it)->getHousehold();
  int hsize = allHosts.get<household>().count( hid );
  if ( hsize == 1 ) {
    allHouseholds.erase( hid );
  } 

  // Update numInfecteds
  if ( demComplete != 0 ) {
    if ( (*it)->isInfected() ) {
      int age = (*it)->getAgeInY();
      int nhood = (*it)->getNeighborhood();
      for ( int s = 0; s < INIT_NUM_STYPES; s++) {
	numInfecteds[ age ][ s ][ nhood ] -= (*it)->isInfectedZ( s );
      }
    }
  }
  allHosts.erase( it );
}

void Simulation::ageHost( int id ) {
  HostsByID::iterator it = allHosts.find( id );
  if ( demComplete != 0 ) {
    if ( (*it)->isInfected() ) {
      int age = (*it)->getAgeInY();
      int nhood = (*it)->getNeighborhood();
      for ( int s = 0; s < INIT_NUM_STYPES; s++) {
	numInfecteds[ age ][ s ][ nhood ] -= (*it)->isInfectedZ( s );
	numInfecteds[ age+1 ][ s ][ nhood ] += (*it)->isInfectedZ( s );
      } // end for each serotype
    } // end for if host infected
  } // end for if in epid part of simulation
  allHosts.modify( it, updateAge() );
}

void Simulation::pairHost( int id ) {
  HostsByID::iterator it = allHosts.find( id );
  if ( (*it)->isPaired() == false ) {

    // Identify absorbing household (= initiating host's household )
    int hhold1 = (*it)->getHousehold();
  
    // Find eligible partner 
    int thisAge = (*it)->getAgeInY();
    int hhold2 = hhold1; // ensure do not pair w/in own household
    int count = 0;
    int partnerAge = 0;
    int attempts = 0;
    int rInd = 0;
    int id2  = 0;
    int famCount = 1;
    bool foundPartner = false;
    while ( count == 0 && attempts < DATE_BUFFER ) { // First attempt to find partner in preferred age range
      partnerAge = calcPartnerAge( thisAge );
      count = allHosts.get<aeh>().count(boost::make_tuple( partnerAge, true ));
      famCount = allHosts.get<aeh>().count(boost::make_tuple( partnerAge, true, hhold1 ));
      if ( count == 0 || famCount == count ) {
	attempts++;
      } else {
	while ( hhold1 == hhold2 && attempts < DATE_BUFFER ) {
	  rInd = ceil( (double)count * r01(rng) );
	  std::pair< HostsByAEH::const_iterator, HostsByAEH::const_iterator > pit = allHosts.get<aeh>().equal_range(boost::make_tuple( partnerAge, true ));
	  HostsByAEH::const_iterator partnerIt = pit.first;
	  for ( int i = 1; i < rInd; i++ ) {
	    partnerIt++;
	  }
	  hhold2 = (*partnerIt)->getHousehold();
	  if ( hhold1 != hhold2 ) { // Found partner; will update marital status
	    id2 = (*partnerIt)->getID();
	    foundPartner = true;
	  }
	  attempts++;
	} // end while hhold1==hhold2 && attempts < MAX_DATES
      } // end if (count == 0 ){} else{}
    } // end while count==0 && attempts < MAX_DATES

    // If no one found, search outside original range
    if ( foundPartner == false ) {
      count = allHosts.get<eh>().count( boost::make_tuple( true ));
      int famCount = allHosts.get<eh>().count( boost::make_tuple(true, hhold1 ));
      attempts = 0;
      hhold2 = hhold1;

      // Check if there are any not in household
      if ( count - famCount > 0 ) {
	while ( hhold1 == hhold2 ) {
	  rInd = ceil( (double)count * r01(rng) );
	  std::pair< HostsByEH::const_iterator, HostsByEH::const_iterator > pit = allHosts.get<eh>().equal_range( boost::make_tuple( true ));
	  HostsByEH::const_iterator partnerIt = pit.first;
	  for ( int i = 1; i < rInd; i++ ) {
	    partnerIt++;
	  }
	  hhold2 = (*partnerIt)->getHousehold();
	  if ( hhold1 != hhold2 ) { // Found partner; will update marital status
	    id2 = (*partnerIt)->getID();
	    foundPartner = true;
	  }
	} // end while (hhold1==hhold2)
      } // end if count > 0
    } // end if count == 0

    if ( foundPartner == true ) {
      partner2Hosts( id, id2 );
    } // end if partner found (count > 0)
  } // end if host eligible
}

void Simulation::partner2Hosts( int hid1, int hid2 ) {
      // Update both partners' partner IDs and get household IDs
      HostsByID::iterator it1 = allHosts.find( hid1 ); // Original host
      int hhold1 = (*it1)->getHousehold();
      int nhood1 = (*it1)->getNeighborhood();
      allHosts.modify( it1, updatePartner( hid2 ) ); // Partner host
      HostsByID::iterator it2 = allHosts.find( hid2 ); // Partner host
      int hhold2 = (*it2)->getHousehold();
      int nhood2 = (*it2)->getNeighborhood();
      allHosts.modify( it2, updatePartner( hid1 ) ); 

      // Check if initiating host still lives at home. If so, fledge and give new household & neighborhood.
      if ( (*it1)->hasFledged() == false ) {
	allHosts.modify( it1, updateFledge( true ) ); 
	int famSize = allHosts.get<household>().count( hhold1 );
	if ( famSize != 1 ) {
	  allHosts.modify( it1, updateHousehold( hholdCtr ) );
	  int randNeighborhood = (int)( floor( (double)NUM_NEIGHBORHOODS*r01(rng) ) );
	  allHosts.modify( it1, updateNeighborhood( randNeighborhood ) );
	  if ( (*it1)->isInfected() ) {
	    for ( int s = 0; s < INIT_NUM_STYPES; s++ ) {
	      int infections = (*it1)->isInfectedZ(s);
	      if ( infections > 0 ) {
		numInfecteds[ (*it1)->getAgeInY() ][ s ][ nhood1 ] -= infections;
		numInfecteds[ (*it1)->getAgeInY() ][ s ][ randNeighborhood ] += infections; 
	      }
	    }
	  }
	  allHouseholds.insert( hholdCtr );
	  hholdCtr++;
	  hhold1 = (*it1)->getHousehold();
	  nhood1 = (*it1)->getNeighborhood();
	}
      }

      // Update household of partner and partner's kids, if any
      int numFamily = allHosts.get<household>().count( hhold2 );
      bool fledged = (*it2)->hasFledged();

      // Remove household2 from master list of households, IF host lives alone or has fledged
      if ( numFamily == 1 ) {
	allHouseholds.erase( hhold2 );
      }

      if ( numFamily > 1 && fledged == true ) { // New partner living independently with own family
	// Get IDs of all people whose households need updating
	int familyIDs[ numFamily ];
	for ( int f = 0; f < numFamily; f++ ) {
	  familyIDs[ f ] = 0;
	}
	int indID = 0;
	int f = 0;
	std::pair< HostsByHH::iterator, HostsByHH::iterator > pit = allHosts.get<household>().equal_range( hhold2 );
	for ( HostsByHH::const_iterator fit = pit.first; fit != pit.second; fit++ ) {
	  indID = (*fit )->getID();
	  familyIDs[ f ] = indID;
	  f++;
	} 
	// Now update households and neighborhoods
	HostsByID::iterator it;
	for ( int f = 0; f < numFamily; f++ ) {
	  indID = familyIDs[f];
	  it = allHosts.find( indID ); 
	  int thisAge = (*it)->getAgeInY();
	  int origNhood = (*it)->getNeighborhood();
	  allHosts.modify( it, updateHousehold( hhold1 ) );
	  allHosts.modify( it, updateNeighborhood( nhood1 ) );
	  if ( (*it)->isInfected() ) {
	    for ( int s = 0; s < INIT_NUM_STYPES; s++ ) {
	      int infections = (*it)->isInfectedZ(s);
	      numInfecteds[ thisAge ][ s ][ origNhood ] -= infections;
	      numInfecteds[ thisAge ][ s ][ nhood1 ] += infections;
	    }
	  }
	}
	allHouseholds.erase( hhold2 ); 
      } else {  // New partner still living in family of origin, so just update partner
	int origNhood = (*it2)->getNeighborhood();
	int age = (*it2)->getAgeInY();
	allHosts.modify( it2, updateHousehold( hhold1 ) );
	allHosts.modify( it2, updateNeighborhood( nhood1 ) );
	if ( (*it2)->isInfected() ) {
	  for ( int s = 0; s < INIT_NUM_STYPES; s++ ) {
	    int infections = (*it2)->isInfectedZ(s);
	    numInfecteds[ age ][ s ][ origNhood ] -= infections;
	    numInfecteds[ age ][ s ][ nhood1 ] += infections;
	  }
	}
	allHosts.modify( it2, updateFledge( true ) ); 
      }
} 

int Simulation::calcPartnerAge( int a ) {
  double partnerAge = (double)MATURITY_AGE - 1.0;
  double rand_offset = 0.0;
  while ( partnerAge < MATURITY_AGE ) {
    rand_offset = r01(rng) * 2.0 * STD_AGE_PAIR - STD_AGE_PAIR;
    partnerAge = a + rand_offset;
  }
  return floor(partnerAge);
}

void Simulation::fledgeHost( int id ) {
  HostsByID::iterator it = allHosts.find( id );
  if ( (*it)->isPaired() == false ) {
    // First consider case of orphan living alone - no new household needed
    int currentHH = (*it)->getHousehold();
    int currentNhood = (*it)->getNeighborhood();
    int hhSize = allHosts.get<household>().count( currentHH );
    if ( hhSize == 1 ) {
      allHosts.modify( it, updateFledge( true ) );
    } else {
      allHosts.modify( it, updateHousehold( hholdCtr ) );
      int randNeighborhood = (int)(floor( (double)NUM_NEIGHBORHOODS*r01(rng) ));
      allHosts.modify( it, updateNeighborhood( randNeighborhood )); 
      int age = (*it)->getAgeInY();
      int infections;
      if ( (*it)->isInfected() ) {
	for ( int s = 0; s < INIT_NUM_STYPES; s++ ) {
	  infections = (*it)->isInfectedZ(s);
	  if ( infections > 0 ) {
	    numInfecteds[age][s][ currentNhood ] -= infections;
	    numInfecteds[age][s][ randNeighborhood ] += infections;
	  }
	}
      }
      allHosts.modify( it, updateFledge( true ) );
      allHouseholds.insert( hholdCtr ); 
      hholdCtr++;
    } // end if come from household > 1
  } 
}

void Simulation::birthHost( int id ) {
  HostsByID::iterator it = allHosts.find( id );
  int hhold = (*it)->getHousehold();
  int nhood = (*it)->getNeighborhood();	
  double DOB = t;
  Host * newHostPtr;
  newHostPtr = new Host( t, DOB, idCtr, hhold, nhood, currentEvents, simParsPtr, rng ); 
  allHosts.insert( boost::shared_ptr<Host>(newHostPtr) );
  idCtr++;   
}

void Simulation::infectHost( int id, int s ) { // Could rewrite to take itr instead of id
  HostsByID::iterator it = allHosts.find( id );
  if ( s < HFLU_INDEX || (*it)->isInfectedHflu()==false ) {
    numInfecteds[ (*it)->getAgeInY() ][ s ][ (*it)->getNeighborhood() ]++;
  }
  if ( s < HFLU_INDEX || (*it)->isInfectedHflu() == false ) {
    allHosts.modify( it, updateCarriage( s, t, currentEvents, rng ) );
    allHosts.modify( it, updateInf( true ) );
  }
}

void Simulation::recoverHost( int id, int s ) {
  HostsByID::iterator it = allHosts.find( id );

  // Host clears just one infection. All infections represented in numInfecteds[][].
  if ( s != HFLU_INDEX || ( s == HFLU_INDEX && (*it)->isInfectedZ(s) == 1 ) ) {
    numInfecteds[ (*it)->getAgeInY() ][ s ][ (*it)->getNeighborhood() ]--;
  }
  allHosts.modify( it, updateRecovery( s, t, currentEvents ) );
  if ( (*it)->totStrains() == 0 ) {
    allHosts.modify( it, updateInf( false ) );
  }
}

void Simulation::vaccinateHost( int id ) {
  HostsByID::iterator it = allHosts.find( id );
  allHosts.modify( it, updateVaccination() );	
}

double Simulation::calcPrev() {
  // Get population sizes of kids <5
  int N_total = 0;
  int I_total_pneumo = 0;
  int I_total_hflu = 0;
  HostsByAge& sorted_index = allHosts.get<age>();
  for ( int a = 0; a < 5; a++ ) {
    N_total += sorted_index.count( a );
  }

  // Now count how many kids in each age group are infected -- note that not sufficient to use numInfecteds, which counts 'effective' number of infecteds 
  // (i.e., the number of infections) for purposes of calculating the force of infection
  int hostAge = 0;
  for ( HostsByAge::iterator it = allHosts.get<age>().begin(); it != allHosts.get<age>().end(); it++ ) { // For each host...
    hostAge = (*it)->getAgeInY();
    if ( hostAge < 5 ) {
      I_total_pneumo += (*it)->isInfectedPneumo();
      I_total_hflu += (*it)->isInfectedHflu();
    }
  } // end for each host

  cout << "\tThere are " << N_total << " kids <5 y old; " << I_total_pneumo << " (" << 100.0*(double)I_total_pneumo/(double)N_total << "%) carry pneumo and " 
       << I_total_hflu << " (" << 100.0*(double)I_total_hflu/(double)N_total << "%) carry Hflu" << endl;
  return (  (double)I_total_pneumo/(double)N_total );
}


void Simulation::calcSI() {

#if defined( NO_HHOLDS ) && defined( NO_AGE_ASSORT ) // RANDOM MIXING: NO HOUSEHOLDS AND NO AGE-ASSORTATIVITY
  int N_total = 0;
  int nhood = 0;
  HostsByAge& sorted_index = allHosts.get<age>();
  for ( int n = 0; n < INIT_NUM_AGE_CATS; n++ ) {
    N_total += sorted_index.count( n );
  }
  int I_total[ INIT_NUM_STYPES ];
  for ( int s = 0; s < INIT_NUM_STYPES; s++ ) {
    I_total[s] = 0;
  }
  for ( int s = 0; s < INIT_NUM_STYPES; s++ ) {
    for ( int a = 0; a < INIT_NUM_AGE_CATS; a++ ) {
      I_total[ s ] += numInfecteds[a][s][nhood]; 
    }
  }
  double prInf;
  double rTrans;
  double infectionTime;
  for ( HostsByAge::iterator it = allHosts.get<age>().begin(); it !=  allHosts.get<age>().end(); it++ ) { 
    for ( int s = 0; s < INIT_NUM_STYPES; s++ ) {
      rTrans = simParsPtr->get_serotypePar_ij( BETA_INDEX, s ) * (double)( I_total[s] - (*it)->isInfectedZ(s) )/(double)( N_total - 1 ); // beta*I/N (excluding self)
      prInf = (*it)->getSusc(s) * ( rTrans + simParsPtr->get_serotypePar_ij( IMMIGRATION_INDEX, s ) );
      prInf = 1.0 - exp( -prInf * EPID_DELTA_T ); 
      if ( r01(rng) < prInf ) {
	infectionTime = (double)r01(rng) * (double)EPID_DELTA_T + (double)t;
	if ( infectionTime <= t ) {
		cout << "\tAdding epsilon to event time." << endl;
		infectionTime += pow(10,APPROX_NOW); 
	}
	if ( infectionTime < (*it)->getDOD() ) {
	  addEvent( infectionTime, INFECTION_EVENT, (*it)->getID(), s );
	}
      }
    } // end for each serotype
  } // end for each host
#endif 

#if defined( NO_HHOLDS ) && !defined( NO_AGE_ASSORT ) // AGE-ASSORTATIVE MIXING BUT NO HOUSEHOLDS
  // Get population ages and normalize age-assortative matrix, alpha
  int N_total = 0;
  int N_age[ INIT_NUM_AGE_CATS ];
  double alpha[ INIT_NUM_AGE_CATS ][ INIT_NUM_AGE_CATS ];
  HostsByAge& sorted_index = allHosts.get<age>();
  for ( int n = 0; n < INIT_NUM_AGE_CATS; n++ ) {
    N_age[ n ] = sorted_index.count( n );
    N_total += N_age[ n ];
    for ( int n2 = 0; n2 < INIT_NUM_AGE_CATS; n2++ ) {
      alpha[n][n2] = 0.0;
    }
  }
  double runningSum;
  for ( int i = 0; i < INIT_NUM_AGE_CATS; i++ ) { // for donors of age j
    if ( N_age[i] > 0 ) {
      runningSum = 0.0;
      for ( int j = 0; j < INIT_NUM_AGE_CATS; j++ ) {
	if (( N_age[j] > 1 ) || ((N_age[j]==1) && (i!=j)) ) { // only consider transmission b/w two age groups if both present and not just self
	  runningSum += simParsPtr->get_waifw_ij(i,j);
	} 
      }
      if ( runningSum > 0 ) {
	for ( int j = 0; j < INIT_NUM_AGE_CATS; j++ ) {
	  if ( ( N_age[j] > 1 ) || ( (N_age[j]==1) && (i!=j) ) ) {
	    alpha[i][j] = simParsPtr->get_waifw_ij(i,j) / runningSum;
	  }
	}
      }
    }
  }

  // Calculate force for every host
  int hostID;
  int hostAge;
  int numInfectionsZ;
  int nhood;
  double susc_z;
  double prInf;
  double rTrans;
  double infectionTime;
  for ( HostsByAge::iterator it = allHosts.get<age>().begin(); it !=  allHosts.get<age>().end(); it++ ) { // For each host...
    hostID = (*it)->getID();
    hostAge = (*it)->getAgeInY();
    nhood = (*it)->getNeighborhood();
    for ( int s = 0; s < INIT_NUM_STYPES; s++ ) {
      susc_z = (*it)->getSusc(s);
      numInfectionsZ = (*it)->isInfectedZ(s);
      rTrans = 0.0;
      for ( int a = 0; a < INIT_NUM_AGE_CATS; a++ ) {
	if ( a != hostAge && N_age[a] > 0 ) {
	  rTrans += simParsPtr->get_serotypePar_ij( BETA_INDEX, s ) * (double)( numInfecteds[a][s][nhood] )/(double)N_age[a] * alpha[hostAge][a]; 
	} else if ( N_age[a] > 1 ) { // no force from self if only member of cohort
	  rTrans += simParsPtr->get_serotypePar_ij( BETA_INDEX, s ) * (double)( numInfecteds[a][s][nhood] - numInfectionsZ )/(double)( N_age[a] - 1 ) * alpha[hostAge][a]; 
	}
      } // end for each age	
      assert( rTrans >= 0.0 );
      prInf = susc_z * rTrans + simParsPtr->get_serotypePar_ij( IMMIGRATION_INDEX, s ); 
      prInf = 1.0 - exp( -prInf * EPID_DELTA_T ); 
      if ( r01(rng) < prInf ) {
	infectionTime = r01(rng) * EPID_DELTA_T + t;
	if ( infectionTime <= t ) {
		cout << "\tAdding epsilon to event time." << endl;
		infectionTime += pow(10,APPROX_NOW); 
	}
	if ( infectionTime < (*it)->getDOD() ) {
	  addEvent( infectionTime, INFECTION_EVENT, hostID, s );
	}
      }
    } // end for each serotype
  } // end for each host
#endif


#if !defined( NO_HHOLDS ) && defined( NO_AGE_ASSORT ) // HOUSEHOLDS BUT RANDOM MIXING
  // Count total size of each neighborhood
  int neighborhoodSizes[ NUM_NEIGHBORHOODS ];
  int thisSize;
  for ( int n = 0; n < NUM_NEIGHBORHOODS; n++ ) {
    thisSize = 0;
    for ( int a = 0; a < INIT_NUM_AGE_CATS; a++ ) {
      thisSize += allHosts.get<an>().count( boost::make_tuple( n, a ) );
    }
    neighborhoodSizes[ n ] = thisSize;
  } // end for each neighborhood

  // Calculate force for every host
  int hostID;
  int otherID;
  int otherHhold;
  int hostHhold;
  int hostNhood;
  int hholdSize;
  int numInfectionsZ;
  int numHholdInfecteds;
  double numNhoodInfecteds;
  double susc_z;
  double prInf;
  double rTrans;
  double infectionTime;
  double n_weight;
  for ( HostsByAge::iterator it = allHosts.get<age>().begin(); it !=  allHosts.get<age>().end(); it++ ) {
    hostID = (*it)->getID();
    hostHhold = (*it)->getHousehold();
    hostNhood = (*it)->getNeighborhood();
    hholdSize = allHosts.get<household>().count( hostHhold ) - 1; // do not count self in household 

    for ( int s = 0; s < INIT_NUM_STYPES; s++ ) {
      susc_z = (*it)->getSusc(s);
      numInfectionsZ = (*it)->isInfectedZ(s);
      numHholdInfecteds = 0;
      rTrans = 0.0;

      if ( hholdSize > 0 ) {  
	std::pair< HostsByHH::iterator, HostsByHH::iterator > pit = allHosts.get<household>().equal_range( hostHhold ); 
	for ( HostsByHH::iterator fit = pit.first; fit != pit.second; fit++ ) {
	  numHholdInfecteds += (*fit)->isInfectedZ(s); // adding self in numerator (easier than checking family members' IDs)
	} // end for each host in household
 	rTrans += RHO_H * simParsPtr->get_serotypePar_ij( BETA_INDEX, s ) * (double)( numHholdInfecteds - numInfectionsZ )/(double)hholdSize; // not counting self in denominator
      } // end for household size > 0

      // Do non-household contacts, weighted for each contacted neighborhood
      for ( int n = 0; n < NUM_NEIGHBORHOODS; n++ ) {
	n_weight = simParsPtr->get_normalized_neighbor( hostNhood, n );
	if ( n_weight > 0 ) {
	  numNhoodInfecteds = 0; 
	  for ( int a = 0; a < INIT_NUM_AGE_CATS; a++ ) {
	    numNhoodInfecteds += (double)numInfecteds[ a ][ s ][ n ];
	  } // end for this age
	  if ( n != hostNhood ) {
	    rTrans += ( 1.0 - RHO_H ) * simParsPtr->get_serotypePar_ij( BETA_INDEX, s ) * (double)n_weight*( numNhoodInfecteds )/(double)( neighborhoodSizes[ n ] ); 
	  } else {
	    rTrans += ( 1.0 - RHO_H ) * simParsPtr->get_serotypePar_ij( BETA_INDEX, s ) * (double)n_weight*( numNhoodInfecteds - numHholdInfecteds )/(double)( neighborhoodSizes[ n ] - 1 - hholdSize );
	  }
	} // end for if n_weight > 0 
      } // end for each neighborhood
	
      prInf = susc_z * rTrans + simParsPtr->get_serotypePar_ij( IMMIGRATION_INDEX, s ); 
      prInf = 1.0 - exp( -prInf * EPID_DELTA_T ); 
      if ( r01(rng) < prInf ) {
	infectionTime = r01(rng) * EPID_DELTA_T + t;
	if ( infectionTime <= t ) {
	  cout << "\tAdding epsilon to event time." << endl;
	  infectionTime += pow(10,APPROX_NOW); 
	}
	if ( infectionTime < (*it)->getDOD() ) {
	  addEvent( infectionTime, INFECTION_EVENT, hostID, s );
	}
      } // end for if r01(rng)<prInf
    } // end for each serotype
  } // end for each host
#endif


#if !defined( NO_HHOLDS ) && !defined( NO_AGE_ASSORT ) // HOUSEHOLDS AND AGE-ASSORTATIVE MIXING
  int N_total[ NUM_NEIGHBORHOODS ];
  for ( int n = 0; n < NUM_NEIGHBORHOODS; n++ ) {
    N_total[ n ] = 0;
  }
  int N_age[ INIT_NUM_AGE_CATS ][ NUM_NEIGHBORHOODS ];
  double alpha[ INIT_NUM_AGE_CATS ][ INIT_NUM_AGE_CATS ][ NUM_NEIGHBORHOODS ];
  HostsByAge& sorted_index = allHosts.get<age>();
  for ( int h = 0; h < NUM_NEIGHBORHOODS; h++ ) {
    for ( int n = 0; n < INIT_NUM_AGE_CATS; n++ ) {
      N_age[ n ][ h ] = allHosts.get<an>().count( boost::make_tuple( h, n ) );
      N_total[ h ] += N_age[ n ][ h ];
      for ( int n2 = 0; n2 < INIT_NUM_AGE_CATS; n2++ ) {
	alpha[n][n2][h] = 0.0;
      }
    } // end for age n
  } // end for neighborhood h

  // Calculate alpha_nh for each neighborhood
  double runningSum;
  for ( int n = 0; n < NUM_NEIGHBORHOODS; n++ ) {
    for ( int i = 0; i < INIT_NUM_AGE_CATS; i++ ) { // for donors of age j
      if ( N_age[i][n] > 0 ) {
	runningSum = 0.0;
	for ( int j = 0; j < INIT_NUM_AGE_CATS; j++ ) {
	  if (( N_age[j][n] > 1 ) || ( (N_age[j][n]==1) && (i!=j) )) { // only consider transmission b/w two age groups if both present
	    runningSum += simParsPtr->get_waifw_ij(i,j);
	  }
	}
	if ( runningSum > 0 ) {
	  for ( int j = 0; j < INIT_NUM_AGE_CATS; j++ ) {
	    if (( N_age[j][n] > 1 ) || ( (N_age[j][n]==1) && (i!=j) )) {
	      alpha[i][j][n] = simParsPtr->get_waifw_ij(i,j) / runningSum;
	    }
	  }
	}
      }
    } // end for donors of age j
  } // end for this neighborhood n

  // Calculate force for every host
  int hostID;
  int otherID;
  int hostHhold;
  int otherHhold;
  int hostNhood;
  int hholdSize;
  int hostAge;
  int numInfectionsZ;
  int numHholdInfecteds;
  double numNhoodInfecteds;
  double n_weight;
  double alpha_hh[ INIT_NUM_AGE_CATS ]; // alpha_ij for this host i
  double susc_z;
  double prInf;
  double rTrans;
  double infectionTime;
  double sumAlphas;
  int hholdAges[ INIT_NUM_AGE_CATS ];
  for ( HostsByAge::iterator it = allHosts.get<age>().begin(); it !=  allHosts.get<age>().end(); it++ ) { // For each host...
    hostID = (*it)->getID();
    hostHhold = (*it)->getHousehold();
    hostNhood = (*it)->getNeighborhood();
    hostAge = (*it)->getAgeInY();
    hholdSize = allHosts.get<household>().count( hostHhold ) - 1; // do not count self in household
    sumAlphas = 0.0;
    for ( int a = 0; a < INIT_NUM_AGE_CATS; a++ ) {
       hholdAges[ a ] = 0;
       alpha_hh[ a ] = 0;
    }

    // Update alpha_hh vector for household members
    if ( hholdSize > 0 ) {
      for ( int a = 0; a < INIT_NUM_AGE_CATS; a++ ) { // ...for each age 
	std::pair< HostsByAH::iterator, HostsByAH::iterator > pit = allHosts.get<ah>().equal_range( boost::make_tuple( hostHhold, a ) );
	for ( HostsByAH::iterator fit = pit.first; fit != pit.second; fit++ ) {
	  if ( (*fit )->getID() != hostID ) {
	    hholdAges[ a ]++;
	  }
	} // end for each household member of this age
	if ( hholdAges[ a ] > 0 ) { // if there's a non-self household member of this age
	  sumAlphas += simParsPtr->get_waifw_ij(hostAge,a);
	}
      } // end for each age
      if ( sumAlphas > 0 ) {
	for ( int a = 0; a < INIT_NUM_AGE_CATS; a++ ) {
	  if ( hholdAges[ a ] > 0  ) {
	    alpha_hh[ a ] = simParsPtr->get_waifw_ij(hostAge,a) / sumAlphas;
	  }
	}
      }
    } // end for if hholdSize > 0

    for ( int s = 0; s < INIT_NUM_STYPES; s++ ) { // for each serotype
      susc_z = (*it)->getSusc(s);
      rTrans = 0.0;
      numInfectionsZ = (*it)->isInfectedZ(s);
      for ( int a = 0; a < INIT_NUM_AGE_CATS; a++ ) {
	numHholdInfecteds = 0;
	if ( hholdSize > 0 ) {
	  std::pair< HostsByAH::iterator, HostsByAH::iterator > pit = allHosts.get<ah>().equal_range( boost::make_tuple( hostHhold, a ) );
	  for ( HostsByAH::iterator fit = pit.first; fit != pit.second; fit++ ) {
	    numHholdInfecteds += (*fit)->isInfectedZ(s); // note that this number includes the host
	  } // end for all hosts in household of this age
	  if ( hholdAges[ a ] > 0 && a == hostAge ) { // Calculate household transmission for this age
	    rTrans += RHO_H * simParsPtr->get_serotypePar_ij( BETA_INDEX, s ) * (double)( numHholdInfecteds - numInfectionsZ )/(double)hholdAges[a] * alpha_hh[a];
	  } else if ( hholdAges[ a ] > 0 && a != hostAge ) {
	    rTrans += RHO_H * simParsPtr->get_serotypePar_ij( BETA_INDEX, s ) * (double)( numHholdInfecteds )/(double)hholdAges[a] * alpha_hh[a];
	  }
	} // end if hholdSize > 0

	for ( int n = 0; n < NUM_NEIGHBORHOODS; n++ ) {
	  n_weight = simParsPtr->get_normalized_neighbor( hostNhood, n );
	  if ( n_weight > 0 ) {
	    numNhoodInfecteds = numInfecteds[ a ][ s ][ n ];
	    if ( n == hostNhood && a == hostAge && ( N_age[a][n] > hholdAges[a]+1 ) ) { // same neighborhood and same age as host
	      rTrans += ( 1.0 - RHO_H ) * simParsPtr->get_serotypePar_ij( BETA_INDEX, s ) * n_weight*(double)( numNhoodInfecteds - numHholdInfecteds )/(double)( N_age[a][n] - hholdAges[a] - 1 ) * alpha[hostAge][a][n];
	    } else if ( n == hostNhood && a != hostAge && ( N_age[a][n] > hholdAges[a] ) ) { // same neighborhood but different age
	      rTrans += ( 1.0 - RHO_H ) * simParsPtr->get_serotypePar_ij( BETA_INDEX, s ) * n_weight*(double)( numNhoodInfecteds - numHholdInfecteds )/(double)( N_age[a][n] - hholdAges[a] ) * alpha[hostAge][a][n];
	    } else if ( n != hostNhood && N_age[a][n] > 0 ) { // different neighborhood, but need people there
	      rTrans += ( 1.0 - RHO_H ) * simParsPtr->get_serotypePar_ij( BETA_INDEX, s ) * n_weight*(double)( numNhoodInfecteds )/(double)( N_age[a][n] ) * alpha[hostAge][a][n];
	    }
	  } // end if n_weight > 0
	} // end for each neighborhood
      } // end for every age
	
      prInf = susc_z * rTrans + simParsPtr->get_serotypePar_ij( IMMIGRATION_INDEX, s );
      prInf = 1.0 - exp( -prInf * EPID_DELTA_T ); 
      if ( r01(rng) < prInf ) {
	infectionTime = r01(rng) * EPID_DELTA_T + t;
	if ( infectionTime < (*it)->getDOD() ) {
	  addEvent( infectionTime, INFECTION_EVENT, hostID, s );
	}
      }
    } // end for each serotype
  } // end for each host
#endif
}

string Simulation::d2str( double d ) {
  std::stringstream t;
  t << d;
  return t.str();
}

string Simulation::makeName( string suffix ) {
  string thisCtr = d2str( simID );
  string thisTr = d2str( treatment );
  string thisName = "tr_" + thisTr + "_sim_" + thisCtr + "_" + suffix;
  return thisName;
}

string Simulation::makeBigName( string suffix, int index ) {
  string thisName = "tr_" + d2str( treatment ) + "_sim_" + d2str( simID ) + "_" + suffix + "_" + d2str(index);
  return thisName;
}

string Simulation::makeBiggerName( string suffix1, int index1, string suffix2, int index2 ) {
  string thisName = "tr_" + d2str( treatment ) + "_sim_" + d2str( simID ) + "_" + suffix1 + "_" + d2str(index1) + "_" + suffix2 + "_" + d2str(index2);
  return thisName;
}

void Simulation::addEvent( double et, int eid, int hid, int s ) {
  while ( currentEvents.find( Event(et) ) != currentEvents.end() ) {
    et += pow(10,APPROX_NOW);
    cout << "\tSimulation is adjusting event time to prevent collision (host id " << hid << ", event id " << eid << ", strain " << s << ", event time " << et << ").\n";
  }
  Event thisEvent( et, eid, hid, s ); 
  currentEvents.insert( thisEvent );
}
