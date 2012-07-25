/*
*
* Simulation.h
* Pneumo ABM - S. Cobey
*
*/

#ifndef SIMULATION_H
#define SIMULATION_H

#include <cstdlib>
#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <boost/random.hpp>

#include "Parameters.h"
#include "Host.h"
#include "Event.h"
#include "Rdraws.h"
#include "Containers.h"
#include "SimPars.h"


class Simulation
{
 public:
  Simulation( int trt, int sid, SimPars * spPtr  );
  ~Simulation();

  // MEMBER FUNCTION PROTOTYPES
  void runDemSim( void );
  void runEpidSim( void ); 
  double runTestEpidSim( void );

 private:
  // SIMULATION OBJECTS
  boost::mt19937 rng;
  double t;
  double demOutputStrobe;
  double epidOutputStrobe;
  double demComplete;
  int idCtr;
  int hholdCtr;
  int simID;
  int eventCtr; // for memory management and debugging
  HostContainer allHosts; // shared_ptr to Hosts
  HHSet allHouseholds; // set of household ids
  int numInfecteds[ INIT_NUM_AGE_CATS ][ INIT_NUM_STYPES ][ NUM_NEIGHBORHOODS ]; // counts number infected with each serotype by age
  EventPQ currentEvents; // current events heap
  SimPars * simParsPtr;
  int treatment;

  // SIMULATION FUNCTIONS
  void ageHost( int id );
  int calcPartnerAge( int a );
  void executeEvent( Event & te );
  void seedInfections( void );
  void killHost( int id );
  void pairHost( int id );
  void partner2Hosts( int id1, int id2 );
  void fledgeHost( int id );
  void birthHost( int id );
  void infectHost( int id, int s );
  void recoverHost( int id, int s );
  void vaccinateHost( int id );
  void calcSI( void );
  string d2str( double d );
  string makeName( string s );
  string makeBigName( string, int );
  string makeBiggerName( string, int, string, int);
  void addEvent( double et, int eid, int hid, int s );
  double calcPrev( void );

  // STREAM MANAGEMENT
  void writeDemOutput();
  void writeEpidOutput();
  void initOutput();
  void closeOutput();
  void writeTheta();

  std::ofstream ageDistStream;
  std::ofstream hhDistStream;
  std::ofstream demTimesStream;
  std::ofstream infectionStream;
  std::ofstream epidTimesStream;
  std::ofstream infectedStream;
  std::ofstream coinfectionHistStream;
  std::ofstream coinfectionHFHistStream;
  std::ofstream thetaStream;
  std::ofstream totCarriageStream;

  std::string ageDistFile;
  std::string hhDistFile; 
  std::string demTimesFile;
  std::string epidTimesFile;
  std::string coinfectionHistFile;
  std::string coinfectionHFHistFile;
  std::string thetaFile;
  std::string totCarriageFile;
};

#endif
