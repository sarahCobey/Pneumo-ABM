/*
*
* Host.h
* Pneumo-ABM - S. Cobey
*
*/

#ifndef HOST_H
#define HOST_H


#include "Event.h"
#include "Infection.h"
#include "Parameters.h"
#include "SimPars.h"
#include <boost/random.hpp>


class Host
{
 public:
  Host( double, double, int, int, int, EventPQ &, SimPars * , boost::mt19937 &); 
  ~Host( void );

  // SET FUNCTION PROTOTYPES
  void setHousehold( int );
  void setNeighborhood( int );
  void setGroup( int );
  void setDOB( double );
  void incrementAge();
  void setPartner( int );
  void setFledge( bool );
  void setInf( bool );
  void getVaccinated();

  // GET FUNCTION PROTOTYPES
  int getAgeInY() const;
  int getID() const;
  int getGroup() const;
  int getHousehold() const;
  int getPartner() const;
  bool isPaired() const;
  bool isAdult() const;
  bool isEligible() const;
  bool hasFledged() const;
  int isInfectedZ( int ) const;
  bool isInfected() const;
  bool isInfectedPneumo() const;
  bool isInfectedHflu() const;
  int totStrains() const;
  double getSusc( int ) const;
  double getDOD() const;
  int getNeighborhood() const;
  int getSummedTheta() const;
  double getDOB() const;

  // MEMBER FUNCTION PROTOTYPES
  void calcLifeHist( double, EventPQ &, double, boost::mt19937 & );
  double calcDeath( boost::mt19937 & );
  void addEvent( double, int, int, EventPQ & );
  void addEvent( double, int, int, int, EventPQ & );
  double calcFledge( boost::mt19937 & );
  double calcPairAge( boost::mt19937 & );
  int calcNumBirths( boost::mt19937 & );
  double calcBirthAge( boost::mt19937 & );
  void becomeInfected( int, double, EventPQ &, boost::mt19937 & );
  void recover( int, double, EventPQ & );
  double calcRecovery( int, double, double, boost::mt19937 & );
  double calcRecovery( int, double, Infection &, boost::mt19937 & );
  void calcSusc( double );

 private:
  double DOB;
  double DOD; 
  int age;
  int id;
  int group;
  int household;
  int neighborhood;
  int partner;
  bool fledge;
  InfectionMap carriage;
  int immune[ INIT_NUM_STYPES ]; 
  int carriageSummary[ INIT_NUM_STYPES ];
  double susc[ INIT_NUM_STYPES ];
  bool inf;
  SimPars * simParsPtr;

};

#endif
