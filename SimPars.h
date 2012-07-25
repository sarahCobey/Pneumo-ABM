/*
*
* SimPars.h
* Pneumo-ABM - S. Cobey
*
*/

#ifndef SIMPARS_H
#define SIMPARS_H

#include <cstdlib>
#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>

#include "Parameters.h"

typedef double age_index[ INIT_NUM_AGE_CATS ];
typedef age_index dem_array[ NUM_SOCIODEM_FILES ];
typedef age_index waifw_array[ INIT_NUM_AGE_CATS ];
typedef double kids_index[ PARITY_BUFFER ]; 

typedef double stype_index[ INIT_NUM_STYPES ];
typedef stype_index epid_array[ NUM_EPID_FILES + 1 ];
typedef stype_index xi_array[ INIT_NUM_STYPES ]; 


class SimPars 
{
 public:
  SimPars( int tID, int sID );
  ~SimPars( void );

  // PUBLIC FUNCTION PROTOYPES
  const age_index& get_demPMF_row( int ) const;
  const double get_serotypePar_ij( int, int ) const;
  const double get_XI_ij( int, int ) const;
  const double get_waifw_ij( int, int ) const;
  const kids_index& get_parity( void ) const;
  const double get_neighbor( int, int ) const;
  const double get_normalized_neighbor( int, int ) const;
  const double get_Hflu_prob( int ) const;
  const double get_reductions( int ) const;

 private:
  // SIMPARS OBJECTS
  int treatmentID;
  int simID;
  dem_array demPMFs;
  waifw_array WAIFW; 
  waifw_array normWAIFW;
  epid_array serotypePars;
  xi_array XI; 
  kids_index parity_pmf;
  double neighborhoods[ NUM_NEIGHBORHOODS ][ NUM_NEIGHBORHOODS ];
  double NNN[ NUM_NEIGHBORHOODS ][ NUM_NEIGHBORHOODS ];
  double Hflu_probs[ INIT_NUM_STYPES - 1 ]; 
  double reductions[ INIT_NUM_STYPES ]; // holds (1 - MAX_REDUCTION) for each strain

  // PRIVATE FUNCTION PROTOTYPES
  void initializeDemInput();
  void initializeEpidInput();
  string makeName( int treatmentIdx, int simIdx, string suffix );
  string d2str( double d );
};

#endif
