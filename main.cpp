/*
*
* main.cpp
* Pneumo-ABM - S. Cobey
*
*/

#include <cstdlib>
using namespace std;

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <sstream>

#include "Parameters.h"
#include "Host.h"
#include "Event.h"
#include "Rdraws.h"
#include "Simulation.h"
#include "SimPars.h"

#define HFLU_INDEX (INIT_NUM_STYPES-1)

// Function prototypes
void adjustTreatment( int treatmentNumber, double treatment, int simNumber );
void initializeBeta( int treatmentNumber, double beta, int simNumber );
void printAssumptions( void );
double adjustBeta( double preve, double weight, double sb, int treatmentNumber, int simNumber );
string d2str( double d );
string makeName( int treatmentIdx, int simIdx, string suffix );
void printTotalTime( time_t t1, time_t t2 );

int main()
{
  int treatmentNumber;
  double treatment;
  int simNumber;
  cin >> treatmentNumber;
  cin >> treatment;
  cin >> simNumber;
  cout << "Treatment number " << treatmentNumber << ", treatment value " << treatment << ", simulation number " << simNumber << endl;
  printAssumptions();

  double startingBeta;
  double betaTable[100][100]; // buffer

  // Get beta values
  std::ifstream thisFile;
  thisFile.open( "Betas_used.txt", ios::in );
  if ( !thisFile ) {
    cerr << "Error reading Betas_used.txt." << endl;
    exit(1);
  }
  double thisVal;
  int b = 0;
  while ( !thisFile.eof() ) { // Should instead grab beta closest to treatment value in Betas_used
    thisFile >> thisVal;
    betaTable[ b ][ 1 ] = thisVal;
    b++;
  }

  // Get corresponding treatment value (sigma)
  std::ifstream thisFileT;
  thisFileT.open( "Treatments.txt",ios::in );
  if ( !thisFileT ) {
    cerr << "Error reading Treatments.txt." << endl;
    exit(1);
  }
  double thisValT;
  int bT = 0;
  double diff = 10000.0;
  int bestTreatment = 0;
  while ( !thisFileT.eof() ) {
    thisFileT >> thisValT;
    if ( abs(treatment-thisValT) < diff ) { // Finds index of previous treatment most similar to this one
      bestTreatment = bT;
      diff = abs(treatment-thisValT);
    }
    betaTable[ bT ][ 0 ] = thisValT;
    bT++;
  }

  startingBeta = betaTable[ bestTreatment ][ 1 ];
  cout << "Best match to treatment value " << treatment << " is beta = " << startingBeta << " (associated treatment value is " << betaTable[bestTreatment][0] << ")" << endl;
  initializeBeta( treatmentNumber, startingBeta, simNumber );
  adjustTreatment( treatmentNumber, treatment, simNumber );

#ifdef MATCH_PREVALENCE
  cout << "Treatment #" << treatmentNumber << " and simulation #" << simNumber << ":" << endl;
  int matchAttempts = 0;
  double error = 0.0;
  double usedBeta = startingBeta;
  double prevError = 10.0;
  double oldPrevError = 1.0;
  double thisBeta = 0.0;
  double weight = INIT_WEIGHT;
  while ( matchAttempts < MAX_MATCH_ATTEMPTS && abs(prevError) > PREV_ERROR_THOLD ) {
    SimPars thesePars( treatmentNumber, simNumber );
    SimPars * spPtr = &thesePars;
    Simulation thisSim( treatmentNumber, simNumber, spPtr );
    cout << "  Attempt #" << matchAttempts + 1 << endl;
    if ( prevError == 10.0 ) {
      usedBeta = startingBeta;
    }
    thisSim.runDemSim();
    prevError =  thisSim.runTestEpidSim() - TARGET_PREV;
    cout << "Prevalence error=" << prevError << "; weight=" << weight << endl;
    if ( abs(prevError) > PREV_ERROR_THOLD ) {
      if ( prevError*oldPrevError < 0 ) { // if error changed signs and overshot, reduce weight
	weight *= COOL_DOWN;
      } else if ( abs(TARGET_PREV - prevError)/abs(TARGET_PREV - oldPrevError) > TEMP_THOLD ) { // if climbing too slowly, increase weight
	weight *= WARM_UP;
      }
      thisBeta = adjustBeta( prevError, weight, startingBeta, treatmentNumber, simNumber );
      usedBeta = thisBeta;
      oldPrevError = prevError;
    }
    matchAttempts++;
  }
  error = prevError;

  if ( prevError < PREV_ERROR_THOLD ) {
    cout << "\tBeginning simulation #" << simNumber << endl;
    time_t tic;
    tic = time ( NULL );
    SimPars thesePars( treatmentNumber, simNumber );
    SimPars * spPtr = &thesePars;
    Simulation thisSim( treatmentNumber, simNumber, spPtr );
    thisSim.runDemSim();
    thisSim.runEpidSim();
    time_t toc;
    toc = time ( NULL );
    printTotalTime(tic,toc);
  } else {
    cout << "Acceptable prevalence not found." << endl;
  }

  // Output used beta
  std::ofstream thisBetaStream;
  std::ofstream errorStream;
  string betaFile = makeName( treatmentNumber, simNumber, "beta_used" );
  string errorFile = makeName( treatmentNumber, simNumber, "prevalence_errors" );
  thisBetaStream.open( betaFile.c_str(),ios::out );
  errorStream.open( errorFile.c_str(),ios::out );
  thisBetaStream << treatment << "\t" << usedBeta;
  errorStream << error;
  thisBetaStream.close();
  errorStream.close();

#else // Not trying to match prevalence

  cout << "\tBeginning simulation #" << simNumber << endl;
  time_t tic;
  tic = time ( NULL );
  SimPars thesePars( treatmentNumber, simNumber );
  SimPars * spPtr = &thesePars;
  Simulation thisSim( treatmentNumber, simNumber, spPtr );
  thisSim.runDemSim();
  thisSim.runEpidSim();
  time_t toc;
  toc = time ( NULL );
  printTotalTime(tic,toc);

#endif
  return 0;
} 

double adjustBeta( double preve, double w, double sb, int treatmentNumber, int simNumber ) {
  double newBetas[ INIT_NUM_STYPES ];
  for ( int b = 0; b < INIT_NUM_STYPES; b++ ) {
    newBetas[ b ] = 0.0;
  }

  // Either write starting betas or adjust old betas
  if ( preve == 10.0 ) {
    for ( int b = 0; b < HFLU_INDEX; b++ ) {
      newBetas[ b ] = sb;
    }
    newBetas[ HFLU_INDEX ] = HFLU_BETA;
  } else {
    std::ifstream thisFile3;
    string filename = makeName( treatmentNumber, simNumber, "BETA");
    thisFile3.open( filename.c_str(), ios::in );
    if ( !thisFile3 ) {
      cerr << "Error reading " << filename << endl;
      exit(1);
    }
    double thisVal;
    int b = 0;
    while ( !thisFile3.eof() ) {
      thisFile3 >> thisVal;
      cout << "beta=" << thisVal << ";";
      if ( b < HFLU_INDEX ) {
	newBetas[ b ] = thisVal*(1.0 - w*preve);
      } else {
	newBetas[ b ] = HFLU_BETA;
      }
      cout << "newBetas[" << b << "]=" << newBetas[b] << endl;   
      b++;
    }
    thisFile3.close();
  }

  std::ofstream betaStream;
  string filename2 = makeName( treatmentNumber, simNumber, "BETA");
  betaStream.open( filename2.c_str(),ios::out );
  for ( int b = 0; b < INIT_NUM_STYPES; b++ ) {
    betaStream << newBetas[ b ] << "\t";
  }
  betaStream.close();

  double returnBeta = newBetas[0];
  return returnBeta;
}

void initializeBeta( int treatmentNumber, double beta, int simNumber ) {
  std::ofstream betaStream;
  string filename2 = makeName( treatmentNumber, simNumber, "BETA");
  betaStream.open( filename2.c_str(),ios::out );
  for ( int b = 0; b < HFLU_INDEX; b++ ) {
    betaStream << beta << "\t";
  }
  betaStream << HFLU_BETA;
  betaStream.close();
}

void adjustTreatment( int treatmentNumber, double treatment, int simNumber ) {
  double thisXI[ INIT_NUM_STYPES ][ INIT_NUM_STYPES ];
  for ( int i = 0; i < INIT_NUM_STYPES; i++) {
    for ( int j = 0; j < INIT_NUM_STYPES; j++ ) {
      thisXI[i][j] = 0;
    }
  }

  for ( int i = 0; i < HFLU_INDEX; i++ ) {
	thisXI[i][i] = treatment;
  }
  thisXI[ HFLU_INDEX ][ HFLU_INDEX ] = HFLU_SIGMA;

  std::ofstream xiStream;
  string XIFile = makeName( treatmentNumber, simNumber, "XI" );
  xiStream.open( XIFile.c_str(),ios::out );
  for ( int i = 0; i < INIT_NUM_STYPES; i++ ) {
    for ( int j = 0; j < INIT_NUM_STYPES; j++ ) {
      xiStream << thisXI[i][j] ;
      if ( (j+1) % INIT_NUM_STYPES == 0 ) {
	xiStream << endl;
      } else {
	xiStream << "\t";
      }
    }
  }
  xiStream.close();
}

string d2str( double d ) {
  std::stringstream t;
  t << d;
  return t.str();
}

string makeName( int treatmentIdx, int simIdx, string suffix ) {
  string thisCtr = d2str( simIdx );
  string thisTr = d2str( treatmentIdx );
  string thisName = "tr_" + thisTr + "_sim_" + thisCtr + "_" + suffix;
  return thisName;
}

void printAssumptions() {
  cout << "Model assumptions:" << endl;
#ifdef NO_HHOLDS
  cout << "\t--no household structure" <<endl;
#else 
  cout << "\t--households present" << endl;
#endif

#ifdef NO_AGE_ASSORT
  cout << "\t--no age-assortative mixing" << endl;
#else
  cout << "\t--age-assortative mixing" << endl;
#endif

#ifdef MATCH_PREVALENCE
  cout << "\t--matching prevalence" << endl;
#else
  cout << "\t--not matching prevalence"<< endl;
#endif
}

void printTotalTime( time_t t1, time_t t2 ) {
    double ttSec = t2-t1;
    if ( ttSec < 60 ) {
      cout << "  Total simulation time: " << ttSec << " seconds" << endl;
    } else if ( ttSec < 3600 ) {
      int nmin = floor(ttSec/60.0);
      cout << "  Total simulation time: " << nmin << " min " << ttSec - nmin*60 << " seconds " << endl; 
    } else {
      int nhour = floor(ttSec/3600.0);
      int nmin = floor((ttSec - nhour*3600)/60);
      int nsec = ttSec - nhour*3600 - nmin*60;
      cout << "  Total simulation time: " << nhour << " h " << nmin << " m " << " s " << endl;
    }
}

