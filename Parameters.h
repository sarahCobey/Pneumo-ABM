/*
*
* Parameters.h
* Pneumo-ABM - S. Cobey
*
*/

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <cmath>

// All temporal units are days
// Haemophilus influenzae is the last serotype with index INIT_NUM_STYPES-1

// MODEL OPTIONS
//#define MATCH_PREVALENCE // turn off to run with transmission rates in Betas_used.txt and Treatments.txt
#define NO_HHOLDS // if defined, contact rates are independent of household status
#define NO_AGE_ASSORT // if defined, contact rates are independent of host age
#define SIM_PCV // if on, introduces vaccine

// SIMULATION PARAMETERS
// ...input:
const double DEM_SIM_LENGTH = (double)300*365.0; // period of demographic burn-in (no epid dynamics); should be a multiple of strobing interval
const double EPID_SIM_LENGTH = (double)150.0*365.0; // duration of epid dynamics; total sim time is EPID_SIM_LENGTH + DEM_SIM_LENGTH
const double EPID_DELTA_T = (double)1.0; // time step for calculating force of colonization (in days); keep low (~1 day) to avoid error
const double APPROX_NOW = -10.0; // pow(10,APPROX_NOW) used for adjusting times to prevent event collisions

// ...fitting of transmission rate (beta) (parameters below are ignored unless MATCH_PREVALENCE defined)
const double TARGET_PREV = 0.40; // target prevalence in kids <5 y
const double PREV_ERROR_THOLD = 0.01; // allowed error to fit target prevalence (0.01 = target prevalence +/- 1%) 
const int NUM_TEST_SAMPLES = 20; // number of epid strobes used to determine if target criterion has been met (prevalence is average of samples) 
const double COOL_DOWN = 0.8; // decrease in jump size if vacillating around target prevalence prevalence
const double WARM_UP = 1.1; // increase in jump size if undershooting target prevalence
const double TEMP_THOLD = 0.7; // if prevalence is less that TEMP_THOLD fraction of target, jump size (temperature of search) increases (WARM_UP)
const double TEST_EPID_SIM_LENGTH = (double)150*365.0; // duration of epid dynamics for simulations to fit prevalence
const double INIT_WEIGHT = 0.8; // initial coefficient for jump size (jump size = INIT_WEIGHT * difference in prevalence)

// ...output:
const bool CHECK_INPUT = true; // performs non-comprehensive logic checks on input
const double STROBE_EPID = (double)365.0; // interval (in time) between epid observations - Matlab post-processing assumes 365
const double STROBE_DEM = (double)365.0; // interval (in time) between demographic observations - see above
const double PROGRESS_INTERVAL = 10.0; // % interval at which to report progress to screen for each component
const int COCOL_AGE_LIMIT = 5; // in *years*; Haemophilus-pneumo and pneumo-pneumo co-colonization stats printed for hosts <COCOL_AGE_LIMIT

// SOCIODEMOGRAPHIC PARAMETERS
const int N0 = 1000; // initial population size
const double MATURITY_AGE = (double)15.0; 
const int TSTEPS_AGE = 365; // EPID_DELTA_T per age
const int INIT_NUM_AGE_CATS = 111; // if NO_AGE_ASSORT *not* defined, number of age categories (assume categories are YEARS), older ages borrow rates from NUM_AGES
const char * const WAIFW_FILENAME = { "WAIFW.txt" }; // unless NO_AGE_ASSORT defined, age-associated contact matrix (will be normalized)
const int NUM_NEIGHBORHOODS = 1; // do not change--multiple neighborhoods deprecated
const int NUM_SOCIODEM_FILES = 5; // used below
const char * const SOCIODEM_FILENAMES[ NUM_SOCIODEM_FILES ] = { 
  "LSPAN_PMF.txt", // probability mass function of lifespans, indexed by age (in y)
  "FLEDGE_PMF.txt",  // pmf of leaving home of origin (if not already paired), indexed by age (in y)
  "PAIR_PMF.txt",  // pmf of initiating a pairing, by age (in y)
  "BIRTH_AGE_PMF.txt",  // pmf of age at reproduction (in y)
  "INIT_AGE_PMF.txt" // initial age distribution (in y)
};
enum { LSPAN_INDEX, FLEDGE_INDEX, PAIR_INDEX, BIRTH_AGE_INDEX, INIT_AGE_INDEX }; // MUST match file order above
const double INIT_FRAC_PAIRED = 0.40; // fraction of mature adults that are partnered 
const double PROB_PAIR = 0.9; // fraction of the population that should be paired at some point in life
const double STD_AGE_PAIR = 3.0; // half of acceptable range in ages between partners
const char * const PARITY_FILENAME = { "PARITY_PMF.txt" }; // probability of having [i] kids
const char * const NEIGHBORHOOD_FILENAME = { "NEIGHBORHOODS.txt" }; // adjacency matrix for neighborhoods
const char * const HFPROB_FILENAME = { "HFLU_PROBS.txt" }; // probability of serotype being immediately cleared in presence of Hflu
const int HHOLD_SIZE_BUFFER = 30; // buffer for maximum household size
const int PARITY_BUFFER = 50; // max number of offspring permitted
const int DATE_BUFFER = 100; // max number of searches for partner w/in ideal age range; can usually ignore
const double ERR_EPSILON = 0.00007; // acceptable remainder in sum of PMFs; if less than this amount, will add to highest rate class

// EPIDEMIOLOGICAL PARAMETERS
// ...initialization
const int INIT_NUM_STYPES = 26; // = initial number of pneumo serotypes + Haemophilus influenzae
const int NUM_EPID_FILES = 3;
const char * const EPID_FILENAMES[ NUM_EPID_FILES ] = { 
  "INIT_INFECTEDS.txt", // initial fraction of population colonized with each serotype and H. influenzae
  "IMMIG_RATES.txt", // external immigration rate for each serotype and H. influenzae
  "M_DURATION_INFECTION.txt", // mean strain-specific intrinsic duration of infection
};
enum { INIT_INFECTEDS_INDEX, IMMIGRATION_INDEX, MEAN_DURATION_INDEX, BETA_INDEX }; // MUST match file order above		

// ...transmission
const double BASE_DURATION = 25.0; // minimum mean intrinsic duration of carriage
const int COINFECTION_BUFFER = 100; // maximum possible number of coinfections
const double WAIFW_NONZERO = 0.000001; // unless NO_AGE_ASSORT defined, replacement for WAIFW entry ij = 0 (allows w/in household contacts)
const double RHO_H = 0.4; // unless NO_HHOLDS defined, fraction of transmission within household
const double HFLU_BETA = 0; // transmission rate of H. influenzae

// ...immunity
const char * const XI_FILENAME = { "XI.txt" }; // if VARY_SSI not defined, strain-specific cross-immunity (from col j to row i)
const double REC_EPS = 0.25; // coefficient (epsilon) for reduction in duration from current & past carriage
const double MAX_REDUCTION = 0.25; // FOI reduced by (1-MAX_REDUCTION) for res. stype 0; interpolated linearly to zero for others; if 0, no competition from resident
const double HFLU_SIGMA = 0.3; // reduction in susceptibility to H. flu if host has carried it before
const double RSCC_HFLU = 1.0; // reduction in susceptibility to H. flu if carrying H. flu (ensures single-strain dynamics w/o co-colonization)
const double EPSILON = 0.0001; // future time (in days) to 'instantaneous' recovery 

// ...vaccination (ignored unless SIM_PCV defined)
const double VACCINATION_START = (double)100.0*365.0; // days after start of epid simulation to begin vaccination
const double VACCINE_EFFICACY = 0.6; // percent reduction in susceptibility to serotypes in vaccine (susceptibility is max(1-vaccine efficacy,XI))
const double VACCINE_AGE = 30*6; // in days--host age at which vaccinated
const int NUM_VACCINE_SEROTYPES = 5; // must match next two variables
const int VACCINE_SEROTYPES[ NUM_VACCINE_SEROTYPES ] = { 0, 3, 7, 10, 20 }; // must match nonzero entries below
const bool IN_VACCINE[ 25 ] = { 
	1, 0, 0, 1, 0,
	0, 0, 1, 0, 0,
	1, 0, 0, 0, 0,
	0, 0, 0, 0, 0,
	1, 0, 0, 0, 0
};

#endif
