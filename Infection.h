/*
*
* Infection.h
* Pneumo-ABM - S. Cobey
*
*/

#ifndef INFECTION_H
#define INFECTION_H

#include <cstdlib>
#include <boost/unordered_map.hpp>

struct Infection {
public:
  explicit Infection( double it, double rt ) : infT( it ), recT( rt ) {}
  double infT;
  double recT;
};

typedef boost::unordered_multimap< int, Infection > InfectionMap;

#endif

