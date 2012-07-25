/*
*
* Containers.h
* Pneumo-ABM - S. Cobey
*
*/

#ifndef CONTAINERS_H
#define CONTAINERS_H

#include <cstdlib>
#include <set>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/mem_fun.hpp>
#include <boost/multi_index/composite_key.hpp> 
#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/tuple/tuple_io.hpp>

#include "Host.h"
#include "Event.h"

using namespace boost::multi_index; 

#define BOOST_MULTI_INDEX_ENABLE_INVARIANT_CHECKING
#define BOOST_MULTI_INDEX_ENABLE_SAFE_MODE

struct age{};
struct household{};
struct aeh{};
struct eh{};
struct ah{};
struct an{};
struct inf{};

typedef multi_index_container<
  boost::shared_ptr< Host >,
  indexed_by< 
    hashed_unique< const_mem_fun<Host,int,&Host::getID> >, // 0 - ID index
    ordered_non_unique< tag<age>,const_mem_fun<Host,int,&Host::getAgeInY> >, // 1 - Age index
    hashed_non_unique< tag<household>,const_mem_fun<Host,int,&Host::getHousehold> >, // 2 - Household index
    ordered_non_unique< // 3 - Eligible by age & household
      tag<aeh>,
      composite_key<
        Host,
        const_mem_fun<Host,int,&Host::getAgeInY>,
        const_mem_fun<Host,bool,&Host::isEligible>,
        const_mem_fun<Host,int,&Host::getHousehold>
	>
      >,
    ordered_non_unique< // 4 - Eligible by household (all single adults)
      tag<eh>,
      composite_key<
	Host,
	const_mem_fun<Host,bool,&Host::isEligible>,
	const_mem_fun<Host,int,&Host::getHousehold>
	>
      >,
    ordered_non_unique< // 5 - Neighborhood & age
      tag<an>,
      composite_key<
        Host,
        const_mem_fun<Host,int,&Host::getNeighborhood>,
        const_mem_fun<Host,int,&Host::getAgeInY>
        >
      >,
    ordered_non_unique< // 6 - Household & age
      tag<ah>,
      composite_key<
	Host,
	const_mem_fun<Host,int,&Host::getHousehold>,
	const_mem_fun<Host,int,&Host::getAgeInY>
	>
      >,
    ordered_non_unique< tag<inf>,const_mem_fun<Host,bool,&Host::isInfected> > // 7 - Carriage status
   > // end indexed_by
  > HostContainer; 

typedef HostContainer::nth_index<0>::type HostsByID;
typedef HostContainer::nth_index<1>::type HostsByAge;
typedef HostContainer::nth_index<2>::type HostsByHH;
typedef HostContainer::nth_index<3>::type HostsByAEH;
typedef HostContainer::nth_index<4>::type HostsByEH;
typedef HostContainer::nth_index<5>::type HostsByAN;
typedef HostContainer::nth_index<6>::type HostsByAH;
typedef HostContainer::nth_index<7>::type HostsByInf;

struct updateAge {
  void operator() (boost::shared_ptr<Host> ptr) {
    ptr->incrementAge(); 
  }
};

struct updateVaccination {
  void operator() (boost::shared_ptr<Host> ptr) {
    ptr->getVaccinated();
  }
};

struct updateHousehold {   
  updateHousehold( int h ) : h_(h) {}
  void operator() (boost::shared_ptr<Host> ptr) {
    ptr->setHousehold( h_ );
  }
private:
  int h_;
};

struct updatePartner{
  updatePartner( int b ) : b_(b) {}
  void operator() (boost::shared_ptr<Host> ptr ) {
    ptr->setPartner( b_ );
  }
private:
  int b_;
};

struct updateFledge{ 
  updateFledge( bool b ) : b_(b) {}
  void operator() (boost::shared_ptr<Host> ptr ) {
    ptr->setFledge( b_ );
  }
private:
  bool b_;
};

struct updateCarriage{
updateCarriage( int s, double t, EventPQ & ce, boost::mt19937& rng ) : s_(s), t_(t), ce_(ce), rng_(rng) {}
  void operator() (boost::shared_ptr<Host> ptr ) {
    ptr->becomeInfected( s_, t_, ce_, rng_ );
  }
private:
  int s_;
  double t_;
  EventPQ & ce_;
  boost::mt19937 & rng_;
};

struct updateRecovery{
updateRecovery( int s, double t, EventPQ & ce ) : s_(s), t_(t), ce_(ce) {}
  void operator() (boost::shared_ptr<Host> ptr ) {
    ptr->recover( s_, t_, ce_ );
  }
private:
  int s_;
  double t_;
  EventPQ & ce_;
};

struct updateInf {
updateInf( bool b ) : b_(b) {}
  void operator() (boost::shared_ptr<Host> ptr ) {
    ptr->setInf( b_ );
  }
private:
  bool b_;
};

struct updateNeighborhood{
updateNeighborhood( int n ) : n_(n) {}
  void operator() (boost::shared_ptr<Host> ptr) {
    ptr->setNeighborhood( n_ );
  }
private:
  int n_;
};

typedef std::set<int> HHSet;

#endif
