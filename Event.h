/*
*
* Event.h
* Pneumo-ABM - S. Cobey
*
*/

#ifndef EVENT_H
#define EVENT_H

#include <cstdlib>
#include <set>
#include <functional>

#define DEATH_EVENT 0
#define FLEDGE_EVENT 1
#define PAIR_EVENT 2
#define DIVORCE_EVENT 3
#define BIRTH_EVENT 4
#define BIRTHDAY 5
#define INFECTION_EVENT 6
#define RECOVERY_EVENT 7
#define VACCINATION 8

struct Event {

 public:
  explicit Event(double t) : time(t), eventID(), hostID(), s() {}
Event(double t, int eid, int hid, int stype) : time(t), eventID( eid ), hostID( hid ), s(stype) {}

  bool operator < ( const Event & rhs ) const {
    return ( time < rhs.time );
  }

  double time;
  int eventID;
  int hostID;
  int s; // strain, used in infection events
};

typedef std::multiset< Event, std::less< Event > > EventPQ;

#endif

