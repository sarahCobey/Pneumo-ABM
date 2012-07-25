/*
* Event.cpp
* Pneumo-ABM - S. Cobey
*
*/

#include <cstdlib>
#include <set>

#include "Event.h"

Event::Event( double t, int eid, int hid ) {
  time = t;
  eventID = eid;
  hostID = hid;
}

Event::~Event() {
}

void Event::add( EventPQ& q ) {
  mq = &q;
  mIt = q.insert( this );
}

void Event::remove() {
  mq.erase( mIt );
  mq = 0;
  mIt = EventPQ::iterator();
}
