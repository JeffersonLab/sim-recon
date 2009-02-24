// DPhoton member functions

#include "DPi0.h"


DPi0::DPi0()
{
   fTags[0] = 0;
   fTags[1] = 0;
}

DPi0::~DPi0()
{
}


// set Pi0 bits with respect to the photon detection 
void DPi0::setChildrenTag(unsigned int tag1, unsigned int tag2 )
{
   fTags[0] = tag1;
   fTags[1] = tag2;
}

// set Pi0 bits with respect to the photon detection 
void DPi0::setChildrenID(oid_t id1, oid_t id2 )
{
   fIDs[0] = id1;
   fIDs[1] = id2;
}

