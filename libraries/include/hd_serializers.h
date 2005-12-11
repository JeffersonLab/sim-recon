//
// Auto-generated serializer methods:
// This file was generated from the file hd_ana_d.xml
// on Fri Dec  9 22:27:09 EST 2005
//
// Command line options: hd_ana_d.xml h=hd_classes.h o=hd_serializers 
//


// -------------------- typeid2name ----------------
#include <typeinfo>
const char* JILtypeid2name(const std::type_info *t);


// -------------------- JILMakeObject ----------------
#include <string>
class JILObjectRecord;
class JILStream;
JILObjectRecord* JILMakeObject(const char* type_name, JILStream *s, std::string tag);

#ifndef _JIL_SERIALIZERS_H_
#define _JIL_SERIALIZERS_H_

#include <JILStream.h>
#include "hd_classes.h"


// ---------------- Function Prototypes ------------------
JILStream& operator<<(JILStream &s, const DBCALHit &c);
JILStream& operator>>(JILStream &s, DBCALHit* &c);
JILStream& operator<<(JILStream &s, const DCDCHit &c);
JILStream& operator>>(JILStream &s, DCDCHit* &c);
JILStream& operator<<(JILStream &s, const DFCALGeometry &c);
JILStream& operator>>(JILStream &s, DFCALGeometry* &c);
JILStream& operator<<(JILStream &s, const DFCALHit &c);
JILStream& operator>>(JILStream &s, DFCALHit* &c);
JILStream& operator<<(JILStream &s, const DFCALMCResponse &c);
JILStream& operator>>(JILStream &s, DFCALMCResponse* &c);
JILStream& operator<<(JILStream &s, const DFCALShower &c);
JILStream& operator>>(JILStream &s, DFCALShower* &c);
JILStream& operator<<(JILStream &s, const DFDCHit &c);
JILStream& operator>>(JILStream &s, DFDCHit* &c);
JILStream& operator<<(JILStream &s, const DGeometry &c);
JILStream& operator>>(JILStream &s, DGeometry* &c);
JILStream& operator<<(JILStream &s, const DHDDMForwardShower &c);
JILStream& operator>>(JILStream &s, DHDDMForwardShower* &c);
JILStream& operator<<(JILStream &s, const DHDDMTOFHit &c);
JILStream& operator>>(JILStream &s, DHDDMTOFHit* &c);
JILStream& operator<<(JILStream &s, const DHDDMTOFTruth &c);
JILStream& operator>>(JILStream &s, DHDDMTOFTruth* &c);
JILStream& operator<<(JILStream &s, const DMCThrown &c);
JILStream& operator>>(JILStream &s, DMCThrown* &c);
JILStream& operator<<(JILStream &s, const DMCTrackHit &c);
JILStream& operator>>(JILStream &s, DMCTrackHit* &c);
JILStream& operator<<(JILStream &s, const DObject &c);
JILStream& operator>>(JILStream &s, DObject* &c);
JILStream& operator<<(JILStream &s, const DParameter &c);
JILStream& operator>>(JILStream &s, DParameter* &c);
JILStream& operator<<(JILStream &s, const DQuickFit &c);
JILStream& operator>>(JILStream &s, DQuickFit* &c);
JILStream& operator<<(JILStream &s, const DTOFGeometry &c);
JILStream& operator>>(JILStream &s, DTOFGeometry* &c);
JILStream& operator<<(JILStream &s, const DTOFHit &c);
JILStream& operator>>(JILStream &s, DTOFHit* &c);
JILStream& operator<<(JILStream &s, const DTOFMCResponse &c);
JILStream& operator>>(JILStream &s, DTOFMCResponse* &c);
JILStream& operator<<(JILStream &s, const DTOFPoint &c);
JILStream& operator>>(JILStream &s, DTOFPoint* &c);
JILStream& operator<<(JILStream &s, const DTrack &c);
JILStream& operator>>(JILStream &s, DTrack* &c);
JILStream& operator<<(JILStream &s, const DTrackCandidate &c);
JILStream& operator>>(JILStream &s, DTrackCandidate* &c);
JILStream& operator<<(JILStream &s, const DTrackEfficiency &c);
JILStream& operator>>(JILStream &s, DTrackEfficiency* &c);
JILStream& operator<<(JILStream &s, const DTrackHit &c);
JILStream& operator>>(JILStream &s, DTrackHit* &c);

// -------------------- JILMyDictionary ----------------
const char* JILMyDictionary(void);

#endif // _JIL_SERIALIZERS_H_

