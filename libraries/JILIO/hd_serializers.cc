//
// Auto-generated serializer methods:
// This file was generated from the file hd_ana_d.xml
// on Fri Dec  9 12:52:49 EST 2005
//
// Command line options: hd_ana_d.xml h=hd_classes.h o=hd_serializers 
//

#include "hd_serializers.h"


// -------------------- typeid2name ----------------
const char* JILtypeid2name(const type_info *t)
{
	if(t == &typeid(DBCALHit))return "DBCALHit";
	if(t == &typeid(DCDCHit))return "DCDCHit";
	if(t == &typeid(DFCALGeometry))return "DFCALGeometry";
	if(t == &typeid(DFCALHit))return "DFCALHit";
	if(t == &typeid(DFCALMCResponse))return "DFCALMCResponse";
	if(t == &typeid(DFCALShower))return "DFCALShower";
	if(t == &typeid(DFDCHit))return "DFDCHit";
	if(t == &typeid(DGeometry))return "DGeometry";
	if(t == &typeid(DHDDMForwardShower))return "DHDDMForwardShower";
	if(t == &typeid(DHDDMTOFHit))return "DHDDMTOFHit";
	if(t == &typeid(DHDDMTOFTruth))return "DHDDMTOFTruth";
	if(t == &typeid(DMCThrown))return "DMCThrown";
	if(t == &typeid(DMCTrackHit))return "DMCTrackHit";
	if(t == &typeid(DObject))return "DObject";
	if(t == &typeid(DParameter))return "DParameter";
	if(t == &typeid(DQuickFit))return "DQuickFit";
	if(t == &typeid(DTOFGeometry))return "DTOFGeometry";
	if(t == &typeid(DTOFHit))return "DTOFHit";
	if(t == &typeid(DTOFMCResponse))return "DTOFMCResponse";
	if(t == &typeid(DTOFPoint))return "DTOFPoint";
	if(t == &typeid(DTrack))return "DTrack";
	if(t == &typeid(DTrackCandidate))return "DTrackCandidate";
	if(t == &typeid(DTrackEfficiency))return "DTrackEfficiency";
	if(t == &typeid(DTrackHit))return "DTrackHit";

	if(t == &typeid(short))return "short";
	if(t == &typeid(int))return "int";
	if(t == &typeid(long))return "long";
	if(t == &typeid(unsigned short))return "ushort";
	if(t == &typeid(unsigned int))return "uint";
	if(t == &typeid(unsigned long))return "ulong";
	if(t == &typeid(float))return "float";
	if(t == &typeid(double))return "double";
	if(t == &typeid(string))return "string";

	return "unknown";
}


// -------------------- JILMakeObject ----------------
JILObjectRecord* JILMakeObject(const char* type_name, JILStream *s, std::string tag)
{
	if(!strcmp(type_name,JILtypeid2name(&typeid(DBCALHit))))return new JILObjectRecordT<DBCALHit >(s, tag);
	if(!strcmp(type_name,JILtypeid2name(&typeid(DCDCHit))))return new JILObjectRecordT<DCDCHit >(s, tag);
	if(!strcmp(type_name,JILtypeid2name(&typeid(DFCALGeometry))))return new JILObjectRecordT<DFCALGeometry >(s, tag);
	if(!strcmp(type_name,JILtypeid2name(&typeid(DFCALHit))))return new JILObjectRecordT<DFCALHit >(s, tag);
	if(!strcmp(type_name,JILtypeid2name(&typeid(DFCALMCResponse))))return new JILObjectRecordT<DFCALMCResponse >(s, tag);
	if(!strcmp(type_name,JILtypeid2name(&typeid(DFCALShower))))return new JILObjectRecordT<DFCALShower >(s, tag);
	if(!strcmp(type_name,JILtypeid2name(&typeid(DFDCHit))))return new JILObjectRecordT<DFDCHit >(s, tag);
	if(!strcmp(type_name,JILtypeid2name(&typeid(DGeometry))))return new JILObjectRecordT<DGeometry >(s, tag);
	if(!strcmp(type_name,JILtypeid2name(&typeid(DHDDMForwardShower))))return new JILObjectRecordT<DHDDMForwardShower >(s, tag);
	if(!strcmp(type_name,JILtypeid2name(&typeid(DHDDMTOFHit))))return new JILObjectRecordT<DHDDMTOFHit >(s, tag);
	if(!strcmp(type_name,JILtypeid2name(&typeid(DHDDMTOFTruth))))return new JILObjectRecordT<DHDDMTOFTruth >(s, tag);
	if(!strcmp(type_name,JILtypeid2name(&typeid(DMCThrown))))return new JILObjectRecordT<DMCThrown >(s, tag);
	if(!strcmp(type_name,JILtypeid2name(&typeid(DMCTrackHit))))return new JILObjectRecordT<DMCTrackHit >(s, tag);
	if(!strcmp(type_name,JILtypeid2name(&typeid(DObject))))return new JILObjectRecordT<DObject >(s, tag);
	if(!strcmp(type_name,JILtypeid2name(&typeid(DParameter))))return new JILObjectRecordT<DParameter >(s, tag);
	if(!strcmp(type_name,JILtypeid2name(&typeid(DQuickFit))))return new JILObjectRecordT<DQuickFit >(s, tag);
	if(!strcmp(type_name,JILtypeid2name(&typeid(DTOFGeometry))))return new JILObjectRecordT<DTOFGeometry >(s, tag);
	if(!strcmp(type_name,JILtypeid2name(&typeid(DTOFHit))))return new JILObjectRecordT<DTOFHit >(s, tag);
	if(!strcmp(type_name,JILtypeid2name(&typeid(DTOFMCResponse))))return new JILObjectRecordT<DTOFMCResponse >(s, tag);
	if(!strcmp(type_name,JILtypeid2name(&typeid(DTOFPoint))))return new JILObjectRecordT<DTOFPoint >(s, tag);
	if(!strcmp(type_name,JILtypeid2name(&typeid(DTrack))))return new JILObjectRecordT<DTrack >(s, tag);
	if(!strcmp(type_name,JILtypeid2name(&typeid(DTrackCandidate))))return new JILObjectRecordT<DTrackCandidate >(s, tag);
	if(!strcmp(type_name,JILtypeid2name(&typeid(DTrackEfficiency))))return new JILObjectRecordT<DTrackEfficiency >(s, tag);
	if(!strcmp(type_name,JILtypeid2name(&typeid(DTrackHit))))return new JILObjectRecordT<DTrackHit >(s, tag);

	return NULL;
}


//----------------- DBCALHit -----------------
JILStream& operator<<(JILStream &s, const DBCALHit &c){

	// Send object type to stream
	s<<&typeid(DBCALHit);


	// Write out base class first
	s<<*((DObject*)&c);

	s<<c.module;  // int   (public)
	s<<c.layer;  // int   (public)
	s<<c.sector;  // int   (public)
	s<<c.end;  // DBCALHit::END_t   (public)
	s<<c.E;  // float   (public)
	s<<c.t;  // float   (public)

	// Send end of object flag to stream
	s<<JILStream::END_OBJECT;

	return s;
}

JILStream& operator>>(JILStream &s, DBCALHit* &c){
	// Create object if it doesn't already exist
	if(c == NULL)c = new DBCALHit();

	// Read in base class first
	s>>(DObject*&)c;

	s>>c->module;  // int   (public)
	s>>c->layer;  // int   (public)
	s>>c->sector;  // int   (public)
	s>>c->end;  // DBCALHit::END_t   (public)
	s>>c->E;  // float   (public)
	s>>c->t;  // float   (public)

	return s;
}

//----------------- DCDCHit -----------------
JILStream& operator<<(JILStream &s, const DCDCHit &c){

	// Send object type to stream
	s<<&typeid(DCDCHit);


	// Write out base class first
	s<<*((DObject*)&c);

	s<<c.radius;  // float   (public)
	s<<c.phim;  // float   (public)
	s<<c.dE;  // float   (public)
	s<<c.t;  // float   (public)

	// Send end of object flag to stream
	s<<JILStream::END_OBJECT;

	return s;
}

JILStream& operator>>(JILStream &s, DCDCHit* &c){
	// Create object if it doesn't already exist
	if(c == NULL)c = new DCDCHit();

	// Read in base class first
	s>>(DObject*&)c;

	s>>c->radius;  // float   (public)
	s>>c->phim;  // float   (public)
	s>>c->dE;  // float   (public)
	s>>c->t;  // float   (public)

	return s;
}

//----------------- DFCALGeometry -----------------
JILStream& operator<<(JILStream &s, const DFCALGeometry &c){

	// Send object type to stream
	s<<&typeid(DFCALGeometry);


	// Write out base class first
	s<<*((DObject*)&c);

	// for(unsigned int i0=0; i0<53; i0++)
	//	s.WriteArray(c.m_activeBlock[i0], 53);  // bool[53][53]   (private)
	// s<<c.m_positionOnFace;  // TVector2   (private)
	// for(unsigned int i0=0; i0<53; i0++)
	//	s.WriteArray(c.m_channelNumber[i0], 53);  // int[53][53]   (private)
	// s.WriteArray(c.m_row, 2809);  // int[2809]   (private)
	// s.WriteArray(c.m_column, 2809);  // int[2809]   (private)
	// s<<c.m_numActiveBlocks;  // int   (private)

	// Send end of object flag to stream
	s<<JILStream::END_OBJECT;

	return s;
}

JILStream& operator>>(JILStream &s, DFCALGeometry* &c){
	// Create object if it doesn't already exist
	if(c == NULL)c = new DFCALGeometry();

	// Read in base class first
	s>>(DObject*&)c;

	// for(unsigned int i0=0; i0<53; i0++)
	//	s.ReadArray(c->m_activeBlock[i0], 53);  // bool[53][53]   (private)
	// s>>c->m_positionOnFace;  // TVector2   (private)
	// for(unsigned int i0=0; i0<53; i0++)
	//	s.ReadArray(c->m_channelNumber[i0], 53);  // int[53][53]   (private)
	// s.ReadArray(c->m_row, 2809);  // int[2809]   (private)
	// s.ReadArray(c->m_column, 2809);  // int[2809]   (private)
	// s>>c->m_numActiveBlocks;  // int   (private)

	return s;
}

//----------------- DFCALHit -----------------
JILStream& operator<<(JILStream &s, const DFCALHit &c){

	// Send object type to stream
	s<<&typeid(DFCALHit);


	// Write out base class first
	s<<*((DObject*)&c);

	s<<c.x;  // float   (public)
	s<<c.y;  // float   (public)
	s<<c.E;  // float   (public)
	s<<c.t;  // float   (public)

	// Send end of object flag to stream
	s<<JILStream::END_OBJECT;

	return s;
}

JILStream& operator>>(JILStream &s, DFCALHit* &c){
	// Create object if it doesn't already exist
	if(c == NULL)c = new DFCALHit();

	// Read in base class first
	s>>(DObject*&)c;

	s>>c->x;  // float   (public)
	s>>c->y;  // float   (public)
	s>>c->E;  // float   (public)
	s>>c->t;  // float   (public)

	return s;
}

//----------------- DFCALMCResponse -----------------
JILStream& operator<<(JILStream &s, const DFCALMCResponse &c){

	// Send object type to stream
	s<<&typeid(DFCALMCResponse);


	// Write out base class first
	s<<*((DObject*)&c);

	// s<<c.m_channel;  // int   (private)
	// s<<c.m_E;  // double   (private)
	// s<<c.m_t;  // double   (private)

	// Send end of object flag to stream
	s<<JILStream::END_OBJECT;

	return s;
}

JILStream& operator>>(JILStream &s, DFCALMCResponse* &c){
	// Create object if it doesn't already exist
	if(c == NULL)c = new DFCALMCResponse();

	// Read in base class first
	s>>(DObject*&)c;

	// s>>c->m_channel;  // int   (private)
	// s>>c->m_E;  // double   (private)
	// s>>c->m_t;  // double   (private)

	return s;
}

//----------------- DFCALShower -----------------
JILStream& operator<<(JILStream &s, const DFCALShower &c){

	// Send object type to stream
	s<<&typeid(DFCALShower);


	// Write out base class first
	s<<*((DObject*)&c);

	s<<c.x;  // float   (public)
	s<<c.y;  // float   (public)
	s<<c.E;  // float   (public)
	s<<c.t;  // float   (public)

	// Send end of object flag to stream
	s<<JILStream::END_OBJECT;

	return s;
}

JILStream& operator>>(JILStream &s, DFCALShower* &c){
	// Create object if it doesn't already exist
	if(c == NULL)c = new DFCALShower();

	// Read in base class first
	s>>(DObject*&)c;

	s>>c->x;  // float   (public)
	s>>c->y;  // float   (public)
	s>>c->E;  // float   (public)
	s>>c->t;  // float   (public)

	return s;
}

//----------------- DFDCHit -----------------
JILStream& operator<<(JILStream &s, const DFDCHit &c){

	// Send object type to stream
	s<<&typeid(DFDCHit);


	// Write out base class first
	s<<*((DObject*)&c);

	s<<c.layer;  // int   (public)
	s<<c.module;  // int   (public)
	s<<c.tau;  // float   (public)
	s<<c.z;  // float   (public)
	s<<c.u;  // float   (public)
	s<<c.dE;  // float   (public)
	s<<c.t;  // float   (public)
	s<<c.type;  // int   (public)

	// Send end of object flag to stream
	s<<JILStream::END_OBJECT;

	return s;
}

JILStream& operator>>(JILStream &s, DFDCHit* &c){
	// Create object if it doesn't already exist
	if(c == NULL)c = new DFDCHit();

	// Read in base class first
	s>>(DObject*&)c;

	s>>c->layer;  // int   (public)
	s>>c->module;  // int   (public)
	s>>c->tau;  // float   (public)
	s>>c->z;  // float   (public)
	s>>c->u;  // float   (public)
	s>>c->dE;  // float   (public)
	s>>c->t;  // float   (public)
	s>>c->type;  // int   (public)

	return s;
}

//----------------- DGeometry -----------------
JILStream& operator<<(JILStream &s, const DGeometry &c){

	// Send object type to stream
	s<<&typeid(DGeometry);

	// s<<c.min_run_number;  // int   (protected)
	// s<<c.max_run_number;  // int   (protected)

	// Send end of object flag to stream
	s<<JILStream::END_OBJECT;

	return s;
}

JILStream& operator>>(JILStream &s, DGeometry* &c){
	// Create object if it doesn't already exist
	if(c == NULL)c = new DGeometry();

	// s>>c->min_run_number;  // int   (protected)
	// s>>c->max_run_number;  // int   (protected)

	return s;
}

//----------------- DHDDMForwardShower -----------------
JILStream& operator<<(JILStream &s, const DHDDMForwardShower &c){

	// Send object type to stream
	s<<&typeid(DHDDMForwardShower);


	// Write out base class first
	s<<*((DObject*)&c);

	// s<<c.m_x;  // float   (private)
	// s<<c.m_y;  // float   (private)
	// s<<c.m_E;  // float   (private)
	// s<<c.m_t;  // float   (private)

	// Send end of object flag to stream
	s<<JILStream::END_OBJECT;

	return s;
}

JILStream& operator>>(JILStream &s, DHDDMForwardShower* &c){
	// Create object if it doesn't already exist
	if(c == NULL)c = new DHDDMForwardShower();

	// Read in base class first
	s>>(DObject*&)c;

	// s>>c->m_x;  // float   (private)
	// s>>c->m_y;  // float   (private)
	// s>>c->m_E;  // float   (private)
	// s>>c->m_t;  // float   (private)

	return s;
}

//----------------- DHDDMTOFHit -----------------
JILStream& operator<<(JILStream &s, const DHDDMTOFHit &c){

	// Send object type to stream
	s<<&typeid(DHDDMTOFHit);


	// Write out base class first
	s<<*((DObject*)&c);

	s<<c.orientation;  // int   (public)
	s<<c.end;  // int   (public)
	s<<c.y;  // float   (public)
	s<<c.t;  // float   (public)
	s<<c.E;  // float   (public)

	// Send end of object flag to stream
	s<<JILStream::END_OBJECT;

	return s;
}

JILStream& operator>>(JILStream &s, DHDDMTOFHit* &c){
	// Create object if it doesn't already exist
	if(c == NULL)c = new DHDDMTOFHit();

	// Read in base class first
	s>>(DObject*&)c;

	s>>c->orientation;  // int   (public)
	s>>c->end;  // int   (public)
	s>>c->y;  // float   (public)
	s>>c->t;  // float   (public)
	s>>c->E;  // float   (public)

	return s;
}

//----------------- DHDDMTOFTruth -----------------
JILStream& operator<<(JILStream &s, const DHDDMTOFTruth &c){

	// Send object type to stream
	s<<&typeid(DHDDMTOFTruth);


	// Write out base class first
	s<<*((DObject*)&c);

	s<<c.orientation;  // int   (public)
	s<<c.track;  // int   (public)
	s<<c.primary;  // int   (public)
	s<<c.x;  // float   (public)
	s<<c.y;  // float   (public)
	s<<c.z;  // float   (public)
	s<<c.t;  // float   (public)

	// Send end of object flag to stream
	s<<JILStream::END_OBJECT;

	return s;
}

JILStream& operator>>(JILStream &s, DHDDMTOFTruth* &c){
	// Create object if it doesn't already exist
	if(c == NULL)c = new DHDDMTOFTruth();

	// Read in base class first
	s>>(DObject*&)c;

	s>>c->orientation;  // int   (public)
	s>>c->track;  // int   (public)
	s>>c->primary;  // int   (public)
	s>>c->x;  // float   (public)
	s>>c->y;  // float   (public)
	s>>c->z;  // float   (public)
	s>>c->t;  // float   (public)

	return s;
}

//----------------- DMCThrown -----------------
JILStream& operator<<(JILStream &s, const DMCThrown &c){

	// Send object type to stream
	s<<&typeid(DMCThrown);


	// Write out base class first
	s<<*((DObject*)&c);

	s<<c.type;  // int   (public)
	s<<c.q;  // float   (public)
	s<<c.p;  // float   (public)
	s<<c.E;  // float   (public)
	s<<c.theta;  // float   (public)
	s<<c.phi;  // float   (public)
	s<<c.x;  // float   (public)
	s<<c.y;  // float   (public)
	s<<c.z;  // float   (public)
	s<<c.mass;  // float   (public)

	// Send end of object flag to stream
	s<<JILStream::END_OBJECT;

	return s;
}

JILStream& operator>>(JILStream &s, DMCThrown* &c){
	// Create object if it doesn't already exist
	if(c == NULL)c = new DMCThrown();

	// Read in base class first
	s>>(DObject*&)c;

	s>>c->type;  // int   (public)
	s>>c->q;  // float   (public)
	s>>c->p;  // float   (public)
	s>>c->E;  // float   (public)
	s>>c->theta;  // float   (public)
	s>>c->phi;  // float   (public)
	s>>c->x;  // float   (public)
	s>>c->y;  // float   (public)
	s>>c->z;  // float   (public)
	s>>c->mass;  // float   (public)

	return s;
}

//----------------- DMCTrackHit -----------------
JILStream& operator<<(JILStream &s, const DMCTrackHit &c){

	// Send object type to stream
	s<<&typeid(DMCTrackHit);


	// Write out base class first
	s<<*((DObject*)&c);

	s<<c.r;  // float   (public)
	s<<c.phi;  // float   (public)
	s<<c.z;  // float   (public)
	s<<c.track;  // int   (public)
	s<<c.primary;  // int   (public)
	s<<(int)c.system;  // DetectorSystem_t   (public)

	// Send end of object flag to stream
	s<<JILStream::END_OBJECT;

	return s;
}

JILStream& operator>>(JILStream &s, DMCTrackHit* &c){
	// Create object if it doesn't already exist
	if(c == NULL)c = new DMCTrackHit();

	// Read in base class first
	s>>(DObject*&)c;

	s>>c->r;  // float   (public)
	s>>c->phi;  // float   (public)
	s>>c->z;  // float   (public)
	s>>c->track;  // int   (public)
	s>>c->primary;  // int   (public)
	s>>*(int*)&c->system;  // DetectorSystem_t   (public)

	return s;
}

//----------------- DObject -----------------
JILStream& operator<<(JILStream &s, const DObject &c){

	// Send object type to stream
	s<<&typeid(DObject);

	s<<c.id;  // identifier_t   (public)

	// Send end of object flag to stream
	s<<JILStream::END_OBJECT;

	return s;
}

JILStream& operator>>(JILStream &s, DObject* &c){
	// Create object if it doesn't already exist
	if(c == NULL)c = new DObject();

	s>>c->id;  // identifier_t   (public)

	return s;
}

//----------------- DParameter -----------------
JILStream& operator<<(JILStream &s, const DParameter &c){
	c.Serialize(s);
	return s;
}

void DParameter::Serialize(JILStream &s) const
{

	// Send object type to stream
	s<<&typeid(DParameter);

	s<<key;  // string   (protected)
	s<<value;  // string   (protected)
	s<<isdefault;  // bool   (protected)
	s<<hasdefault;  // bool   (protected)
	s<<printme;  // bool   (protected)
	s<<type;  // DParameter::dataType_t   (protected)

	// Send end of object flag to stream
	s<<JILStream::END_OBJECT;

}

JILStream& operator>>(JILStream &s, DParameter* &c){
	// Create object if it doesn't already exist
	if(c == NULL)c = new DParameter();

	c->Deserialize(s);
	return s;
}

void DParameter::Deserialize(JILStream &s)
{
	s>>key;  // string   (protected)
	s>>value;  // string   (protected)
	s>>isdefault;  // bool   (protected)
	s>>hasdefault;  // bool   (protected)
	s>>printme;  // bool   (protected)
	s>>type;  // DParameter::dataType_t   (protected)

}

//----------------- DQuickFit -----------------
JILStream& operator<<(JILStream &s, const DQuickFit &c){

	// Send object type to stream
	s<<&typeid(DQuickFit);

	s<<c.x0;  // float   (public)
	s<<c.y0;  // float   (public)
	s<<c.q;  // float   (public)
	s<<c.p;  // float   (public)
	s<<c.p_trans;  // float   (public)
	s<<c.phi;  // float   (public)
	s<<c.theta;  // float   (public)
	s<<c.z_vertex;  // float   (public)
	s<<c.chisq;  // float   (public)
	s<<c.dphidz;  // float   (public)
	s<<c.chisq_source;  // DQuickFit::ChiSqSourceType_t   (public)
	// s<<c.hits;  // vector   (protected)
	// s<<c.bfield;  // DMagneticFieldMap   (protected)
	// s<<c.Bz_avg;  // float   (protected)
	// s<<c.z_mean;  // float   (protected)
	// s<<c.phi_mean;  // float   (protected)

	// Send end of object flag to stream
	s<<JILStream::END_OBJECT;

	return s;
}

JILStream& operator>>(JILStream &s, DQuickFit* &c){
	// Create object if it doesn't already exist
	if(c == NULL)c = new DQuickFit();

	s>>c->x0;  // float   (public)
	s>>c->y0;  // float   (public)
	s>>c->q;  // float   (public)
	s>>c->p;  // float   (public)
	s>>c->p_trans;  // float   (public)
	s>>c->phi;  // float   (public)
	s>>c->theta;  // float   (public)
	s>>c->z_vertex;  // float   (public)
	s>>c->chisq;  // float   (public)
	s>>c->dphidz;  // float   (public)
	s>>c->chisq_source;  // DQuickFit::ChiSqSourceType_t   (public)
	// s>>c->hits;  // vector   (protected)
	// if(s.GetPointerFromStreamT(c->bfield))s>>c->bfield;  // DMagneticFieldMap   (protected)
	// s>>c->Bz_avg;  // float   (protected)
	// s>>c->z_mean;  // float   (protected)
	// s>>c->phi_mean;  // float   (protected)

	return s;
}

//----------------- DTOFGeometry -----------------
JILStream& operator<<(JILStream &s, const DTOFGeometry &c){

	// Send object type to stream
	s<<&typeid(DTOFGeometry);


	// Write out base class first
	s<<*((DObject*)&c);

	s<<c.NLONGBARS;  // int   (public)
	s<<c.NSHORTBARS;  // int   (public)
	s<<c.LONGBARLENGTH;  // float   (public)
	s<<c.SHORTBARLENGTH;  // float   (public)
	s<<c.BARWIDTH;  // float   (public)

	// Send end of object flag to stream
	s<<JILStream::END_OBJECT;

	return s;
}

JILStream& operator>>(JILStream &s, DTOFGeometry* &c){
	// Create object if it doesn't already exist
	if(c == NULL)c = new DTOFGeometry();

	// Read in base class first
	s>>(DObject*&)c;

	s>>c->NLONGBARS;  // int   (public)
	s>>c->NSHORTBARS;  // int   (public)
	s>>c->LONGBARLENGTH;  // float   (public)
	s>>c->SHORTBARLENGTH;  // float   (public)
	s>>c->BARWIDTH;  // float   (public)

	return s;
}

//----------------- DTOFHit -----------------
JILStream& operator<<(JILStream &s, const DTOFHit &c){

	// Send object type to stream
	s<<&typeid(DTOFHit);


	// Write out base class first
	s<<*((DObject*)&c);

	s<<c.orientation;  // int   (public)
	s<<c.end;  // int   (public)
	s<<c.y;  // float   (public)
	s<<c.t;  // float   (public)
	s<<c.E;  // float   (public)

	// Send end of object flag to stream
	s<<JILStream::END_OBJECT;

	return s;
}

JILStream& operator>>(JILStream &s, DTOFHit* &c){
	// Create object if it doesn't already exist
	if(c == NULL)c = new DTOFHit();

	// Read in base class first
	s>>(DObject*&)c;

	s>>c->orientation;  // int   (public)
	s>>c->end;  // int   (public)
	s>>c->y;  // float   (public)
	s>>c->t;  // float   (public)
	s>>c->E;  // float   (public)

	return s;
}

//----------------- DTOFMCResponse -----------------
JILStream& operator<<(JILStream &s, const DTOFMCResponse &c){

	// Send object type to stream
	s<<&typeid(DTOFMCResponse);


	// Write out base class first
	s<<*((DObject*)&c);

	s<<c.orientation;  // int   (public)
	s<<c.end;  // int   (public)
	s<<c.y;  // float   (public)
	s<<c.t;  // float   (public)
	s<<c.E;  // float   (public)

	// Send end of object flag to stream
	s<<JILStream::END_OBJECT;

	return s;
}

JILStream& operator>>(JILStream &s, DTOFMCResponse* &c){
	// Create object if it doesn't already exist
	if(c == NULL)c = new DTOFMCResponse();

	// Read in base class first
	s>>(DObject*&)c;

	s>>c->orientation;  // int   (public)
	s>>c->end;  // int   (public)
	s>>c->y;  // float   (public)
	s>>c->t;  // float   (public)
	s>>c->E;  // float   (public)

	return s;
}

//----------------- DTOFPoint -----------------
JILStream& operator<<(JILStream &s, const DTOFPoint &c){

	// Send object type to stream
	s<<&typeid(DTOFPoint);


	// Write out base class first
	s<<*((DObject*)&c);

	s<<c.trackid;  // int   (public)
	s<<c.x;  // float   (public)
	s<<c.y;  // float   (public)
	s<<c.z;  // float   (public)
	s<<c.t;  // float   (public)
	s<<c.dedx;  // float   (public)
	s<<c.nhits;  // int   (public)
	s.WriteArray(c.hits, 16);  // unsigned int[16]   (public)
	s<<c.chisq;  // float   (public)

	// Send end of object flag to stream
	s<<JILStream::END_OBJECT;

	return s;
}

JILStream& operator>>(JILStream &s, DTOFPoint* &c){
	// Create object if it doesn't already exist
	if(c == NULL)c = new DTOFPoint();

	// Read in base class first
	s>>(DObject*&)c;

	s>>c->trackid;  // int   (public)
	s>>c->x;  // float   (public)
	s>>c->y;  // float   (public)
	s>>c->z;  // float   (public)
	s>>c->t;  // float   (public)
	s>>c->dedx;  // float   (public)
	s>>c->nhits;  // int   (public)
	s.ReadArray(c->hits, 16);  // unsigned int[16]   (public)
	s>>c->chisq;  // float   (public)

	return s;
}

//----------------- DTrack -----------------
JILStream& operator<<(JILStream &s, const DTrack &c){

	// Send object type to stream
	s<<&typeid(DTrack);


	// Write out base class first
	s<<*((DObject*)&c);

	s<<c.q;  // float   (public)
	s<<c.p;  // float   (public)
	s<<c.theta;  // float   (public)
	s<<c.phi;  // float   (public)
	s<<c.x;  // float   (public)
	s<<c.y;  // float   (public)
	s<<c.z;  // float   (public)
	s<<c.candidateid;  // identifier_t   (public)

	// Send end of object flag to stream
	s<<JILStream::END_OBJECT;

	return s;
}

JILStream& operator>>(JILStream &s, DTrack* &c){
	// Create object if it doesn't already exist
	if(c == NULL)c = new DTrack();

	// Read in base class first
	s>>(DObject*&)c;

	s>>c->q;  // float   (public)
	s>>c->p;  // float   (public)
	s>>c->theta;  // float   (public)
	s>>c->phi;  // float   (public)
	s>>c->x;  // float   (public)
	s>>c->y;  // float   (public)
	s>>c->z;  // float   (public)
	s>>c->candidateid;  // identifier_t   (public)

	return s;
}

//----------------- DTrackCandidate -----------------
JILStream& operator<<(JILStream &s, const DTrackCandidate &c){

	// Send object type to stream
	s<<&typeid(DTrackCandidate);


	// Write out base class first
	s<<*((DObject*)&c);

	s<<c.hitid;  // vector   (public)
	s<<c.x0;  // float   (public)
	s<<c.y0;  // float   (public)
	s<<c.z_vertex;  // float   (public)
	s<<c.dphidz;  // float   (public)
	s<<c.q;  // float   (public)
	s<<c.p;  // float   (public)
	s<<c.p_trans;  // float   (public)
	s<<c.phi;  // float   (public)
	s<<c.theta;  // float   (public)

	// Send end of object flag to stream
	s<<JILStream::END_OBJECT;

	return s;
}

JILStream& operator>>(JILStream &s, DTrackCandidate* &c){
	// Create object if it doesn't already exist
	if(c == NULL)c = new DTrackCandidate();

	// Read in base class first
	s>>(DObject*&)c;

	s>>c->hitid;  // vector   (public)
	s>>c->x0;  // float   (public)
	s>>c->y0;  // float   (public)
	s>>c->z_vertex;  // float   (public)
	s>>c->dphidz;  // float   (public)
	s>>c->q;  // float   (public)
	s>>c->p;  // float   (public)
	s>>c->p_trans;  // float   (public)
	s>>c->phi;  // float   (public)
	s>>c->theta;  // float   (public)

	return s;
}

//----------------- DTrackEfficiency -----------------
JILStream& operator<<(JILStream &s, const DTrackEfficiency &c){

	// Send object type to stream
	s<<&typeid(DTrackEfficiency);


	// Write out base class first
	s<<*((DObject*)&c);

	s<<c.Nhits_thrown;  // int   (public)
	s<<c.Nhits_found;  // int   (public)
	s<<c.Nhits_thrown_and_found;  // int   (public)
	s<<c.Nhits_found_different;  // int   (public)
	s<<c.Nhits_thrown_unused;  // int   (public)
	s<<c.fittable;  // int   (public)
	s<<c.trackid;  // identifier_t   (public)

	// Send end of object flag to stream
	s<<JILStream::END_OBJECT;

	return s;
}

JILStream& operator>>(JILStream &s, DTrackEfficiency* &c){
	// Create object if it doesn't already exist
	if(c == NULL)c = new DTrackEfficiency();

	// Read in base class first
	s>>(DObject*&)c;

	s>>c->Nhits_thrown;  // int   (public)
	s>>c->Nhits_found;  // int   (public)
	s>>c->Nhits_thrown_and_found;  // int   (public)
	s>>c->Nhits_found_different;  // int   (public)
	s>>c->Nhits_thrown_unused;  // int   (public)
	s>>c->fittable;  // int   (public)
	s>>c->trackid;  // identifier_t   (public)

	return s;
}

//----------------- DTrackHit -----------------
JILStream& operator<<(JILStream &s, const DTrackHit &c){

	// Send object type to stream
	s<<&typeid(DTrackHit);


	// Write out base class first
	s<<*((DObject*)&c);

	s<<c.x;  // float   (public)
	s<<c.y;  // float   (public)
	s<<c.z;  // float   (public)
	s<<c.r;  // float   (public)
	s<<c.phi;  // float   (public)
	s<<(int)c.system;  // DetectorSystem_t   (public)

	// Send end of object flag to stream
	s<<JILStream::END_OBJECT;

	return s;
}

JILStream& operator>>(JILStream &s, DTrackHit* &c){
	// Create object if it doesn't already exist
	if(c == NULL)c = new DTrackHit();

	// Read in base class first
	s>>(DObject*&)c;

	s>>c->x;  // float   (public)
	s>>c->y;  // float   (public)
	s>>c->z;  // float   (public)
	s>>c->r;  // float   (public)
	s>>c->phi;  // float   (public)
	s>>*(int*)&c->system;  // DetectorSystem_t   (public)

	return s;
}

// -------------------- JILMyDictionary ----------------
const char* JILMyDictionary(void)
{
	return "<JILDictionary>\n"
"	<class name=\"DBCALHit\" baseclass=\"DObject\">\n"
"		<typedef type=\"int\" name=\"module\" section=\"public\"/>\n"
"		<typedef type=\"int\" name=\"layer\" section=\"public\"/>\n"
"		<typedef type=\"int\" name=\"sector\" section=\"public\"/>\n"
"		<typedef type=\"DBCALHit::END_t\" name=\"end\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"E\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"t\" section=\"public\"/>\n"
"	</class>\n"
"\n"
"	<class name=\"DCDCHit\" baseclass=\"DObject\">\n"
"		<typedef type=\"float\" name=\"radius\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"phim\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"dE\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"t\" section=\"public\"/>\n"
"	</class>\n"
"\n"
"	<class name=\"DFCALGeometry\" baseclass=\"DObject\">\n"
"		<typedef type=\"bool\" name=\"m_activeBlock\" section=\"private\" size=\"53,53\"/>\n"
"		<typedef type=\"TVector2\" name=\"m_positionOnFace\" section=\"private\" size=\",\"/>\n"
"		<typedef type=\"int\" name=\"m_channelNumber\" section=\"private\" size=\"53,53\"/>\n"
"		<typedef type=\"int\" name=\"m_row\" section=\"private\" size=\"2809\"/>\n"
"		<typedef type=\"int\" name=\"m_column\" section=\"private\" size=\"2809\"/>\n"
"		<typedef type=\"int\" name=\"m_numActiveBlocks\" section=\"private\"/>\n"
"	</class>\n"
"\n"
"	<class name=\"DFCALHit\" baseclass=\"DObject\">\n"
"		<typedef type=\"float\" name=\"x\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"y\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"E\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"t\" section=\"public\"/>\n"
"	</class>\n"
"\n"
"	<class name=\"DFCALMCResponse\" baseclass=\"DObject\">\n"
"		<typedef type=\"int\" name=\"m_channel\" section=\"private\"/>\n"
"		<typedef type=\"double\" name=\"m_E\" section=\"private\"/>\n"
"		<typedef type=\"double\" name=\"m_t\" section=\"private\"/>\n"
"	</class>\n"
"\n"
"	<class name=\"DFCALShower\" baseclass=\"DObject\">\n"
"		<typedef type=\"float\" name=\"x\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"y\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"E\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"t\" section=\"public\"/>\n"
"	</class>\n"
"\n"
"	<class name=\"DFDCHit\" baseclass=\"DObject\">\n"
"		<typedef type=\"int\" name=\"layer\" section=\"public\"/>\n"
"		<typedef type=\"int\" name=\"module\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"tau\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"z\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"u\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"dE\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"t\" section=\"public\"/>\n"
"		<typedef type=\"int\" name=\"type\" section=\"public\"/>\n"
"	</class>\n"
"\n"
"	<class name=\"DGeometry\">\n"
"		<typedef type=\"int\" name=\"min_run_number\" section=\"protected\" unsigned=\"true\"/>\n"
"		<typedef type=\"int\" name=\"max_run_number\" section=\"protected\" unsigned=\"true\"/>\n"
"	</class>\n"
"\n"
"	<class name=\"DHDDMForwardShower\" baseclass=\"DObject\">\n"
"		<typedef type=\"float\" name=\"m_x\" section=\"private\"/>\n"
"		<typedef type=\"float\" name=\"m_y\" section=\"private\"/>\n"
"		<typedef type=\"float\" name=\"m_E\" section=\"private\"/>\n"
"		<typedef type=\"float\" name=\"m_t\" section=\"private\"/>\n"
"	</class>\n"
"\n"
"	<class name=\"DHDDMTOFHit\" baseclass=\"DObject\">\n"
"		<typedef type=\"int\" name=\"orientation\" section=\"public\"/>\n"
"		<typedef type=\"int\" name=\"end\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"y\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"t\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"E\" section=\"public\"/>\n"
"	</class>\n"
"\n"
"	<class name=\"DHDDMTOFTruth\" baseclass=\"DObject\">\n"
"		<typedef type=\"int\" name=\"orientation\" section=\"public\"/>\n"
"		<typedef type=\"int\" name=\"track\" section=\"public\"/>\n"
"		<typedef type=\"int\" name=\"primary\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"x\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"y\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"z\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"t\" section=\"public\"/>\n"
"	</class>\n"
"\n"
"	<class name=\"DMCThrown\" baseclass=\"DObject\">\n"
"		<typedef type=\"int\" name=\"type\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"q\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"p\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"E\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"theta\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"phi\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"x\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"y\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"z\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"mass\" section=\"public\"/>\n"
"	</class>\n"
"\n"
"	<class name=\"DMCTrackHit\" baseclass=\"DObject\">\n"
"		<typedef type=\"float\" name=\"r\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"phi\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"z\" section=\"public\"/>\n"
"		<typedef type=\"int\" name=\"track\" section=\"public\"/>\n"
"		<typedef type=\"int\" name=\"primary\" section=\"public\"/>\n"
"		<typedef type=\"DetectorSystem_t\" name=\"system\" section=\"public\"/>\n"
"	</class>\n"
"\n"
"	<class name=\"DObject\">\n"
"		<typedef type=\"identifier_t\" name=\"id\" section=\"public\"/>\n"
"	</class>\n"
"\n"
"	<class name=\"DParameter\" hasSerializer=\"true\" hasDeserializer=\"true\">\n"
"		<typedef type=\"string\" name=\"key\" section=\"protected\"/>\n"
"		<typedef type=\"string\" name=\"value\" section=\"protected\"/>\n"
"		<typedef type=\"bool\" name=\"isdefault\" section=\"protected\"/>\n"
"		<typedef type=\"bool\" name=\"hasdefault\" section=\"protected\"/>\n"
"		<typedef type=\"bool\" name=\"printme\" section=\"protected\"/>\n"
"		<typedef type=\"DParameter::dataType_t\" name=\"type\" section=\"protected\"/>\n"
"	</class>\n"
"\n"
"	<class name=\"DQuickFit\">\n"
"		<typedef type=\"float\" name=\"x0\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"y0\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"q\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"p\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"p_trans\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"phi\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"theta\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"z_vertex\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"chisq\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"dphidz\" section=\"public\"/>\n"
"		<typedef type=\"DQuickFit::ChiSqSourceType_t\" name=\"chisq_source\" section=\"public\"/>\n"
"		<typedef type=\"vector\" name=\"hits\" section=\"protected\">\n"
"			<typedef type=\"DQFHit_t\" pointer=\"true\"/>\n"
"		</typedef>\n"
"		<typedef type=\"DMagneticFieldMap\" name=\"bfield\" section=\"protected\" pointer=\"true\"/>\n"
"		<typedef type=\"float\" name=\"Bz_avg\" section=\"protected\"/>\n"
"		<typedef type=\"float\" name=\"z_mean\" section=\"protected\"/>\n"
"		<typedef type=\"float\" name=\"phi_mean\" section=\"protected\"/>\n"
"	</class>\n"
"\n"
"	<class name=\"DTOFGeometry\" baseclass=\"DObject\">\n"
"		<typedef type=\"int\" name=\"NLONGBARS\" section=\"public\"/>\n"
"		<typedef type=\"int\" name=\"NSHORTBARS\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"LONGBARLENGTH\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"SHORTBARLENGTH\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"BARWIDTH\" section=\"public\"/>\n"
"	</class>\n"
"\n"
"	<class name=\"DTOFHit\" baseclass=\"DObject\">\n"
"		<typedef type=\"int\" name=\"orientation\" section=\"public\"/>\n"
"		<typedef type=\"int\" name=\"end\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"y\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"t\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"E\" section=\"public\"/>\n"
"	</class>\n"
"\n"
"	<class name=\"DTOFMCResponse\" baseclass=\"DObject\">\n"
"		<typedef type=\"int\" name=\"orientation\" section=\"public\"/>\n"
"		<typedef type=\"int\" name=\"end\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"y\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"t\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"E\" section=\"public\"/>\n"
"	</class>\n"
"\n"
"	<class name=\"DTOFPoint\" baseclass=\"DObject\">\n"
"		<typedef type=\"int\" name=\"trackid\" section=\"public\" unsigned=\"true\"/>\n"
"		<typedef type=\"float\" name=\"x\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"y\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"z\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"t\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"dedx\" section=\"public\"/>\n"
"		<typedef type=\"int\" name=\"nhits\" section=\"public\" unsigned=\"true\"/>\n"
"		<typedef type=\"unsigned int\" name=\"hits\" section=\"public\" size=\"16\"/>\n"
"		<typedef type=\"float\" name=\"chisq\" section=\"public\"/>\n"
"	</class>\n"
"\n"
"	<class name=\"DTrack\" baseclass=\"DObject\">\n"
"		<typedef type=\"float\" name=\"q\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"p\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"theta\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"phi\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"x\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"y\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"z\" section=\"public\"/>\n"
"		<typedef type=\"identifier_t\" name=\"candidateid\" section=\"public\"/>\n"
"	</class>\n"
"\n"
"	<class name=\"DTrackCandidate\" baseclass=\"DObject\">\n"
"		<typedef type=\"vector\" name=\"hitid\" section=\"public\">\n"
"			<typedef type=\"identifier_t\"/>\n"
"		</typedef>\n"
"		<typedef type=\"float\" name=\"x0\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"y0\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"z_vertex\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"dphidz\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"q\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"p\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"p_trans\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"phi\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"theta\" section=\"public\"/>\n"
"	</class>\n"
"\n"
"	<class name=\"DTrackEfficiency\" baseclass=\"DObject\">\n"
"		<typedef type=\"int\" name=\"Nhits_thrown\" section=\"public\"/>\n"
"		<typedef type=\"int\" name=\"Nhits_found\" section=\"public\"/>\n"
"		<typedef type=\"int\" name=\"Nhits_thrown_and_found\" section=\"public\"/>\n"
"		<typedef type=\"int\" name=\"Nhits_found_different\" section=\"public\"/>\n"
"		<typedef type=\"int\" name=\"Nhits_thrown_unused\" section=\"public\"/>\n"
"		<typedef type=\"int\" name=\"fittable\" section=\"public\"/>\n"
"		<typedef type=\"identifier_t\" name=\"trackid\" section=\"public\"/>\n"
"	</class>\n"
"\n"
"	<class name=\"DTrackHit\" baseclass=\"DObject\">\n"
"		<typedef type=\"float\" name=\"x\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"y\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"z\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"r\" section=\"public\"/>\n"
"		<typedef type=\"float\" name=\"phi\" section=\"public\"/>\n"
"		<typedef type=\"DetectorSystem_t\" name=\"system\" section=\"public\"/>\n"
"	</class>\n"
"\n"
"	<enum name=\"DetectorSystem_t\"/>\n"
"</JILDictionary>\n"
;
}

