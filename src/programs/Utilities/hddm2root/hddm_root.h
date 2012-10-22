#include<vector>
using namespace std;

#include<RTypes.h>
#include<TObject.h>

typedef int Particle_t;

class Cerenkov_t;
class ComptonEMcal_t;
class HDDM_t;
class TrackingErrorMatrix_t;
class barrelEMcal_t;
class bcalCell_t;
class bcalHit_t;
class bcalIncidentParticle_t;
class bcalSiPMDownHit_t;
class bcalSiPMSpectrum_t;
class bcalSiPMUpHit_t;
class bcalTDCHit_t;
class bcalTruthShower_t;
class bcalfADCCell_t;
class bcalfADCDownHit_t;
class bcalfADCUpHit_t;
class beam_t;
class ccalBlock_t;
class ccalHit_t;
class ccalTruthHit_t;
class ccalTruthShower_t;
class cdcStrawHit_t;
class cdcStrawTruthHit_t;
class cdcStraw_t;
class cdcTruthPoint_t;
class centralDC_t;
class cereHit_t;
class cereRichHit_t;
class cereSection_t;
class cereTruthPoint_t;
class errorMatrix_t;
class fcalBlock_t;
class fcalHit_t;
class fcalTruthHit_t;
class fcalTruthShower_t;
class fdcAnodeHit_t;
class fdcAnodeTruthHit_t;
class fdcAnodeWire_t;
class fdcCathodeHit_t;
class fdcCathodeStrip_t;
class fdcCathodeTruthHit_t;
class fdcChamber_t;
class fdcTruthPoint_t;
class forwardDC_t;
class forwardEMcal_t;
class forwardTOF_t;
class ftofCounter_t;
class ftofMCHit_t;
class ftofNorthHit_t;
class ftofNorthTruthHit_t;
class ftofSouthHit_t;
class ftofSouthTruthHit_t;
class ftofTruthPoint_t;
class gapEMcal_t;
class gcalCell_t;
class gcalHit_t;
class gcalTruthShower_t;
class hitView_t;
class mcTrajectoryPoint_t;
class mcTrajectory_t;
class microChannel_t;
class momentum_t;
class origin_t;
class physicsEvent_t;
class product_t;
class properties_t;
class random_t;
class reaction_t;
class reconView_t;
class startCntr_t;
class stcHit_t;
class stcPaddle_t;
class stcTruthHit_t;
class stcTruthPoint_t;
class taggerHit_t;
class tagger_t;
class target_t;
class tracktimebased_t;
class upstreamEMveto_t;
class upvLeftHit_t;
class upvPaddle_t;
class upvRightHit_t;
class upvTruthShower_t;
class vertex_t;

class fdcAnodeHit_t:TObject{
	public:
		float d;
		float dE;
		int itrack;
		int ptype;
		float t;

		ClassDef(fdcAnodeHit_t,1)
};

class fdcAnodeTruthHit_t:TObject{
	public:
		float d;
		float dE;
		int itrack;
		int ptype;
		float t;
		float t_unsmeared;

		ClassDef(fdcAnodeTruthHit_t,1)
};

class fdcCathodeHit_t:TObject{
	public:
		int itrack;
		int ptype;
		float q;
		float t;

		ClassDef(fdcCathodeHit_t,1)
};

class fdcCathodeTruthHit_t:TObject{
	public:
		int itrack;
		int ptype;
		float q;
		float t;

		ClassDef(fdcCathodeTruthHit_t,1)
};

class ftofMCHit_t:TObject{
	public:
		float E;
		float dist;
		int itrack;
		int ptype;
		float px;
		float py;
		float pz;
		float x;
		float y;
		float z;

		ClassDef(ftofMCHit_t,1)
};

class bcalHit_t:TObject{
	public:
		float E;
		float t;
		float zLocal;

		ClassDef(bcalHit_t,1)
};

class bcalSiPMDownHit_t:TObject{
	public:
		float E;
		float t;

		ClassDef(bcalSiPMDownHit_t,1)
};

class bcalSiPMUpHit_t:TObject{
	public:
		float E;
		float t;

		ClassDef(bcalSiPMUpHit_t,1)
};

class bcalfADCDownHit_t:TObject{
	public:
		float E;
		float t;

		ClassDef(bcalfADCDownHit_t,1)
};

class bcalfADCUpHit_t:TObject{
	public:
		float E;
		float t;

		ClassDef(bcalfADCUpHit_t,1)
};

class ccalHit_t:TObject{
	public:
		float E;
		float t;

		ClassDef(ccalHit_t,1)
};

class ccalTruthHit_t:TObject{
	public:
		float E;
		float t;

		ClassDef(ccalTruthHit_t,1)
};

class cdcStrawHit_t:TObject{
	public:
		float d;
		float dE;
		int itrack;
		int ptype;
		float t;

		ClassDef(cdcStrawHit_t,1)
};

class cdcStrawTruthHit_t:TObject{
	public:
		float d;
		float dE;
		int itrack;
		int ptype;
		float t;

		ClassDef(cdcStrawTruthHit_t,1)
};

class cereHit_t:TObject{
	public:
		float pe;
		float t;

		ClassDef(cereHit_t,1)
};

class fcalHit_t:TObject{
	public:
		float E;
		float t;

		ClassDef(fcalHit_t,1)
};

class fcalTruthHit_t:TObject{
	public:
		float E;
		float t;

		ClassDef(fcalTruthHit_t,1)
};

class fdcAnodeWire_t:TObject{
	public:
		vector<fdcAnodeHit_t> fdcAnodeHits;
		vector<fdcAnodeTruthHit_t> fdcAnodeTruthHits;
		int wire;

		ClassDef(fdcAnodeWire_t,1)
};

class fdcCathodeStrip_t:TObject{
	public:
		vector<fdcCathodeHit_t> fdcCathodeHits;
		vector<fdcCathodeTruthHit_t> fdcCathodeTruthHits;
		int plane;
		int strip;

		ClassDef(fdcCathodeStrip_t,1)
};

class fdcTruthPoint_t:TObject{
	public:
		float E;
		float dEdx;
		float dradius;
		bool primary;
		int ptype;
		float px;
		float py;
		float pz;
		float t;
		int track;
		float x;
		float y;
		float z;

		ClassDef(fdcTruthPoint_t,1)
};

class ftofNorthHit_t:TObject{
	public:
		float dE;
		float t;

		ClassDef(ftofNorthHit_t,1)
};

class ftofNorthTruthHit_t:TObject{
	public:
		float dE;
		vector<ftofMCHit_t> ftofMCHits;
		float t;

		ClassDef(ftofNorthTruthHit_t,1)
};

class ftofSouthHit_t:TObject{
	public:
		float dE;
		float t;

		ClassDef(ftofSouthHit_t,1)
};

class ftofSouthTruthHit_t:TObject{
	public:
		float dE;
		vector<ftofMCHit_t> ftofMCHits;
		float t;

		ClassDef(ftofSouthTruthHit_t,1)
};

class gcalHit_t:TObject{
	public:
		float E;
		float t;
		float zLocal;

		ClassDef(gcalHit_t,1)
};

class stcHit_t:TObject{
	public:
		float dE;
		float t;

		ClassDef(stcHit_t,1)
};

class stcTruthHit_t:TObject{
	public:
		float dE;
		float t;

		ClassDef(stcTruthHit_t,1)
};

class taggerHit_t:TObject{
	public:
		float t;

		ClassDef(taggerHit_t,1)
};

class upvLeftHit_t:TObject{
	public:
		float E;
		float t;

		ClassDef(upvLeftHit_t,1)
};

class upvRightHit_t:TObject{
	public:
		float E;
		float t;

		ClassDef(upvRightHit_t,1)
};

class TrackingErrorMatrix_t:TObject{
	public:
		int Ncols;
		int Nrows;
		string type;
		string vals;

		ClassDef(TrackingErrorMatrix_t,1)
};

class bcalCell_t:TObject{
	public:
		vector<bcalHit_t> bcalHits;
		vector<bcalSiPMDownHit_t> bcalSiPMDownHits;
		vector<bcalSiPMUpHit_t> bcalSiPMUpHits;
		int layer;
		int module;
		int sector;

		ClassDef(bcalCell_t,1)
};

class bcalIncidentParticle_t:TObject{
	public:
		int id;
		int ptype;
		float px;
		float py;
		float pz;
		float x;
		float y;
		float z;

		ClassDef(bcalIncidentParticle_t,1)
};

class bcalSiPMSpectrum_t:TObject{
	public:
		float Etruth;
		float bin_width;
		int end;
		int incident_id;
		int layer;
		int module;
		int sector;
		float tstart;
		string vals;

		ClassDef(bcalSiPMSpectrum_t,1)
};

class bcalTDCHit_t:TObject{
	public:
		int end;
		int layer;
		int module;
		int sector;
		float t;

		ClassDef(bcalTDCHit_t,1)
};

class bcalTruthShower_t:TObject{
	public:
		float E;
		float phi;
		bool primary;
		int ptype;
		float px;
		float py;
		float pz;
		float r;
		float t;
		int track;
		float z;

		ClassDef(bcalTruthShower_t,1)
};

class bcalfADCCell_t:TObject{
	public:
		vector<bcalfADCDownHit_t> bcalfADCDownHits;
		vector<bcalfADCUpHit_t> bcalfADCUpHits;
		int layer;
		int module;
		int sector;

		ClassDef(bcalfADCCell_t,1)
};

class ccalBlock_t:TObject{
	public:
		vector<ccalHit_t> ccalHits;
		vector<ccalTruthHit_t> ccalTruthHits;
		int column;
		int row;

		ClassDef(ccalBlock_t,1)
};

class ccalTruthShower_t:TObject{
	public:
		float E;
		bool primary;
		int ptype;
		float px;
		float py;
		float pz;
		float t;
		int track;
		float x;
		float y;
		float z;

		ClassDef(ccalTruthShower_t,1)
};

class cdcStraw_t:TObject{
	public:
		vector<cdcStrawHit_t> cdcStrawHits;
		vector<cdcStrawTruthHit_t> cdcStrawTruthHits;
		int ring;
		int straw;

		ClassDef(cdcStraw_t,1)
};

class cdcTruthPoint_t:TObject{
	public:
		float dEdx;
		float dradius;
		float phi;
		bool primary;
		int ptype;
		float px;
		float py;
		float pz;
		float r;
		float t;
		int track;
		float z;

		ClassDef(cdcTruthPoint_t,1)
};

class cereRichHit_t:TObject{
	public:
		float t;
		float x;
		float y;
		float z;

		ClassDef(cereRichHit_t,1)
};

class cereSection_t:TObject{
	public:
		vector<cereHit_t> cereHits;
		int sector;

		ClassDef(cereSection_t,1)
};

class cereTruthPoint_t:TObject{
	public:
		float E;
		bool primary;
		int ptype;
		float px;
		float py;
		float pz;
		float t;
		int track;
		float x;
		float y;
		float z;

		ClassDef(cereTruthPoint_t,1)
};

class errorMatrix_t:TObject{
	public:
		int Ncols;
		int Nrows;
		string type;
		string vals;

		ClassDef(errorMatrix_t,1)
};

class fcalBlock_t:TObject{
	public:
		int column;
		vector<fcalHit_t> fcalHits;
		vector<fcalTruthHit_t> fcalTruthHits;
		int row;

		ClassDef(fcalBlock_t,1)
};

class fcalTruthShower_t:TObject{
	public:
		float E;
		bool primary;
		int ptype;
		float px;
		float py;
		float pz;
		float t;
		int track;
		float x;
		float y;
		float z;

		ClassDef(fcalTruthShower_t,1)
};

class fdcChamber_t:TObject{
	public:
		vector<fdcAnodeWire_t> fdcAnodeWires;
		vector<fdcCathodeStrip_t> fdcCathodeStrips;
		vector<fdcTruthPoint_t> fdcTruthPoints;
		int layer;
		int module;

		ClassDef(fdcChamber_t,1)
};

class ftofCounter_t:TObject{
	public:
		int bar;
		vector<ftofNorthHit_t> ftofNorthHits;
		vector<ftofNorthTruthHit_t> ftofNorthTruthHits;
		vector<ftofSouthHit_t> ftofSouthHits;
		vector<ftofSouthTruthHit_t> ftofSouthTruthHits;
		int plane;

		ClassDef(ftofCounter_t,1)
};

class ftofTruthPoint_t:TObject{
	public:
		float E;
		bool primary;
		int ptype;
		float px;
		float py;
		float pz;
		float t;
		int track;
		float x;
		float y;
		float z;

		ClassDef(ftofTruthPoint_t,1)
};

class gcalCell_t:TObject{
	public:
		vector<gcalHit_t> gcalHits;
		int module;

		ClassDef(gcalCell_t,1)
};

class gcalTruthShower_t:TObject{
	public:
		float E;
		float phi;
		bool primary;
		int ptype;
		float px;
		float py;
		float pz;
		float r;
		float t;
		int track;
		float z;

		ClassDef(gcalTruthShower_t,1)
};

class mcTrajectoryPoint_t:TObject{
	public:
		float E;
		float dE;
		int mech;
		int part;
		int primary_track;
		float px;
		float py;
		float pz;
		float radlen;
		float step;
		float t;
		int track;
		float x;
		float y;
		float z;

		ClassDef(mcTrajectoryPoint_t,1)
};

class microChannel_t:TObject{
	public:
		float E;
		int column;
		int row;
		vector<taggerHit_t> taggerHits;

		ClassDef(microChannel_t,1)
};

class momentum_t:TObject{
	public:
		float E;
		float px;
		float py;
		float pz;

		ClassDef(momentum_t,1)
};

class origin_t:TObject{
	public:
		float t;
		float vx;
		float vy;
		float vz;

		ClassDef(origin_t,1)
};

class product_t:TObject{
	public:
		int decayVertex;
		int id;
		int mech;
		momentum_t momentum;
		int parentid;
		int pdgtype;
		properties_t properties;
		Particle_t type;

		ClassDef(product_t,1)
};

class properties_t:TObject{
	public:
		int charge;
		float mass;

		ClassDef(properties_t,1)
};

class stcPaddle_t:TObject{
	public:
		int sector;
		vector<stcHit_t> stcHits;
		vector<stcTruthHit_t> stcTruthHits;

		ClassDef(stcPaddle_t,1)
};

class stcTruthPoint_t:TObject{
	public:
		float E;
		float dEdx;
		float phi;
		bool primary;
		int ptype;
		float px;
		float py;
		float pz;
		float r;
		int sector;
		float t;
		int track;
		float z;

		ClassDef(stcTruthPoint_t,1)
};

class upvPaddle_t:TObject{
	public:
		int layer;
		int row;
		vector<upvLeftHit_t> upvLeftHits;
		vector<upvRightHit_t> upvRightHits;

		ClassDef(upvPaddle_t,1)
};

class upvTruthShower_t:TObject{
	public:
		float E;
		bool primary;
		int ptype;
		float px;
		float py;
		float pz;
		float t;
		int track;
		float x;
		float y;
		float z;

		ClassDef(upvTruthShower_t,1)
};

class Cerenkov_t:TObject{
	public:
		vector<cereRichHit_t> cereRichHits;
		vector<cereSection_t> cereSections;
		vector<cereTruthPoint_t> cereTruthPoints;

		ClassDef(Cerenkov_t,1)
};

class ComptonEMcal_t:TObject{
	public:
		vector<ccalBlock_t> ccalBlocks;
		vector<ccalTruthShower_t> ccalTruthShowers;

		ClassDef(ComptonEMcal_t,1)
};

class barrelEMcal_t:TObject{
	public:
		vector<bcalCell_t> bcalCells;
		vector<bcalIncidentParticle_t> bcalIncidentParticles;
		vector<bcalSiPMSpectrum_t> bcalSiPMSpectrums;
		vector<bcalTDCHit_t> bcalTDCHits;
		vector<bcalTruthShower_t> bcalTruthShowers;
		vector<bcalfADCCell_t> bcalfADCCells;

		ClassDef(barrelEMcal_t,1)
};

class beam_t:TObject{
	public:
		momentum_t momentum;
		properties_t properties;
		Particle_t type;

		ClassDef(beam_t,1)
};

class centralDC_t:TObject{
	public:
		vector<cdcStraw_t> cdcStraws;
		vector<cdcTruthPoint_t> cdcTruthPoints;

		ClassDef(centralDC_t,1)
};

class forwardDC_t:TObject{
	public:
		vector<fdcChamber_t> fdcChambers;

		ClassDef(forwardDC_t,1)
};

class forwardEMcal_t:TObject{
	public:
		vector<fcalBlock_t> fcalBlocks;
		vector<fcalTruthShower_t> fcalTruthShowers;

		ClassDef(forwardEMcal_t,1)
};

class forwardTOF_t:TObject{
	public:
		vector<ftofCounter_t> ftofCounters;
		vector<ftofTruthPoint_t> ftofTruthPoints;

		ClassDef(forwardTOF_t,1)
};

class gapEMcal_t:TObject{
	public:
		vector<gcalCell_t> gcalCells;
		vector<gcalTruthShower_t> gcalTruthShowers;

		ClassDef(gapEMcal_t,1)
};

class mcTrajectory_t:TObject{
	public:
		vector<mcTrajectoryPoint_t> mcTrajectoryPoints;

		ClassDef(mcTrajectory_t,1)
};

class random_t:TObject{
	public:
		int seed1;
		int seed2;
		int seed_mcsmear1;
		int seed_mcsmear2;
		int seed_mcsmear3;

		ClassDef(random_t,1)
};

class startCntr_t:TObject{
	public:
		vector<stcPaddle_t> stcPaddles;
		vector<stcTruthPoint_t> stcTruthPoints;

		ClassDef(startCntr_t,1)
};

class tagger_t:TObject{
	public:
		vector<microChannel_t> microChannels;

		ClassDef(tagger_t,1)
};

class target_t:TObject{
	public:
		momentum_t momentum;
		properties_t properties;
		Particle_t type;

		ClassDef(target_t,1)
};

class tracktimebased_t:TObject{
	public:
		float FOM;
		int Ndof;
		TrackingErrorMatrix_t TrackingErrorMatrix;
		int candidateid;
		float chisq;
		errorMatrix_t errorMatrix;
		int id;
		momentum_t momentum;
		origin_t origin;
		properties_t properties;
		int trackid;

		ClassDef(tracktimebased_t,1)
};

class upstreamEMveto_t:TObject{
	public:
		vector<upvPaddle_t> upvPaddles;
		vector<upvTruthShower_t> upvTruthShowers;

		ClassDef(upstreamEMveto_t,1)
};

class vertex_t:TObject{
	public:
		origin_t origin;
		vector<product_t> products;

		ClassDef(vertex_t,1)
};

class hitView_t:TObject{
	public:
		Cerenkov_t Cerenkov;
		ComptonEMcal_t ComptonEMcal;
		barrelEMcal_t barrelEMcal;
		centralDC_t centralDC;
		forwardDC_t forwardDC;
		forwardEMcal_t forwardEMcal;
		forwardTOF_t forwardTOF;
		gapEMcal_t gapEMcal;
		mcTrajectory_t mcTrajectory;
		startCntr_t startCntr;
		tagger_t tagger;
		upstreamEMveto_t upstreamEMveto;

		ClassDef(hitView_t,1)
};

class reaction_t:TObject{
	public:
		beam_t beam;
		random_t random;
		target_t target;
		int type;
		vector<vertex_t> vertexs;
		float weight;

		ClassDef(reaction_t,1)
};

class reconView_t:TObject{
	public:
		vector<tracktimebased_t> tracktimebaseds;

		ClassDef(reconView_t,1)
};

class physicsEvent_t:TObject{
	public:
		int eventNo;
		hitView_t hitView;
		vector<reaction_t> reactions;
		reconView_t reconView;
		int runNo;

		ClassDef(physicsEvent_t,1)
};

