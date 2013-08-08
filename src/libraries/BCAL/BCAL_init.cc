// $Id$

#include <JANA/JEventLoop.h>
#include "DBCALGeometry_factory.h"
#include "DBCALShower_factory_IU.h"
#include "DBCALShower_factory_KLOE.h"
#include "DBCALShower_factory_KLOE_OLDSMEAR.h"
#include "DBCALShower_factory.h"
#include "DBCALCluster_factory.h"
#include "DBCALCluster_factory_SINGLE.h"
#include "DBCALPoint_factory.h"
#include "DBCALPoint_factory_OLDSMEAR.h"
#include "DBCALUnifiedHit_factory.h"
#include "DBCALDigiHit.h"
#include "DBCALHit_factory.h"
#include "DBCALIncidentParticle.h"
#include "DBCALTDCDigiHit.h"
#include "DBCALTDCHit_factory.h"
#include "DBCALSiPMHit.h"
#include "DBCALSiPMSpectrum.h"
#include "DBCALTruthCell.h"
#include "DBCALClump.h"
#include "DBCALClump_factory.h"
#include "DBCALShower_factory_JLAB.h"

#include "DBCALTruthShower.h"

// These come from the event source, not from any algorithm
typedef JFactory<DBCALDigiHit> DBCALDigiHit_factory;
typedef JFactory<DBCALTDCDigiHit> DBCALTDCDigiHit_factory;
typedef JFactory<DBCALIncidentParticle> DBCALIncidentParticle_factory;
typedef JFactory<DBCALSiPMHit> DBCALSiPMHit_factory;
typedef JFactory<DBCALSiPMSpectrum> DBCALSiPMSpectrum_factory;
typedef JFactory<DBCALTruthShower> DBCALTruthShower_factory;
typedef JFactory<DBCALTruthCell> DBCALTruthCell_factory;

jerror_t BCAL_init(JEventLoop *loop)
{
	/// Create and register BCAL data factories
	loop->AddFactory(new DBCALDigiHit_factory());
	loop->AddFactory(new DBCALTDCDigiHit_factory());
	loop->AddFactory(new DBCALHit_factory());
	loop->AddFactory(new DBCALIncidentParticle_factory());
	loop->AddFactory(new DBCALTDCHit_factory());
	loop->AddFactory(new DBCALSiPMHit_factory());
	loop->AddFactory(new DBCALSiPMSpectrum_factory());
	loop->AddFactory(new DBCALGeometry_factory());
	loop->AddFactory(new DBCALShower_factory_IU());
	loop->AddFactory(new DBCALShower_factory_KLOE());
	loop->AddFactory(new DBCALShower_factory_KLOE_OLDSMEAR());
	loop->AddFactory(new DBCALShower_factory());
	loop->AddFactory(new DBCALCluster_factory());
	loop->AddFactory(new DBCALCluster_factory_SINGLE());
	loop->AddFactory(new DBCALTruthShower_factory());
	loop->AddFactory(new DBCALTruthCell_factory());
	loop->AddFactory(new DBCALPoint_factory());
	loop->AddFactory(new DBCALPoint_factory_OLDSMEAR());
	loop->AddFactory(new DBCALUnifiedHit_factory());
        loop->AddFactory(new DBCALClump_factory());
        loop->AddFactory(new DBCALShower_factory_JLAB());
   
	return NOERROR;
}
