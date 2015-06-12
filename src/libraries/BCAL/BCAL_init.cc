// $Id$

#include <JANA/JEventLoop.h>
#include "DBCALGeometry_factory.h"
#include "DBCALShower_factory_IU.h"
#include "DBCALShower_factory_KLOE.h"
#include "DBCALShower_factory.h"
#include "DBCALCluster_factory.h"
#include "DBCALCluster_factory_SINGLE.h"
#include "DBCALPoint_factory.h"
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


jerror_t BCAL_init(JEventLoop *loop)
{
	/// Create and register BCAL data factories
	loop->AddFactory(new JFactory<DBCALDigiHit>());
	loop->AddFactory(new JFactory<DBCALTDCDigiHit>());
	loop->AddFactory(new DBCALHit_factory());
	loop->AddFactory(new JFactory<DBCALIncidentParticle>());
	loop->AddFactory(new DBCALTDCHit_factory());
	loop->AddFactory(new JFactory<DBCALSiPMHit>());
	loop->AddFactory(new JFactory<DBCALSiPMSpectrum>());
	loop->AddFactory(new JFactory<DBCALSiPMSpectrum>("TRUTH"));
	loop->AddFactory(new DBCALGeometry_factory());
	loop->AddFactory(new DBCALShower_factory_IU());
	loop->AddFactory(new DBCALShower_factory_KLOE());
	loop->AddFactory(new DBCALShower_factory());
	loop->AddFactory(new DBCALCluster_factory());
	loop->AddFactory(new DBCALCluster_factory_SINGLE());
	loop->AddFactory(new JFactory<DBCALTruthShower>());
	loop->AddFactory(new JFactory<DBCALTruthCell>());
	loop->AddFactory(new DBCALPoint_factory());
	loop->AddFactory(new DBCALUnifiedHit_factory());
	loop->AddFactory(new DBCALClump_factory());
	loop->AddFactory(new DBCALShower_factory_JLAB());

	return NOERROR;
}
