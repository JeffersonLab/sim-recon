// $Id$
//
//    File: DBeamPhoton_factory.cc
// Created: Thu Dec  3 17:27:55 EST 2009
// Creator: staylor (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#include <iostream>
#include <iomanip>
using namespace std;

#include "DBeamPhoton_factory.h"
using namespace jana;

//------------------
// init
//------------------
jerror_t DBeamPhoton_factory::init(void)
{
    DELTA_T_DOUBLES_MAX = 1.5; // ns
    gPARMS->SetDefaultParameter("BeamPhoton:DELTA_T_DOUBLES_MAX", DELTA_T_DOUBLES_MAX,
    "Maximum time difference in ns between a TAGM-TAGH pair of beam photons"
    " for them to be merged into a single photon");
    DELTA_E_DOUBLES_MAX = 0.05; // GeV
    gPARMS->SetDefaultParameter("BeamPhoton:DELTA_E_DOUBLES_MAX", DELTA_E_DOUBLES_MAX,
    "Maximum energy difference in GeV between a TAGM-TAGH pair of beam photons"
    " for them to be merged into a single photon");
    return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DBeamPhoton_factory::brun(jana::JEventLoop *locEventLoop, int32_t runnumber)
{
    DApplication* dapp = dynamic_cast<DApplication*>(locEventLoop->GetJApplication());
    DGeometry* locGeometry = dapp->GetDGeometry(locEventLoop->GetJEvent().GetRunNumber());
    dTargetCenterZ = 0.0;
    locGeometry->GetTargetZ(dTargetCenterZ);

    return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DBeamPhoton_factory::evnt(jana::JEventLoop *locEventLoop, uint64_t eventnumber)
{
    vector<const DTAGMHit*> tagm_hits;
    locEventLoop->Get(tagm_hits);

    for (unsigned int ih=0; ih < tagm_hits.size(); ++ih)
    {
        if (!tagm_hits[ih]->has_fADC) continue; // Skip TDC-only hits (i.e. hits with no ADC info.)
        if (tagm_hits[ih]->row > 0) continue; // Skip individual fiber readouts
        DBeamPhoton *gamma = new DBeamPhoton;
        Set_BeamPhoton(gamma, tagm_hits[ih]);
        _data.push_back(gamma);
    }

    vector<const DTAGHHit*> tagh_hits;
    locEventLoop->Get(tagh_hits);

    for (unsigned int ih=0; ih < tagh_hits.size(); ++ih)
    {
        if (!tagh_hits[ih]->has_fADC) continue; // Skip TDC-only hits (i.e. hits with no ADC info.)
        DBeamPhoton *gamma = nullptr;
        for (unsigned int jh=0; jh < _data.size(); ++jh)
        {
            if (fabs(_data[jh]->momentum().Mag() - tagh_hits[ih]->E) < DELTA_E_DOUBLES_MAX
            && fabs(_data[jh]->time() - tagh_hits[ih]->t) < DELTA_T_DOUBLES_MAX)
            {
                gamma = _data[jh];
                if (_data[jh]->momentum().Mag() < tagh_hits[ih]->E)
                {
                    gamma->Reset();
                    Set_BeamPhoton(gamma, tagh_hits[ih]);
                }
            }
        }
        if (gamma == nullptr)
        {
            gamma = new DBeamPhoton;
            Set_BeamPhoton(gamma, tagh_hits[ih]);
            _data.push_back(gamma);
        }
    }

    return NOERROR;
}

void DBeamPhoton_factory::Set_BeamPhoton(DBeamPhoton* gamma, const DTAGMHit* hit)
{
    DVector3 pos(0.0, 0.0, dTargetCenterZ);
    DVector3 mom(0.0, 0.0, hit->E);
    gamma->setPID(Gamma);
    gamma->setMomentum(mom);
    gamma->setPosition(pos);
    gamma->setCharge(0);
    gamma->setMass(0);
    gamma->setTime(hit->t);
    gamma->setT0(hit->t, 0.200, SYS_TAGM);
    gamma->dCounter = hit->column;
    gamma->AddAssociatedObject(hit);
}

void DBeamPhoton_factory::Set_BeamPhoton(DBeamPhoton* gamma, const DTAGHHit* hit)
{
    DVector3 pos(0.0, 0.0, dTargetCenterZ);
    DVector3 mom(0.0, 0.0, hit->E);
    gamma->setPID(Gamma);
    gamma->setMomentum(mom);
    gamma->setPosition(pos);
    gamma->setCharge(0);
    gamma->setMass(0);
    gamma->setTime(hit->t);
    gamma->setT0(hit->t, 0.350, SYS_TAGH);
    gamma->dCounter = hit->counter_id;
    gamma->AddAssociatedObject(hit);
}

//------------------
// erun
//------------------
jerror_t DBeamPhoton_factory::erun(void)
{
    return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DBeamPhoton_factory::fini(void)
{
    return NOERROR;
}

