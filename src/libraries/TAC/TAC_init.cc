/*
 * TAC_init.cc
 *
 *  Created on: Mar 28, 2017
 *      Author: Hovanes Egiyan
 */

#include <JANA/JEventLoop.h>
using namespace jana;

#include "DTACDigiHit.h"
#include "DTACTDCDigiHit.h"
#include "DTACHit_factory.h"
#include "DRebuildFromRawFADC_factory.h"
#include "HitRebuilderByFit.h"
#include "WaveformSpikeFunctor.h"
#include "WaveformErfcFunctor.h"

jerror_t TAC_init(JEventLoop *loop) {
	/// Create and register TAC data factories
	loop->AddFactory(new JFactory<DTACDigiHit>());
	loop->AddFactory(new JFactory<DTACTDCDigiHit>());
	loop->AddFactory(new DTACHit_factory());
	loop->AddFactory(
			new DRebuildFromRawFADC_factory<DTACHit_factory,
					HitRebuilderByFit<WaveformSpikeFunctor>>());
	loop->AddFactory(
			new DRebuildFromRawFADC_factory<DTACHit_factory,
					HitRebuilderByFit<WaveformErfcFunctor>>());

	loop->AddFactory(new DRebuildFromRawFADC_factory<>());
	loop->AddFactory(new JFactory<DTACHit>("TRUTH"));

	return NOERROR;
}

