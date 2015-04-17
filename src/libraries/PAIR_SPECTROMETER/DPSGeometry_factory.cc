#include "DPSGeometry_factory.h"
#include "DPSGeometry.h"

//------------------
// brun
//------------------
jerror_t DPSGeometry_factory::brun(JEventLoop *loop, int runnumber)
{
  if(!_data.empty())
    {
      //for change in run #
      delete _data[0];
      _data.clear();
    }

  flags = PERSISTANT;
  _data.push_back( new DPSGeometry(loop, factory_tag, runnumber) );
   
  return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DPSGeometry_factory::erun(void)
{
  for (unsigned int i=0; i < _data.size(); i++)
    delete _data[i];
  _data.clear();
   
  return NOERROR;
}
