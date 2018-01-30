#include "DLumi_factory.h"
#include "DLumi.h"

//------------------
// brun
//------------------
jerror_t DLumi_factory::brun(JEventLoop *loop, int32_t runnumber)
{

  std::cout << "DLumi:  BRUN " << std::endl;

  if(!_data.empty())
    {
      //for change in run #
      delete _data[0];
      _data.clear();
    }

  flags = PERSISTANT;
  _data.push_back( new DLumi(loop) );
   
  return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DLumi_factory::erun(void)
{
  for (unsigned int i = 0; i < _data.size(); i++)
    delete _data[i];
  _data.clear();
   
  return NOERROR;
}
