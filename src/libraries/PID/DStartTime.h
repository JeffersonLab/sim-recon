#ifndef _DStartTime_
#define _DStartTime_

#include <string>
#include <utility>
#include <vector>

#include <JANA/JObject.h>
#include <JANA/JFactory.h>

class DStartTime : public jana::JObject {
 public:
  JOBJECT_PUBLIC(DStartTime);
  DStartTime(float start_time = 0.0);
  virtual ~DStartTime();

  float getStartTime();

  void toStrings(std::vector<std::pair<std::string, std::string> > &items) const;
  
 private:
  float _t_start;    // Event start time as determined by various means.
  
};

inline void DStartTime::toStrings(std::vector<std::pair<std::string, std::string> > &items) const {
  AddString(items, "start_time", "%f", _t_start);
}

inline float DStartTime::getStartTime() {
  return _t_start;
}



#endif // _DStartTime_
