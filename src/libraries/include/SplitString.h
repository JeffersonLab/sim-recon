#ifndef _SplitString_
#define _SplitString_

#include <sstream>

template<typename T>
void SplitString(string str, vector<T> &vec, const string &delim=" ")
{
	if(str.length()==0)return;
	unsigned int cutAt;
	while( (cutAt = str.find(delim)) != (unsigned int)str.npos ){
		if(cutAt > 0){
			std::stringstream ss(str.substr(0,cutAt));
			T t;
			ss>>t;
			vec.push_back(t);
		}
		str = str.substr(cutAt+1);
	}
	if(str.length() > 0){
		std::stringstream ss(str);
		T t;
		ss>>t;
		vec.push_back(t);
	}
}

#endif
