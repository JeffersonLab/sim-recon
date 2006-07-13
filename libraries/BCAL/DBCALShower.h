#ifndef _DBCALShower_
#define _DBCALShower_

#include "JANA/JObject.h"
#include "JANA/JFactory.h"

class DBCALShower:public JObject{
	public:
		HDCLASSDEF(DBCALShower);

    float E;    
    float Ecorr;
    float x;
    float y;
    float z;   
    float t;
    int N_cell;
    int total_layer_cluster;
    float Apx_x;
    float Apx_y;
    float Apx_z;
    float error_Apx_x;
    float error_Apx_y;
    float error_Apx_z;
    float Cx;
    float Cy;
    float Cz;
    float error_Cx;
    float error_Cy;
    float error_Cz;
    float t_rms_a;
    float t_rms_b;

};

#endif // _DBCALShower_

