// TOF Geometry class
#include <DTOFGeometry.h>
#include <map>
#include <vector>

// Constructor
DTOFGeometry::DTOFGeometry(const DGeometry* locGeometry) {
	vector<double> tmp;
	locGeometry->GetTOFZ(tmp);

  	// Store the z position for both planes
  	vector<double>tof_face;
  	locGeometry->Get("//section/composition/posXYZ[@volume='ForwardTOF']/@X_Y_Z",tof_face);
  	vector<double>tof_plane0;
  	locGeometry->Get("//composition[@name='ForwardTOF']/posXYZ[@volume='forwardTOF']/@X_Y_Z/plane[@value='0']", tof_plane0);
  	vector<double>tof_plane1;
  	locGeometry->Get("//composition[@name='ForwardTOF']/posXYZ[@volume='forwardTOF']/@X_Y_Z/plane[@value='1']", tof_plane1);
  	CenterVPlane=tof_face[2]+tof_plane1[2];
  	CenterHPlane=tof_face[2]+tof_plane0[2];
  	// also save position midway between the two planes
  	CenterMPlane=0.5*(CenterHPlane+CenterVPlane);
	
	// load more geometry parameters 
	if(!locGeometry->GetTOFPaddlePerpPositions(YPOS))
		jerr << "Problem loading TOF Geometry!" << endl;

	map<string,double> paddle_params;
	if(!locGeometry->GetTOFPaddleParameters(paddle_params))
		jerr << "Problem loading TOF Geometry!" << endl;
	else {
		NLONGBARS = paddle_params["NLONGBARS"];
		NSHORTBARS = paddle_params["NSHORTBARS"];
		BARWIDTH = paddle_params["BARWIDTH"];

		LONGBARLENGTH = paddle_params["LONGBARLENGTH"];
		HALFLONGBARLENGTH = paddle_params["HALFLONGBARLENGTH"];
		SHORTBARLENGTH = paddle_params["SHORTBARLENGTH"];
		HALFSHORTBARLENGTH = paddle_params["HALFSHORTBARLENGTH"];

		FirstShortBar = paddle_params["FIRSTSHORTBAR"];
		LastShortBar = paddle_params["LASTSHORTBAR"];
	}

    // save the number of bars that exist along one side of the detector
	NINSTALLBARS = NLONGBARS + NSHORTBARS/2;

	// hardcode these for now, they are not likely to change anytime soon
  	NLAYERS        = 2;
  	NENDS          = 2;

}

float DTOFGeometry::bar2y(int bar, int end)  const 
///> convert bar number to the
///> position of the center of the
///> bar in local coordinations
{
    float y;
    y = YPOS.at(bar);

    // handle position of short bars
    if (bar>=FirstShortBar && bar<=LastShortBar && end != 0) y *= -1.0;

    return y;
}
  
  
int DTOFGeometry::y2bar(double y) const   ///> convert local position y to bar number
///> (where y is the position perpendicular to the bar length)
{
	if(y < (YPOS.at(1) - BARWIDTH/2.0))
		return 0;
	if(y > (YPOS.at(NINSTALLBARS) + BARWIDTH/2.0))
		return 0;

	int jm=1;
	if (y>YPOS.at(1))
	{
		int jl=-1;
		int ju=NINSTALLBARS;
		while(ju - jl > 1)
		{
			jm = (ju + jl) >> 1;
			if (y >= YPOS.at(jm))
				jl = jm;
			else
				ju = jm;
		}
		if (fabs(y - YPOS.at(jm - 1)) < fabs(y - YPOS.at(jm)))
			return jm - 1;
		if (fabs(y - YPOS.at(jm + 1)) < fabs(y - YPOS.at(jm)))
			return jm + 1;
	}

	return jm;
}
