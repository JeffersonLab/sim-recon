	typedef struct {
	   int col;	
	   int row;
	   float x;
	   float y;
	   float E;
	   float t;
	} userhit_t;

	typedef struct {
	   int nhits;
	   userhit_t hit[1];
	} userhits_t;

	typedef struct {
	   float x;
	   float y;
	   float z;
	} vector3_t;
		   
