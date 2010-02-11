#define NWGEAN 9000000
#define NWPAW 6000000

typedef struct{
	int GEANT[NWGEAN];
}GCBANK_t;

typedef struct{
	int PAW[NWPAW];
}PAWC_t;

typedef struct {
    float EINI;
	 float TPI0;
	 float FIPI0G;
	 float EPI0SP;
	 float TRECSP;
	 float TKINRM;
} kinem1_t;

typedef struct {
    float EPI0LF;
	 float PPI0LF[3];
	 float EG1LF;
	 float PG1LF[3];
	 float EG2LF;
	 float PG2LF[3];
} kinem3_t;

#ifdef __cplusplus
extern "C" {
#endif
	void gukine_(void);
	void gzebra_(int *nwgean);
	void hlimit_(int *nwpaw);
	void gpaw_(int *nwgean, int *nwpaw);
	void uginit_(void);
	
	extern GCBANK_t gcbank_;
	extern PAWC_t pawc_;
	extern int quest_[100];
	extern kinem1_t kinem1_;
	extern kinem3_t kinem3_;
#ifdef __cplusplus
}
#endif
