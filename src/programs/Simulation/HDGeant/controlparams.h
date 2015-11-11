
// This needs to be kept in sync with controlparams.inc!!

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	int writenohits;
	int shoersincol;
	int driftclusters;
	float tgwidth[2];
	int runtime_geom;
}controlparams_t;
extern controlparams_t controlparams_;

#ifdef __cplusplus
} // extern "C"
#endif
