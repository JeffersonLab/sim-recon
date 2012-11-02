
// This needs to be kept in sync with controlparams.inc!!

extern "C" {
typedef struct {
	int writenohits;
	int shoersincol;
	int driftclusters;
	float tgwidth[2];
	int runtime_geom;
}controlparams_t;
extern controlparams_t controlparams_;

} // extern "C"

