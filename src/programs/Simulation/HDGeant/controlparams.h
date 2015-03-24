
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
  int get_next_evt;
  float trigger_time_signa_ns;
  int runno_ff;
}controlparams_t;
extern controlparams_t controlparams_;

#ifdef __cplusplus
} // extern "C"
#endif
