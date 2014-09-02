#include <string>

typedef struct{
        int geantid;
        int mech; /* what do the values of this correspond to */
        int kfid;
        int parent;
        int firstdaughter;
        int lastdaughter;
}keve_t;

typedef struct{
	float px;
	float py;
	float pz;
	float en;
}peve_t;


void open_hddm_output(std::string outputfile);
void close_hddm_output(void);
void write_hddm_event(int *iev, int *iproc, keve_t *kin,  peve_t *pin, int *ntra, keve_t *keve, peve_t *peve);

