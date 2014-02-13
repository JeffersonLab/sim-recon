/* 
 * prototypes for c_cern package
 *
 * clas c interface to hbook and hplot
 *
 *  jm, ejw, 8-jan-1998
 *
 */


#ifdef __cplusplus
extern "C" {
#endif


void hlimit(int size);
void hrput(int id, char *file, char *opt);
void hf1(int id,float data,float weight);
void hf2(int id,float data1,float data2,float weight);
void hdelet(int id);
void hrin( int id, int icycle , int iofset );
void hfn(int id, float data[]);
void hfithn(int id, char *chfun, char *chopt, int np, float par[],
        float step[], float pmin[], float pmax[], float sigpar[], float *chi2);
void hunpak(int histo,float contents[],char choice[], int num);
void hidopt( int id, char *chopt);
void hpak(int histo,float contents[]);
void hrget( int id, char *chfile, char *chopt);
void hldir(char dir[],char flag[]);
void hmdir(char dir[],char flag[]);
void hcdir(char dir[],char flag[]);
void hropen(int lun, char *name, char*filename, char*status, int stor,
        int istat);
void hrout(int num, int icycle, char*opt);
void hrend(char*filename);
void hreset(int no, char* opt);
void hbook2(int no, char*name, int xbins, float xmin, float xmax, int ybins,
                 float ymin, float ymax, float weight);
void hbook1(int no, char*name, int nbins, float min, float max, float v);
void hfill(int no, float xvalue, float yvalue, float weight);


/* hplot routines located in file clas_hplot.c - Keep these in a separate file!!!! */

void hplint(int no);
void hplego(int no, float theta, float phi);
void hplcon(int histonum, int x, int y);
void hplzon(int x, int y, int num, char *opt);
void hplot(int no, char*chopt, char*chcase, int num);
void iuwk(int num1, int num2);
void hplcap(int unit);
void hplzon(int nx, int ny, int fistplotted, char *option);


#ifdef __cplusplus
}
#endif
