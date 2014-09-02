

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
using namespace std;

#include <signal.h>
#include <time.h>
#include <stdlib.h>


extern vector<char*> INFILENAMES;
extern char *OUTFILENAME;
extern int QUIT;
extern bool HDDM_USE_COMPRESSION;
extern bool HDDM_USE_INTEGRITY_CHECKS;

#define _DBG_ cout<<__FILE__<<":"<<__LINE__<<" "
#define _DBG__ cout<<__FILE__<<":"<<__LINE__<<endl


void Process_s(unsigned int &NEvents, unsigned int &NEvents_read);
void Process_r(unsigned int &NEvents, unsigned int &NEvents_read);
