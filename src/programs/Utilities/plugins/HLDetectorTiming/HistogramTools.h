#ifndef __HistogramTools__
#define __HistogramTools__

#include <iostream>
#include <stdexcept>

#include <JANA/JApplication.h>
#include <TH1I.h>
#include <TH2I.h>
#include <TDirectory.h>

using namespace std;
using namespace jana;

#ifndef ansi_escape
#define ansi_escape         ((char)0x1b)
#define ansi_bold           ansi_escape<<"[1m"
#define ansi_italic             ansi_escape<<"[3m"
#define ansi_underline      ansi_escape<<"[4m"
#define ansi_blink          ansi_escape<<"[5m"
#define ansi_rapid_blink    ansi_escape<<"[6m"
#define ansi_reverse            ansi_escape<<"[7m"
#define ansi_black          ansi_escape<<"[30m"
#define ansi_red                ansi_escape<<"[31m"
#define ansi_green          ansi_escape<<"[32m"
#define ansi_blue               ansi_escape<<"[34m"
#define ansi_normal         ansi_escape<<"[0m"
#define ansi_up(A)          ansi_escape<<"["<<(A)<<"A"
#define ansi_down(A)            ansi_escape<<"["<<(A)<<"B"
#define ansi_forward(A)     ansi_escape<<"["<<(A)<<"C"
#define ansi_back(A)            ansi_escape<<"["<<(A)<<"D"
#endif // ansi_escape

vector <TDirectory *>& GetAllDirectories(void);
map<TString, TH1I*>& Get1DMap(void);
map<TString, TH2I*>& Get2DMap(void);

void Fill1DHistogram (const char * plugin, const char * directoryName, const char * name, const double value = 0.0, const char * title = "", int nBins = 1, double xmin = 0, double xmax = 1, bool print = false);

void Fill2DHistogram (const char * plugin, const char * directoryName, const char * name, const double valueX = 0.0, const double valueY = 0.0, const char * title = "", int nBinsX = 1, double xmin = 0, double xmax = 1, int nBinsY = 1, double ymin = 0, double ymax = 1, bool print = false);

TH1I * Get1DHistogram(const char *, const char *, const char *);

TH2I * Get2DHistogram(const char *, const char *, const char *);

void SortDirectories(void);

#endif

