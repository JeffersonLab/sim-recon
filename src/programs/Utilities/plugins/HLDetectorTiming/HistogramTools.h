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

// Where all the histograms are sitting, global, but locked on access
static map<TString, TH1I*> TH1IMap;
static map<TString, TH2I*> TH2IMap;

static vector <TDirectory *> allDirectories; //Used to sort the directories later alphabetically 

void Fill1DHistogram (const char * plugin, const char * directoryName, const char * name, const double value = 0.0, const char * title = "", int nBins = 1, double xmin = 0, double xmax = 1, bool print = true){
    japp->RootWriteLock();
    TH1I * histogram;
    TString fullName = TString(plugin) + "/" + TString(directoryName) + "/" + TString(name);
    try {
        histogram = TH1IMap.at(fullName);
    }
    catch (const std::out_of_range& oor) {
        if (print) std::cerr << ansi_green << plugin << " ===> Making New 1D Histogram " << name << ansi_normal << endl;
        TDirectory *homedir = gDirectory;
        TDirectory *temp;
        temp = gDirectory->mkdir(plugin);
        if (temp) allDirectories.push_back(temp);
        gDirectory->cd(plugin);
        allDirectories.push_back(gDirectory->mkdir(directoryName));
        gDirectory->cd(directoryName);
        TH1IMap[fullName] = new TH1I( name, title, nBins, xmin, xmax);
        TH1IMap[fullName]->Fill(value);
        homedir->cd();
        japp->RootUnLock();
        return;
    }
    histogram->Fill(value);
    japp->RootUnLock();
    return;
}

void Fill2DHistogram (const char * plugin, const char * directoryName, const char * name, const double valueX = 0.0, const double valueY = 0.0, const char * title = "", int nBinsX = 1, double xmin = 0, double xmax = 1, int nBinsY = 1, double ymin = 0, double ymax = 1, bool print = true){
    japp->RootWriteLock();
    TH2I * histogram;
    TString fullName = TString(plugin) + "/" + TString(directoryName) + "/" + TString(name);
    try {
        histogram = TH2IMap.at(fullName);
    }
    catch (const std::out_of_range& oor) {
        if (print) std::cerr << ansi_green << plugin << " ===> Making New 2D Histogram " << name << ansi_normal << endl;
        TDirectory *homedir = gDirectory;
        TDirectory *temp;
        temp = gDirectory->mkdir(plugin);
        if (temp) allDirectories.push_back(temp);
        gDirectory->cd(plugin);
        allDirectories.push_back(gDirectory->mkdir(directoryName));
        gDirectory->cd(directoryName);
        TH2IMap[fullName] = new TH2I( name, title, nBinsX, xmin, xmax, nBinsY, ymin, ymax);
        TH2IMap[fullName]->Fill(valueX, valueY);
        homedir->cd();
        japp->RootUnLock();
        return;
    }
    histogram->Fill(valueX, valueY);
    japp->RootUnLock();
    return;
}

TH1I * Get1DHistogram(const char * plugin, const char * directoryName, const char * name){
    TH1I * histogram;
    TString fullName = TString(plugin) + "/" + TString(directoryName) + "/" + TString(name);
    japp->RootWriteLock(); // Lock this thread down while we search for the histogram
    try {
        histogram = TH1IMap.at(fullName);
    }
    catch (const std::out_of_range& oor) {
        cout << ansi_red << ansi_blink << "Specified histogram " << fullName.Data() << " does not exist" << ansi_normal << endl;
        // Let's be nice and find the closest match
        TString matchName = "";
        int closestMatch = 1000000; //Just initialize high
        map<TString, TH1I*>::const_iterator iter;
        for(iter = TH1IMap.begin(); iter != TH1IMap.end(); iter++){
            int match = fullName.CompareTo((*iter).first);
            if( match < closestMatch ){
                matchName = (*iter).first;
                closestMatch = match;
            }
        }
        cout << ansi_red << "The closest match is " << matchName.Data() << ansi_normal << endl;
        japp->RootUnLock();
        return NULL;
    }
    japp->RootUnLock();
    return histogram;
}

TH2I * Get2DHistogram(const char * plugin, const char * directoryName, const char * name){
    TH2I * histogram;
    TString fullName = TString(plugin) + "/" + TString(directoryName) + "/" + TString(name);
    japp->RootWriteLock(); // Lock this thread down while we search for the histogram
    try {
        histogram = TH2IMap.at(fullName);
    }
    catch (const std::out_of_range& oor) {
        cout << ansi_red << ansi_blink << "Specified histogram " << fullName.Data() << " does not exist" << ansi_normal << endl;
        // Let's be nice and find the closest match
        TString matchName = "";
        int closestMatch = 1000000; //Just initialize high
        map<TString, TH2I*>::const_iterator iter;
        for(iter = TH2IMap.begin(); iter != TH2IMap.end(); iter++){
            int match = fullName.CompareTo((*iter).first);
            if( match < closestMatch ){
                matchName = (*iter).first;
                closestMatch = match;
            }
        }
        cout << ansi_red << "The closest match is " << matchName.Data() << ansi_normal << endl;
        japp->RootUnLock();
        return NULL;
    }
    japp->RootUnLock();
    return histogram;
}

void SortDirectories(){
    for (unsigned int i=0; i < allDirectories.size(); i++){
        if (allDirectories[i] == 0) continue;
        japp->RootWriteLock();
        allDirectories[i]->GetList()->Sort();
        japp->RootUnLock();
    }
}
#endif

