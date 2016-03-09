#ifndef __HistogramTools__
#define __HistogramTools__

#include <iostream>
#include <stdexcept>

#include <JANA/JApplication.h>
#include <TH1I.h>
#include <TH2I.h>
#include <TProfile.h>
#include <TProfile2D.h>
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

vector <TDirectory *>& GetAllDirectories(void){
    static vector <TDirectory *> allDirectories;
    return allDirectories;
}

map<TString, TH1I*>& Get1DMap(void){
    static map<TString, TH1I*> TH1IMap;
    return TH1IMap;
}

map<TString, TH2I*>& Get2DMap(void){
    static map<TString, TH2I*> TH2IMap;
    return TH2IMap;
}

map<TString, TProfile*>& Get1DProfileMap(void){
    static map<TString, TProfile*> TProfile1DMap;
    return TProfile1DMap;
}

map<TString, TProfile2D*>& Get2DProfileMap(void){
    static map<TString, TProfile2D*> TProfile2DMap;
    return TProfile2DMap;
}

void Fill1DHistogram (const char * plugin, const char * directoryName, const char * name, const double value, const char * title , int nBins, double xmin, double xmax, bool print = false){
    japp->RootWriteLock();
    TH1I * histogram;
    TString fullName = TString(plugin) + "/" + TString(directoryName) + "/" + TString(name);
    try {
        histogram = Get1DMap().at(fullName);
    }
    catch (const std::out_of_range& oor) {
        if (print) std::cerr << ansi_green << plugin << " ===> Making New 1D Histogram " << name << ansi_normal << endl;
        TDirectory *homedir = gDirectory;
        TDirectory *temp;
        temp = gDirectory->mkdir(plugin);
        if(temp) GetAllDirectories().push_back(temp);
        gDirectory->cd(plugin);
        GetAllDirectories().push_back(gDirectory->mkdir(directoryName));
        gDirectory->cd(directoryName);
        Get1DMap()[fullName] = new TH1I( name, title, nBins, xmin, xmax);
        Get1DMap()[fullName]->Fill(value);
        homedir->cd();
        japp->RootUnLock();
        return;
    }
    histogram->Fill(value);
    japp->RootUnLock();
    return;
}

void Fill2DHistogram (const char * plugin, const char * directoryName, const char * name, const double valueX , const double valueY , const char * title , int nBinsX, double xmin, double xmax, int nBinsY, double ymin, double ymax, bool print = false){
    japp->RootWriteLock();
    TH2I * histogram;
    TString fullName = TString(plugin) + "/" + TString(directoryName) + "/" + TString(name);
    try {
        histogram = Get2DMap().at(fullName);
    }
    catch (const std::out_of_range& oor) {
        if (print) std::cerr << ansi_green << plugin << " ===> Making New 2D Histogram " << name << ansi_normal << endl;
        TDirectory *homedir = gDirectory;
        TDirectory *temp;
        temp = gDirectory->mkdir(plugin);
        if(temp) GetAllDirectories().push_back(temp);
        gDirectory->cd(plugin);
        GetAllDirectories().push_back(gDirectory->mkdir(directoryName));
        gDirectory->cd(directoryName);
        Get2DMap()[fullName] = new TH2I( name, title, nBinsX, xmin, xmax, nBinsY, ymin, ymax);
        Get2DMap()[fullName]->Fill(valueX, valueY);
        homedir->cd();
        japp->RootUnLock();
        return;
    }
    histogram->Fill(valueX, valueY);
    japp->RootUnLock();
    return;
}

void Fill1DProfile (const char * plugin, const char * directoryName, const char * name, const double valueX , const double valueY , const char * title , int nBinsX, double xmin, double xmax, bool print = false){
    japp->RootWriteLock();
    TProfile * profile;
    TString fullName = TString(plugin) + "/" + TString(directoryName) + "/" + TString(name);
    try {
        profile = Get1DProfileMap().at(fullName);
    }
    catch (const std::out_of_range& oor) {
        if (print) std::cerr << ansi_green << plugin << " ===> Making New 1D Profile " << name << ansi_normal << endl;
        TDirectory *homedir = gDirectory;
        TDirectory *temp;
        temp = gDirectory->mkdir(plugin);
        if(temp) GetAllDirectories().push_back(temp);
        gDirectory->cd(plugin);
        GetAllDirectories().push_back(gDirectory->mkdir(directoryName));
        gDirectory->cd(directoryName);
        Get1DProfileMap()[fullName] = new TProfile( name, title, nBinsX, xmin, xmax);
        Get1DProfileMap()[fullName]->Fill(valueX, valueY);
        homedir->cd();
        japp->RootUnLock();
        return;
    }
    profile->Fill(valueX, valueY);
    japp->RootUnLock();
    return;
}

void Fill2DProfile (const char * plugin, const char * directoryName, const char * name, const double valueX , const double valueY , const double valueZ, const char * title , int nBinsX, double xmin, double xmax, int nBinsY, double ymin, double ymax, bool print = false){
    japp->RootWriteLock();
    TProfile2D * profile;
    TString fullName = TString(plugin) + "/" + TString(directoryName) + "/" + TString(name);
    try {
        profile = Get2DProfileMap().at(fullName);
    }
    catch (const std::out_of_range& oor) {
        if (print) std::cerr << ansi_green << plugin << " ===> Making New 2D Profile " << name << ansi_normal << endl;
        TDirectory *homedir = gDirectory;
        TDirectory *temp;
        temp = gDirectory->mkdir(plugin);
        if(temp) GetAllDirectories().push_back(temp);
        gDirectory->cd(plugin);
        GetAllDirectories().push_back(gDirectory->mkdir(directoryName));
        gDirectory->cd(directoryName);
        Get2DProfileMap()[fullName] = new TProfile2D( name, title, nBinsX, xmin, xmax, nBinsY, ymin, ymax);
        Get2DProfileMap()[fullName]->Fill(valueX, valueY, valueZ);
        homedir->cd();
        japp->RootUnLock();
        return;
    }
    profile->Fill(valueX, valueY, valueZ);
    japp->RootUnLock();
    return;
}

TH1I * Get1DHistogram(const char * plugin, const char * directoryName, const char * name){
    TH1I * histogram;
    TString fullName = TString(plugin) + "/" + TString(directoryName) + "/" + TString(name);
    japp->RootWriteLock(); // Lock this thread down while we search for the histogram
    try {
        histogram = Get1DMap().at(fullName);
    }
    catch (const std::out_of_range& oor) {
        cout << ansi_red << ansi_blink << "Specified histogram " << fullName.Data() << " does not exist" << ansi_normal << endl;
        // Let's be nice and find the closest match
        TString matchName = "";
        int closestMatch = 1000000; //Just initialize high
        map<TString, TH1I*>::const_iterator iter;
        for(iter = Get1DMap().begin(); iter != Get1DMap().end(); iter++){
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
        histogram = Get2DMap().at(fullName);
    }
    catch (const std::out_of_range& oor) {
        cout << ansi_red << ansi_blink << "Specified histogram " << fullName.Data() << " does not exist" << ansi_normal << endl;
        // Let's be nice and find the closest match
        TString matchName = "";
        int closestMatch = 1000000; //Just initialize high
        map<TString, TH2I*>::const_iterator iter;
        for(iter = Get2DMap().begin(); iter != Get2DMap().end(); iter++){
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
    for (unsigned int i=0; i < GetAllDirectories().size(); i++){
        if (GetAllDirectories()[i] == 0) continue;
        japp->RootWriteLock();
        GetAllDirectories()[i]->GetList()->Sort();
        japp->RootUnLock();
    }
}

#endif
