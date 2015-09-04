#include "HistogramTools.h"

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

void Fill1DHistogram (const char * plugin, const char * directoryName, const char * name, const double value, const char * title , int nBins, double xmin, double xmax, bool print){
    japp->RootWriteLock();
    TH1I * histogram;
    TString fullName = TString(plugin) + "/" + TString(directoryName) + "/" + TString(name);
    try {
        histogram = Get1DMap().at(fullName);
    }
    catch (const std::out_of_range& oor) {
        if (print) std::cerr << ansi_green << plugin << " ===> Making New 1D Histogram " << name << ansi_normal << endl;
        // Ensure we are at the base directory
        string locOutputFileName = "hd_root.root";
        if(gPARMS->Exists("OUTPUT_FILENAME"))
            gPARMS->GetParameter("OUTPUT_FILENAME", locOutputFileName);
        TFile* locFile = (TFile*)gROOT->FindObject(locOutputFileName.c_str());
        if(locFile == NULL)
        {
            cout << "ERROR: OUTPUT HISTOGRAM FILE " << locOutputFileName << " NOT FOUND IN DParticleCombo_factory_PreKinFit::brun(). ABORTING." << endl;
            abort();
        }
        locFile->cd("");
        TDirectory *homedir = gDirectory;
        TDirectory *temp;
        temp = gDirectory->mkdir(plugin);
        if (temp) GetAllDirectories().push_back(temp);
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

void Fill2DHistogram (const char * plugin, const char * directoryName, const char * name, const double valueX , const double valueY , const char * title , int nBinsX, double xmin, double xmax, int nBinsY, double ymin, double ymax, bool print){
    japp->RootWriteLock();
    TH2I * histogram;
    TString fullName = TString(plugin) + "/" + TString(directoryName) + "/" + TString(name);
    try {
        histogram = Get2DMap().at(fullName);
    }
    catch (const std::out_of_range& oor) {
        if (print) std::cerr << ansi_green << plugin << " ===> Making New 2D Histogram " << name << ansi_normal << endl;
        // Ensure we are at the base directory
        string locOutputFileName = "hd_root.root";
        if(gPARMS->Exists("OUTPUT_FILENAME"))
            gPARMS->GetParameter("OUTPUT_FILENAME", locOutputFileName);
        TFile* locFile = (TFile*)gROOT->FindObject(locOutputFileName.c_str());
        if(locFile == NULL)
        {
            cout << "ERROR: OUTPUT HISTOGRAM FILE " << locOutputFileName << " NOT FOUND IN DParticleCombo_factory_PreKinFit::brun(). ABORTING." << endl;
            abort();
        }
        locFile->cd("");
        TDirectory *homedir = gDirectory;
        TDirectory *temp;
        temp = gDirectory->mkdir(plugin);
        if (temp) GetAllDirectories().push_back(temp);
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

