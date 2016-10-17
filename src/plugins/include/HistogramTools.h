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

map<TString, pair<TH1I*, pthread_rwlock_t*> >& Get1DMap(void){
   static map<TString, pair<TH1I*, pthread_rwlock_t*> > TH1IMap;
   return TH1IMap;
}

map<TString, pair<TH1D*, pthread_rwlock_t*> >& Get1DWeightedMap(void){
   static map<TString, pair<TH1D*, pthread_rwlock_t*> > TH1DMap;
   return TH1DMap;
}

map<TString, pair<TH2I*, pthread_rwlock_t*> >& Get2DMap(void){
   static map<TString, pair<TH2I*, pthread_rwlock_t*> > TH2IMap;
   return TH2IMap;
}

map<TString, pair<TH2D*, pthread_rwlock_t*> >& Get2DWeightedMap(void){
   static map<TString, pair<TH2D*, pthread_rwlock_t*> > TH2DMap;
   return TH2DMap;
}

map<TString, pair<TProfile*, pthread_rwlock_t*> >& Get1DProfileMap(void){
   static map<TString, pair<TProfile*, pthread_rwlock_t*> > TProfile1DMap;
   return TProfile1DMap;
}

map<TString, pair<TProfile2D*, pthread_rwlock_t*> >& Get2DProfileMap(void){
   static map<TString, pair<TProfile2D*, pthread_rwlock_t*> > TProfile2DMap;
   return TProfile2DMap;
}

pthread_rwlock_t* InitializeMapLock(){
   pthread_rwlock_t *thisLock = new pthread_rwlock_t();
   pthread_rwlock_init(thisLock, NULL);
   return thisLock;
}

void Fill1DHistogram (const char * plugin, const char * directoryName, const char * name, const double value, const char * title , int nBins, double xmin, double xmax, bool print = false){
   static pthread_rwlock_t *mapLock = InitializeMapLock(); 
   TH1I * histogram;
   pair<TH1I*, pthread_rwlock_t*> histogramPair;

   char fullNameChar[500];
   sprintf(fullNameChar, "%s/%s/%s", plugin, directoryName, name); 
   TString fullName = TString(fullNameChar);

   try {
      pthread_rwlock_rdlock(mapLock); // Grab the read lock
      histogramPair = Get1DMap().at(fullName); // If the exception is caught, it bails immediately and will not release the lock
      pthread_rwlock_unlock(mapLock);
   }
   catch (const std::out_of_range& oor) {
      // Drop the read lock and grab the write lock
      pthread_rwlock_unlock(mapLock);
      pthread_rwlock_wrlock(mapLock);
      // At this point, more than one thread might have made it through the try and ended up in the catch. 
      // do a single threaded (write locked) "try" again to be sure we aren't duplicating the histogram
      try {
         histogramPair = Get1DMap().at(fullName); 
      }
      catch(const std::out_of_range& oor) {
         // This now must be the first thread to get here, so we should make the histogram, fill and move on
         if (print) std::cerr << ansi_green << plugin << " ===> Making New 1D Histogram " << fullName << ansi_normal << endl;
         // Initialize the histogram lock
         pthread_rwlock_t *histogramLock = new pthread_rwlock_t();
         pthread_rwlock_init(histogramLock, NULL);
         // Get the ROOT lock and create the histogram
         // WARNING: Locking inside a lock is bad practice, but sometimes not easy to avoid.
         // there would be a problem if there was another function that tried to grab the map lock
         // inside a root lock. In this code, this will not happen.
         japp->RootWriteLock();
         TDirectory *homedir = gDirectory;
         TDirectory *temp;
         temp = gDirectory->mkdir(plugin);
         if(temp) GetAllDirectories().push_back(temp);
         gDirectory->cd(plugin);
         GetAllDirectories().push_back(gDirectory->mkdir(directoryName));
         gDirectory->cd(directoryName);
         histogram = new TH1I( name, title, nBins, xmin, xmax);
         histogram->Fill(value);
         homedir->cd();
         japp->RootUnLock();

         Get1DMap()[fullName] = make_pair(histogram, histogramLock);
         pthread_rwlock_unlock(mapLock);
         return;
      }
      // If nothing is caught, the histogram must have been created by another thread 
      // while we were waiting to grab the write lock. Drop the lock and carry on...
      pthread_rwlock_unlock(mapLock);
   }

   histogram = histogramPair.first;
   pthread_rwlock_t *histogramLockPtr = histogramPair.second;
   pthread_rwlock_wrlock(histogramLockPtr);
   histogram->Fill(value);
   pthread_rwlock_unlock(histogramLockPtr);

   return;
}

void Fill1DWeightedHistogram (const char * plugin, const char * directoryName, const char * name, const double value, const double weight, const char * title , int nBins, double xmin, double xmax, bool print = false){
    static pthread_rwlock_t *mapLock = InitializeMapLock(); 
    TH1D * histogram;
    pair<TH1D*, pthread_rwlock_t*> histogramPair;

    char fullNameChar[500];
    sprintf(fullNameChar, "%s/%s/%s", plugin, directoryName, name); 
    TString fullName = TString(fullNameChar);

    try {
        pthread_rwlock_rdlock(mapLock); // Grab the read lock
        histogramPair = Get1DWeightedMap().at(fullName); // If the exception is caught, it bails immediately and will not release the lock
        pthread_rwlock_unlock(mapLock);
    }
    catch (const std::out_of_range& oor) {
        // Drop the read lock and grab the write lock
        pthread_rwlock_unlock(mapLock);
        pthread_rwlock_wrlock(mapLock);
        // At this point, more than one thread might have made it through the try and ended up in the catch. 
        // do a single threaded (write locked) "try" again to be sure we aren't duplicating the histogram
        try {
            histogramPair = Get1DWeightedMap().at(fullName); 
        }
        catch(const std::out_of_range& oor) {
            // This now must be the first thread to get here, so we should make the histogram, fill and move on
            if (print) std::cerr << ansi_green << plugin << " ===> Making New 1D Weighted Histogram " << fullName << ansi_normal << endl;
            // Initialize the histogram lock
            pthread_rwlock_t *histogramLock = new pthread_rwlock_t();
            pthread_rwlock_init(histogramLock, NULL);
            // Get the ROOT lock and create the histogram
            // WARNING: Locking inside a lock is bad practice, but sometimes not easy to avoid.
            // there would be a problem if there was another function that tried to grab the map lock
            // inside a root lock. In this code, this will not happen.
            japp->RootWriteLock();
            TDirectory *homedir = gDirectory;
            TDirectory *temp;
            temp = gDirectory->mkdir(plugin);
            if(temp) GetAllDirectories().push_back(temp);
            gDirectory->cd(plugin);
            GetAllDirectories().push_back(gDirectory->mkdir(directoryName));
            gDirectory->cd(directoryName);
            histogram = new TH1D( name, title, nBins, xmin, xmax);
            histogram->Fill(value,weight);
            homedir->cd();
            japp->RootUnLock();

            Get1DWeightedMap()[fullName] = make_pair(histogram, histogramLock);
            pthread_rwlock_unlock(mapLock);
            return;
        }
        // If nothing is caught, the histogram must have been created by another thread 
        // while we were waiting to grab the write lock. Drop the lock and carry on...
        pthread_rwlock_unlock(mapLock);
    }

    histogram = histogramPair.first;
    pthread_rwlock_t *histogramLockPtr = histogramPair.second;
    pthread_rwlock_wrlock(histogramLockPtr);
    histogram->Fill(value,weight);
    pthread_rwlock_unlock(histogramLockPtr);

    return;
}

void Fill2DHistogram (const char * plugin, const char * directoryName, const char * name, const double valueX , const double valueY , const char * title , int nBinsX, double xmin, double xmax, int nBinsY, double ymin, double ymax, bool print = false){

   static pthread_rwlock_t *mapLock = InitializeMapLock();
   TH2I * histogram;
   pair<TH2I*, pthread_rwlock_t*> histogramPair;

   char fullNameChar[500];
   sprintf(fullNameChar, "%s/%s/%s", plugin, directoryName, name);
   TString fullName = TString(fullNameChar);

   try {
      pthread_rwlock_rdlock(mapLock); // Grab the read lock
      histogramPair = Get2DMap().at(fullName);
      pthread_rwlock_unlock(mapLock); 
   }
   catch (const std::out_of_range& oor) {
      // Drop the read lock and grab the write lock
      pthread_rwlock_unlock(mapLock);
      pthread_rwlock_wrlock(mapLock);
      // At this point, more than one thread might have made it through the try and ended up in the catch.
      // do a single threaded (write locked) "try" again to be sure we aren't duplicating the histogram
      try{
         histogramPair = Get2DMap().at(fullName);
      }
      catch  (const std::out_of_range& oor) {
         if (print) std::cerr << ansi_green << plugin << " ===> Making New 2D Histogram " << name << ansi_normal << endl;
         // Initialize the histogram lock
         pthread_rwlock_t *histogramLock = new pthread_rwlock_t();
         pthread_rwlock_init(histogramLock, NULL);

         // Get the ROOT lock and create the histogram
         // WARNING: Locking inside a lock is bad practice, but sometimes not easy to avoid.
         // there would be a problem if there was another function that tried to grab the map lock
         // inside a root lock. In this code, this will not happen.
         japp->RootWriteLock();
         TDirectory *homedir = gDirectory;
         TDirectory *temp;
         temp = gDirectory->mkdir(plugin);
         if(temp) GetAllDirectories().push_back(temp);
         gDirectory->cd(plugin);
         GetAllDirectories().push_back(gDirectory->mkdir(directoryName));
         gDirectory->cd(directoryName);
         histogram = new TH2I( name, title, nBinsX, xmin, xmax, nBinsY, ymin, ymax);
         histogram->Fill(valueX, valueY);
         homedir->cd();
         japp->RootUnLock();

         Get2DMap()[fullName] = make_pair(histogram,histogramLock);
         pthread_rwlock_unlock(mapLock);
         return;
      }
      // If nothing is caught, the histogram must have been created by another thread
      // while we were waiting to grab the write lock. Drop the lock and carry on...
      pthread_rwlock_unlock(mapLock);
   }

   histogram = histogramPair.first;
   pthread_rwlock_t *histogramLockPtr = histogramPair.second;
   pthread_rwlock_wrlock(histogramLockPtr);
   histogram->Fill(valueX, valueY);
   pthread_rwlock_unlock(histogramLockPtr);

   return;
}

void Fill2DWeightedHistogram (const char * plugin, const char * directoryName, const char * name, const double valueX , const double valueY , const double weight , const char * title , int nBinsX, double xmin, double xmax, int nBinsY, double ymin, double ymax, bool print = false){

    static pthread_rwlock_t *mapLock = InitializeMapLock();
    TH2D * histogram;
    pair<TH2D*, pthread_rwlock_t*> histogramPair;

    char fullNameChar[500];
    sprintf(fullNameChar, "%s/%s/%s", plugin, directoryName, name);
    TString fullName = TString(fullNameChar);

    try {
        pthread_rwlock_rdlock(mapLock); // Grab the read lock
        histogramPair = Get2DWeightedMap().at(fullName);
        pthread_rwlock_unlock(mapLock); 
    }
    catch (const std::out_of_range& oor) {
        // Drop the read lock and grab the write lock
        pthread_rwlock_unlock(mapLock);
        pthread_rwlock_wrlock(mapLock);
        // At this point, more than one thread might have made it through the try and ended up in the catch.
        // do a single threaded (write locked) "try" again to be sure we aren't duplicating the histogram
        try{
            histogramPair = Get2DWeightedMap().at(fullName);
        }
        catch  (const std::out_of_range& oor) {
            if (print) std::cerr << ansi_green << plugin << " ===> Making New 2D Histogram " << name << ansi_normal << endl;
            // Initialize the histogram lock
            pthread_rwlock_t *histogramLock = new pthread_rwlock_t();
            pthread_rwlock_init(histogramLock, NULL);

            // Get the ROOT lock and create the histogram
            // WARNING: Locking inside a lock is bad practice, but sometimes not easy to avoid.
            // there would be a problem if there was another function that tried to grab the map lock
            // inside a root lock. In this code, this will not happen.
            japp->RootWriteLock();
            TDirectory *homedir = gDirectory;
            TDirectory *temp;
            temp = gDirectory->mkdir(plugin);
            if(temp) GetAllDirectories().push_back(temp);
            gDirectory->cd(plugin);
            GetAllDirectories().push_back(gDirectory->mkdir(directoryName));
            gDirectory->cd(directoryName);
            histogram = new TH2D( name, title, nBinsX, xmin, xmax, nBinsY, ymin, ymax);
            histogram->Fill(valueX, valueY, weight);
            homedir->cd();
            japp->RootUnLock();

            Get2DWeightedMap()[fullName] = make_pair(histogram,histogramLock);
            pthread_rwlock_unlock(mapLock);
            return;
        }
        // If nothing is caught, the histogram must have been created by another thread
        // while we were waiting to grab the write lock. Drop the lock and carry on...
        pthread_rwlock_unlock(mapLock);
    }

    histogram = histogramPair.first;
    pthread_rwlock_t *histogramLockPtr = histogramPair.second;
    pthread_rwlock_wrlock(histogramLockPtr);
    histogram->Fill(valueX, valueY, weight);
    pthread_rwlock_unlock(histogramLockPtr);

    return;
}

void Fill1DProfile (const char * plugin, const char * directoryName, const char * name, const double valueX , const double valueY , const char * title , int nBinsX, double xmin, double xmax, bool print = false){

   static pthread_rwlock_t *mapLock = InitializeMapLock();
   TProfile * profile;
   pair<TProfile*, pthread_rwlock_t*> profilePair;

   char fullNameChar[500];
   sprintf(fullNameChar, "%s/%s/%s", plugin, directoryName, name);
   TString fullName = TString(fullNameChar);
   try {
      pthread_rwlock_rdlock(mapLock); // Grab the read lock
      profilePair = Get1DProfileMap().at(fullName);
      pthread_rwlock_unlock(mapLock); 
   }
   catch (const std::out_of_range& oor) {
      // Drop the read lock and grab the write lock
      pthread_rwlock_unlock(mapLock);
      pthread_rwlock_wrlock(mapLock);
      // At this point, more than one thread might have made it through the try and ended up in the catch.
      // do a single threaded (write locked) "try" again to be sure we aren't duplicating the histogram
      try{
         profilePair = Get1DProfileMap().at(fullName);
      }
      catch (const std::out_of_range& oor) {
         if (print) std::cerr << ansi_green << plugin << " ===> Making New 1D Profile " << name << ansi_normal << endl;
         // Initialize the profile lock
         pthread_rwlock_t *profileLock = new pthread_rwlock_t();
         pthread_rwlock_init(profileLock, NULL);

         // Get the ROOT lock and create the histogram
         // WARNING: Locking inside a lock is bad practice, but sometimes not easy to avoid.
         // there would be a problem if there was another function that tried to grab the map lock
         // inside a root lock. In this code, this will not happen.
         japp->RootWriteLock();// Get the ROOT lock and create the histogram
         TDirectory *homedir = gDirectory;
         TDirectory *temp;
         temp = gDirectory->mkdir(plugin);
         if(temp) GetAllDirectories().push_back(temp);
         gDirectory->cd(plugin);
         GetAllDirectories().push_back(gDirectory->mkdir(directoryName));
         gDirectory->cd(directoryName);
         profile = new TProfile( name, title, nBinsX, xmin, xmax);
         profile->Fill(valueX, valueY);
         homedir->cd();
         japp->RootUnLock();

         Get1DProfileMap()[fullName] = make_pair(profile,profileLock);
         pthread_rwlock_unlock(mapLock);
         return;
      }
      // If nothing is caught, the histogram must have been created by another thread
      // while we were waiting to grab the write lock. Drop the lock and carry on...
      pthread_rwlock_unlock(mapLock);
   }

   profile = profilePair.first;
   pthread_rwlock_t *profileLockPtr = profilePair.second;
   pthread_rwlock_wrlock(profileLockPtr);
   profile->Fill(valueX, valueY);
   pthread_rwlock_unlock(profileLockPtr);

   return;
}

void Fill2DProfile (const char * plugin, const char * directoryName, const char * name, const double valueX , const double valueY , const double valueZ, const char * title , int nBinsX, double xmin, double xmax, int nBinsY, double ymin, double ymax, bool print = false){
   static pthread_rwlock_t *mapLock = InitializeMapLock();
   TProfile2D * profile;
   pair<TProfile2D*, pthread_rwlock_t*> profilePair;

   char fullNameChar[500];
   sprintf(fullNameChar, "%s/%s/%s", plugin, directoryName, name);
   TString fullName = TString(fullNameChar);
   try {
      pthread_rwlock_rdlock(mapLock); // Grab the read lock
      profilePair = Get2DProfileMap().at(fullName);
      pthread_rwlock_unlock(mapLock); 
   }
   catch (const std::out_of_range& oor) {
      // Drop the read lock and grab the write lock
      pthread_rwlock_unlock(mapLock);
      pthread_rwlock_wrlock(mapLock);
      // At this point, more than one thread might have made it through the try and ended up in the catch.
      // do a single threaded (write locked) "try" again to be sure we aren't duplicating the histogram
      try{
         profilePair = Get2DProfileMap().at(fullName);
      }
      catch (const std::out_of_range& oor) {
         if (print) std::cerr << ansi_green << plugin << " ===> Making New 2D Profile " << name << ansi_normal << endl;
         // Initialize the profile lock
         pthread_rwlock_t *profileLock = new pthread_rwlock_t();
         pthread_rwlock_init(profileLock, NULL);

         // Get the ROOT lock and create the histogram
         // WARNING: Locking inside a lock is bad practice, but sometimes not easy to avoid. 
         // there would be a problem if there was another function that tried to grab the map lock 
         // inside a root lock. In this code, this will not happen. 
         japp->RootWriteLock();
         TDirectory *homedir = gDirectory;
         TDirectory *temp;
         temp = gDirectory->mkdir(plugin);
         if(temp) GetAllDirectories().push_back(temp);
         gDirectory->cd(plugin);
         GetAllDirectories().push_back(gDirectory->mkdir(directoryName));
         gDirectory->cd(directoryName);
         profile = new TProfile2D( name, title, nBinsX, xmin, xmax, nBinsY, ymin, ymax);
         profile->Fill(valueX, valueY, valueZ);
         homedir->cd();
         japp->RootUnLock();

         Get2DProfileMap()[fullName] = make_pair(profile, profileLock);
         pthread_rwlock_unlock(mapLock);
         return;
      }
      // If nothing is caught, the histogram must have been created by another thread
      // while we were waiting to grab the write lock. Drop the lock and carry on...
      pthread_rwlock_unlock(mapLock);
   }

   profile = profilePair.first;
   pthread_rwlock_t *profileLockPtr = profilePair.second;
   pthread_rwlock_wrlock(profileLockPtr);
   profile->Fill(valueX, valueY, valueZ);
   pthread_rwlock_unlock(profileLockPtr);

   return;
}

void SortDirectories(){
   japp->RootWriteLock();
   for (unsigned int i=0; i < GetAllDirectories().size(); i++){
      if (GetAllDirectories()[i] == 0) continue;
      GetAllDirectories()[i]->GetList()->Sort();
   }
   japp->RootUnLock();
}

#endif
