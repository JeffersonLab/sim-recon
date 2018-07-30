/*
 * This script is designed to take a primary and secondary TTree produced
 * by the ReactionFilter plugin and merge the chi-square and ndf of the
 * secondary TTree into the primary.
 *
 * The TTrees should be similar such that some of the final state particles
 * could be swapped with each other, i.e. K+, K-, p and Pi+, Pi-, p
 *
 * Author: Alex Barnes
*/
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <chrono>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TObjArray.h"
#include "TBranch.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

using namespace std;
using namespace chrono;

vector<const char*> getBranches(TTree *tree)
{
   vector<const char*> locVector;

   TObjArray *arr = tree->GetListOfBranches();

   for (int i = 0; i < arr->GetEntries(); ++i)
   {
      TBranch *b = (TBranch*)arr->At(i);
      const char* b_name = b->GetName();

      // Get branches that vary between TTrees
      if ( strstr( b_name, "__ChargedIndex" ) != NULL ||
           strstr( b_name, "__NeutralIndex" ) != NULL )
      {
         locVector.push_back( b_name );
      }
      else
         continue; // branch didn't match with anything that should be kept
   }

   return locVector;
}

bool searchMap( map<string, Int_t> p_map, string substring )
{
   // This method searches a map to determine if it contains a final state particle matching the substring
   bool isMatch = false;
   for ( auto const &p_iter : p_map )
   {
      if ( p_iter.first.find(substring) != string::npos )
         isMatch = true;
   }
   
   return isMatch;
}

bool compareMaps( map<string, Int_t> p_map, map<string, Int_t> s_map )
{
   
   bool isMatch = false;
   int totalParticles = 0;
   int matchedParticles = 0;

   // check that maps are the same size
   if ( p_map.size() != s_map.size() )
      return isMatch;

   // loop over the primary TTree's map of (branch name, track ID)
   for ( auto const &p_iter : p_map )
   {
      totalParticles++;

      // loop over the secondary TTree's map of (branch name, track ID)
      for ( auto const &s_iter : s_map )
      {
         // compare the trackIDs
         if ( p_iter.second == s_iter.second )
            matchedParticles++;
      }
   }

   isMatch = ( totalParticles == matchedParticles ) ? true : false;

   return isMatch;

}

void mergeTrees(const char* primaryFile, const char* primaryTree, const char* secondaryFile, const char* secondaryTree)
{
   // For measuring performance
   high_resolution_clock::time_point t1 = high_resolution_clock::now();

   TFile f_secondary(secondaryFile);
   TTree *t_secondary = (TTree*)f_secondary.Get(secondaryTree);

   TFile f_orig(primaryFile);
   TTree *t_orig = (TTree*)f_orig.Get(primaryTree);

   TFile f_primary("newtree.root", "recreate");
   TTree *t_primary = (TTree*)t_orig->CloneTree();

   // Get branch information for each TTree
   vector<const char*> branches_primary = getBranches(t_primary);
   vector<const char*> branches_secondary = getBranches(t_secondary);

   // get reaction names and set up arrays for new branches
   UInt_t initArraySize = 100;
   Float_t new_chisq[initArraySize];
   UInt_t new_ndf[initArraySize];
   string p(primaryTree);
   string primaryReactionName = p.substr(0, p.find("_"));
   cout << "Primary reaction: " << primaryReactionName << endl;
   string s(secondaryTree);
   string secondaryReactionName = s.substr(0, s.find("_"));
   cout << "Secondary reaction: " << secondaryReactionName << endl;

   // primary tree
   TTreeReader primaryReader(t_primary);
   cout << "added primary tree" << endl;
   TTreeReaderArray<Float_t> chisq(primaryReader, "ChiSq_KinFit");
   TTreeReaderArray<UInt_t> ndf(primaryReader, "NDF_KinFit");
   TTreeReaderValue<UInt_t> run(primaryReader, "RunNumber");
   TTreeReaderValue<ULong64_t> event(primaryReader, "EventNumber");
   TTreeReaderArray<Int_t> ChargedHypoID(primaryReader, "ChargedHypo__TrackID");
   TTreeReaderArray<Int_t> ChargedHypo__PID(primaryReader, "ChargedHypo__PID");
   TTreeReaderArray<Int_t> NeutralHypoID(primaryReader, "NeutralHypo__NeutralID");
   TTreeReaderArray<Int_t> beam_ID(primaryReader, "ComboBeam__BeamIndex");
   // secondary tree
   TTreeReader secondaryReader(t_secondary);
   cout << "added secondary tree" << endl;
   TTreeReaderArray<Float_t> secondary_chisq(secondaryReader, "ChiSq_KinFit");
   TTreeReaderArray<UInt_t> secondary_ndf(secondaryReader, "NDF_KinFit");
   TTreeReaderValue<UInt_t> secondary_run(secondaryReader, "RunNumber");
   TTreeReaderValue<ULong64_t> secondary_event(secondaryReader, "EventNumber");
   TTreeReaderArray<Int_t> secondary_ChargedHypoID(secondaryReader, "ChargedHypo__TrackID");
   TTreeReaderArray<Int_t> secondary_ChargedHypo__PID(secondaryReader, "ChargedHypo__PID");
   TTreeReaderArray<Int_t> secondary_NeutralHypoID(primaryReader, "NeutralHypo__NeutralID");
   TTreeReaderArray<Int_t> secondary_beam_ID(secondaryReader, "ComboBeam__BeamIndex");

   // The only branches that will vary are the Int arrays
   // Map for the Int array branches
   map< const char*, TTreeReaderArray<Int_t> > primary_ReaderArray__Int;
   map< const char*, TTreeReaderArray<Int_t> > secondary_ReaderArray__Int;

   // Create TTreeReaderArrays for the Int_t branches 
   for ( auto const& p_iter : branches_primary )
   {
      TTreeReaderArray<Int_t> a(primaryReader, p_iter);
      auto P = make_pair( p_iter, a );
      primary_ReaderArray__Int.insert(P);
   }
   for ( auto const& s_iter : branches_secondary )
   {
      TTreeReaderArray<Int_t> a(secondaryReader, s_iter);
      auto P = make_pair( s_iter, a );
      secondary_ReaderArray__Int.insert(P);
   }
   cout << "Added to maps" << endl;
 
   // Get sizes of each map
   const int size_primary = (int)primary_ReaderArray__Int.size();
   const int size_secondary = (int)secondary_ReaderArray__Int.size();
   // Make arrays of the map keys
   const char* keys_primary[ size_primary ];
   const char* keys_secondary[ size_secondary ];
   // Get keys and fill arrays
   int count = 0;
   for ( auto const& p_iter : primary_ReaderArray__Int )
   {
      keys_primary[count] = p_iter.first;
      ++count;
   }
   count = 0;
   for ( auto const& s_iter : secondary_ReaderArray__Int )
   {
      keys_secondary[count] = s_iter.first;
      ++count;
   }

   // Create branches for secondary chisq and NDF
   string branchName("ChiSq_KinFit_");
   branchName += secondaryReactionName;
   TBranch *new_b_chisq = t_primary->Branch(branchName.c_str(), &new_chisq, (branchName + "[NumCombos]/F").c_str() );
   branchName = "NDF_KinFit_";
   branchName += secondaryReactionName;
   TBranch *new_b_ndf = t_primary->Branch(branchName.c_str(), &new_ndf, (branchName + "[NumCombos]/i").c_str() );

   // Make std::map< EventNumber, entry > for secondary TTree
   map< ULong64_t, Long64_t > map_event_secondary;
   while (secondaryReader.Next() )
   {
      ULong64_t thisEvent = *secondary_event;
      Long64_t thisEntry = secondaryReader.GetCurrentEntry();
      map_event_secondary.insert( pair< ULong64_t, Long64_t >(thisEvent, thisEntry ) );
   }

   // Find primary's EventNumber in secondary's TTree and get the desired information
   while (primaryReader.Next() )
   {
      Long64_t pentry = primaryReader.GetCurrentEntry();
      if ( pentry / 1000000 * 1000000 == pentry )
         cout << pentry / 1000000 + 1 << " million events analyzed" << endl;

      // TTreeReader.SetEntry() won't work unless the TTreeReader is restarted to the 0th entry
      secondaryReader.Restart();

      // Reset defaults for next event
      for (UInt_t i = 0; i < initArraySize; ++i)
      {
         new_chisq[i] = -1.0;
         new_ndf[i] = 1;
      }
       
      UInt_t p_numCombos = primary_ReaderArray__Int.begin()->second.GetSize();

      // Need to put in safe guard in case event doesn't exist
      if (map_event_secondary.find( *event ) == map_event_secondary.end() )
      {
         new_b_chisq->Fill();
         new_b_ndf->Fill();
         continue;
      }
      else
      {
         Long64_t entry = map_event_secondary.at( *event );
         secondaryReader.SetEntry( entry );
      }

      // Only get secondary number of combos if event exists in secondary TTree
      UInt_t s_numCombos = secondary_ReaderArray__Int.begin()->second.GetSize();

      // Loop over primary combos
      for (UInt_t i = 0; i < p_numCombos; ++i)
      {
         bool isMatched = false;

         // Get primary trackIDs
         map<string, Int_t> primary_trackIDs;
         map<string, Int_t> primary_neutralIDs;
         for ( int i_keys = 0; i_keys < size_primary; ++i_keys )
         {
            Int_t track;
            Int_t neutral;
            const char* b_name = keys_primary[ i_keys ];
            string s(b_name);
            if ( strstr( b_name, "__ChargedIndex" ) != NULL )
            {
               track = ChargedHypoID[ primary_ReaderArray__Int.at( b_name )[i] ];
               auto track_pair = make_pair( s, track );
               primary_trackIDs.insert( track_pair );
            }
            else if ( strstr( b_name, "__NeutralIndex" ) != NULL )
            {
               neutral = NeutralHypoID[ primary_ReaderArray__Int.at( b_name )[i] ];
               auto neutral_pair = make_pair( s, neutral );
               primary_neutralIDs.insert( neutral_pair );
            }
            else continue;
         }
         Int_t beam = beam_ID[i];

         // Loop over the secondary combos
         for (UInt_t j = 0; j < s_numCombos; ++j)
         {
            map<string, Int_t> secondary_trackIDs;
            map<string, Int_t> secondary_neutralIDs;
            for ( int j_keys = 0; j_keys < size_secondary; ++j_keys )
            {
               Int_t track;
               Int_t neutral;
               const char* b_name = keys_secondary[ j_keys ];
               string s(b_name);
               if ( strstr( b_name, "__ChargedIndex" ) != NULL )
               {
                  track = secondary_ChargedHypoID[ secondary_ReaderArray__Int.at( b_name )[j] ];
                  auto track_pair = make_pair( s, track );
                  secondary_trackIDs.insert(track_pair);
               }
               else if ( strstr( b_name, "__NeutralIndex" ) != NULL )
               {
                  neutral = secondary_NeutralHypoID[ secondary_ReaderArray__Int.at( b_name )[j] ];
                  auto neutral_pair = make_pair( s, neutral );
                  secondary_neutralIDs.insert( neutral_pair );
               }
               else continue;
            }
            Int_t secondary_beam = secondary_beam_ID[j];

            // check for the types of final state particles in the primary map
            bool hasPlus = searchMap( primary_trackIDs, "Plus" );
            bool hasMinus = searchMap( primary_trackIDs, "Minus" );
            bool hasProton = searchMap( primary_trackIDs, "Proton" ); // finds both protons and antiprotons
            bool hasPhoton = searchMap( primary_neutralIDs, "Photon" );

            // if the reaction does not contain a certain type of final state particle, set the match to true
            // otherwise, find the correct match
            bool chargedMatch = (hasPlus || hasMinus || hasProton) ? compareMaps( primary_trackIDs, secondary_trackIDs ) : true;
            bool neutralMatch = hasPhoton ? compareMaps( primary_neutralIDs, secondary_neutralIDs ) : true;

            // match final state particle IDs
            if ( !(chargedMatch && neutralMatch) )
               continue;
            // match beam ID
            if (beam != secondary_beam)
               continue;

            // grab secondary chisq, ndf information
            new_chisq[i] = secondary_chisq[j];
            new_ndf[i] = secondary_ndf[j];
            isMatched = true;
         }
      }

      new_b_chisq->Fill();
      new_b_ndf->Fill();

   }
   t_primary->Write("", TObject::kOverwrite); // save only the new version of the tree

   // For measuring performance
   high_resolution_clock::time_point t2 = high_resolution_clock::now();

   auto duration = duration_cast<microseconds>( t2 - t1 ).count();
   cout << "The process took: " << duration << " microseconds" << endl;
   cout << "The process took: " << duration*1E-6 << " seconds" << endl;
   cout << "The process took: " << duration*1E-6/60. << " minutes" << endl;

}
