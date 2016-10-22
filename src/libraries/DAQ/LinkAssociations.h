// $Id:$
//
//    File: LinkAssociations.h
// Created: Tue May 17 08:15:00 EDT 2016
// Creator: davidl (on Linux gluon49.jlab.org 2.6.32-431.20.3.el6.x86_64)
//

#ifndef _LinkAssociations_
#define _LinkAssociations_

#include <stdint.h>
#include <vector>
using std::vector;

//==============================================================
// Sort Routines:
//
//  These are used to sort object lists based on various
// members. This is needed for some of the more complex
// matching routines below.

//----------------
// SortByROCID
//----------------
template<class T>
inline bool SortByROCID(const T* const &obj1, const T* const &obj2)
{
	return obj1->rocid < obj2->rocid;
}

//----------------
// SortByModule
//----------------
template<class T>
inline bool SortByModule(const T* const &obj1, const T* const &obj2)
{
	if(obj1->rocid < obj2->rocid) return true;
	if(obj1->rocid > obj2->rocid) return false;
	return obj1->slot < obj2->slot;
}

//----------------
// SortByChannel
//----------------
template<class T>
inline bool SortByChannel(const T* const &obj1, const T* const &obj2)
{
	if(obj1->rocid   < obj2->rocid   ) return true;
	if(obj1->rocid   > obj2->rocid   ) return false;
	if(obj1->slot    < obj2->slot    ) return true;
	if(obj1->slot    > obj2->slot    ) return false;
	return obj1->channel < obj2->channel;
}

//----------------
// SortByPulseNumber
//----------------
template<class T>
inline bool SortByPulseNumber(const T* const &obj1, const T* const &obj2)
{
	if(obj1->rocid   < obj2->rocid   ) return true;
	if(obj1->rocid   > obj2->rocid   ) return false;
	if(obj1->slot    < obj2->slot    ) return true;
	if(obj1->slot    > obj2->slot    ) return false;
	if(obj1->channel < obj2->channel ) return true;
	if(obj1->channel > obj2->channel ) return false;
	return obj1->pulse_number < obj2->pulse_number;
}


//==============================================================
// What follows are 2 sets of routines used to add objects to
// other objects' associated objects lists. The first set of routines
// are what should be called externally. These call an appropriate
// matching routine, passing it a lambda function for what to do
// once a match is found. This may be as simple as adding one to the
// other's associated objects list. There are also a few cases where
// we need to do some additional copying such as the pedestal from a
// PulsePedestal object into the pedestal member of a PulseIntegral
// object.


// LinkConfig
template<class T, class U>
inline void LinkConfig(vector<T*> &a, vector<U*> &b)
{ MatchConfigF(a, b, [](T *a, U *b){b->AddAssociatedObject(a);}); }

// LinkConfigSamplesCopy
template<class T, class U>
inline void LinkConfigSamplesCopy(vector<T*> &a, vector<U*> &b)
{ MatchConfigF(a, b, [](T *a, U *b){b->AddAssociatedObject(a); b->nsamples_integral = a->NSA_NSB; b->nsamples_pedestal = a->NPED;}); }

// LinkModule
template<class T, class U>
inline void LinkModule(vector<T*> &a, vector<U*> &b)
{ MatchModuleF(a, b, [](T *a, U *b){b->AddAssociatedObject(a);}); }

// LinkModuleBORSamplesCopy
template<class T, class U>
inline void LinkModuleBORSamplesCopy(vector<T*> &a, vector<U*> &b)
{ MatchModuleF(a, b, [](T *a, U *b){b->AddAssociatedObject(a); b->nsamples_integral = a->NSA+a->NSB; b->nsamples_pedestal = a->NPED;}); }

// LinkChannel
template<class T, class U>
inline void LinkChannel(vector<T*> &a, vector<U*> &b)
{ MatchChannelF(a, b, [](T *a, U *b){b->AddAssociatedObject(a); a->AddAssociatedObject(b);}); }

// LinkPulse
template<class T, class U>
inline void LinkPulse(vector<T*> &a, vector<U*> &b)
{ MatchPulseF(a, b, [](T *a, U *b){b->AddAssociatedObject(a);}); }

// LinkPulsePedCopy
template<class T, class U>
inline void LinkPulsePedCopy(vector<T*> &a, vector<U*> &b)
{ MatchPulseF(a, b, [](T *a, U *b){b->AddAssociatedObject(a); b->pedestal = a->pedestal;}); }

// PulsePedCopy
template<class T, class U>
inline void PulsePedCopy(vector<T*> &a, vector<U*> &b)
{ MatchPulseF(a, b, [](T *a, U *b){b->pedestal = a->pedestal;}); }


//----------------------------
// MatchConfig
//----------------------------
template<class T, class U, typename F>
void MatchConfig(vector<T*> &a, vector<U*> &b, F func)
{
	for(auto bptr : b){
		uint32_t slot_mask = 1 << bptr->slot;
		for(auto aptr : a){
			if(aptr->rocid!=bptr->rocid) continue;
			if(aptr->slot_mask & slot_mask){
				func(aptr, bptr);
			}
		}
	}
}

//----------------------------
// MatchModule
//----------------------------
template<class T, class U, typename F>
void MatchModule(vector<T*> &a, vector<U*> &b, F func)
{
	for(auto aptr : a){
		for(auto bptr : b){
			if(aptr->rocid != bptr->rocid) continue;
			if(aptr->slot  != bptr->slot ) continue;
			func(aptr, bptr);
		}
	}
}

//----------------------------
// MatchChannel
//----------------------------
template<class T, class U, typename F>
void MatchChannel(vector<T*> &a, vector<U*> &b, F func)
{
	for(auto aptr : a){
		for(auto bptr : b){
			if(aptr->rocid  != bptr->rocid  ) continue;
			if(aptr->slot   != bptr->slot   ) continue;
			if(aptr->channel!= bptr->channel) continue;
			func(aptr, bptr);
		}
	}
}

//----------------------------
// MatchPulse
//----------------------------
template<class T, class U, typename F>
void MatchPulse(vector<T*> &a, vector<U*> &b, F func)
{
	for(auto aptr : a){
		for(auto bptr : b){
			if(aptr->rocid        != bptr->rocid       ) continue;
			if(aptr->slot         != bptr->slot        ) continue;
			if(aptr->channel      != bptr->channel     ) continue;
			if(aptr->pulse_number != bptr->pulse_number) continue;
			func(aptr, bptr);
		}
	}
}

//==============================================================
// The following are alternatives to the matching routines
// above. These will step through both vectors simultaneously
// to find matching ranges as opposed to simply comparing each
// element of vector "a" with each element of vector "b". For this
// to work, both the a and b vectors must be sorted. The benefit
// here is a slight performance improvement (though it's a little
// hard to measure.) Providing the sorting has been done, these
// routines should each do the same thing as their counterparts
// that don't have the "F" suffix".


//----------------------------
// MatchConfigF
//----------------------------
template<class T, class U, typename F>
void MatchConfigF(vector<T*> &a, vector<U*> &b, F func)
{
	/// Template routine to loop over two vectors: The first points
	/// to DDAQConfig objects and the second to DDAQAddress objects.
	/// DDAQAddress objects are matched by the rocid and slotmask of
	/// the DDAQConfig object.
	///
	/// Note that this assumes the input vectors have been sorted by
	/// rocid in ascending order.

	// Bail early if nothing to link
	if(b.empty()) return;
	if(a.empty()) return;

	for(uint32_t i=0, j=0; i<b.size(); ){

		uint32_t rocid     = b[i]->rocid;

		// Find start and end of range in b
		uint32_t istart = i;
		uint32_t iend   = i+1; // index of first element outside of ROI
		for(; iend<b.size(); iend++){
			if(b[iend]->rocid != rocid) break;
		}
		i = iend; // setup for next iteration
		
		// Find start of range in a
		for(; j<a.size(); j++){
			if( a[j]->rocid >= rocid ) break;
		}
		if(j>=a.size()) break; // exhausted all a's. we're done
		if( a[j]->rocid > rocid ) continue; // couldn't find rocid in a
		
		// Find end of range in a
		uint32_t jend = j+1;
		for(; jend<a.size(); jend++){
			if(a[jend]->rocid != rocid) break;
		}

		// Loop over all combos of both ranges and make associations
		uint32_t jstart = j;
		for(uint32_t ii=istart; ii<iend; ii++){
			uint32_t slot_mask = 1 << b[ii]->slot;

			for(uint32_t jj=jstart; jj<jend; jj++){
				if(a[jj]->slot_mask & slot_mask){
					func(a[jj], b[ii]);
				}
			}
		}
		j = jend;

		if( i>=b.size() ) break;
		if( j>=a.size() ) break;
	}
}

//----------------------------
// MatchModuleF
//----------------------------
template<class T, class U, typename F>
void MatchModuleF(vector<T*> &a, vector<U*> &b, F func)
{
	/// Template routine to loop over two vectors of pointers to
	/// objects derived from DDAQAddress. This will match any hits
	/// coming from the same DAQ module (channel number is not checked)
	///
	/// Note that this assumes the input vectors have been sorted by
	/// rocid, then slot in ascending order.

	// Bail early if nothing to link
	if(b.empty()) return;
	if(a.empty()) return;

	for(uint32_t i=0, j=0; i<b.size(); ){

		uint32_t rocid = b[i]->rocid;
		uint32_t slot  = b[i]->slot;

		// Find start and end of range in b
		uint32_t istart = i;
		uint32_t iend   = i+1; // index of first element outside of ROI
		for(; iend<b.size(); iend++){
			if(b[iend]->rocid != rocid) break;
			if(b[iend]->slot  != slot ) break;
		}
		i = iend; // setup for next iteration
		
		// Find start of range in a
		for(; j<a.size(); j++){
			if( a[j]->rocid > rocid ) break;
			if( a[j]->rocid == rocid ){
				if( a[j]->slot >= slot ) break;
			}
		}
		if(j>=a.size()) break; // exhausted all a's. we're done

		if( a[j]->rocid > rocid ) continue; // couldn't find rocid in a
		if( a[j]->slot  > slot  ) continue; // couldn't find slot in a
		
		// Find end of range in a
		uint32_t jend = j+1;
		for(; jend<a.size(); jend++){
			if(a[jend]->rocid != rocid) break;
			if(a[jend]->slot  != slot ) break;
		}

		// Loop over all combos of both ranges and make associations
		uint32_t jstart = j;
		for(uint32_t ii=istart; ii<iend; ii++){
			for(uint32_t jj=jstart; jj<jend; jj++){
				func(a[jj], b[ii]);
			}
		}
		j = jend;

		if( i>=b.size() ) break;
		if( j>=a.size() ) break;
	}
}

//----------------------------
// MatchChannelF
//----------------------------
template<class T, class U, typename F>
void MatchChannelF(vector<T*> &a, vector<U*> &b, F func)
{
	/// Template routine to loop over two vectors of pointers to
	/// objects derived from DDAQAddress. This will match any hits
	/// coming from the same DAQ module and channel.
	///
	/// Note that this assumes the input vectors have been sorted by
	/// rocid, then slot in ascending order.

	// Bail early if nothing to link
	if(b.empty()) return;
	if(a.empty()) return;

	for(uint32_t i=0, j=0; i<b.size(); ){

		uint32_t rocid        = b[i]->rocid;
		uint32_t slot         = b[i]->slot;
		uint32_t channel      = b[i]->channel;

		// Find start and end of range in b
		uint32_t istart = i;
		uint32_t iend   = i+1; // index of first element outside of ROI
		for(; iend<b.size(); iend++){
			if(b[iend]->rocid         != rocid        ) break;
			if(b[iend]->slot          != slot         ) break;
			if(b[iend]->channel       != channel      ) break;
		}
		i = iend; // setup for next iteration
		
		// Find start of range in a
		for(; j<a.size(); j++){
			if( a[j]->rocid > rocid ) break;
			if( a[j]->rocid == rocid ){
				if( a[j]->slot > slot ) break;
				if( a[j]->slot == slot ){
					if( a[j]->channel >= channel ) break;
				}
			}
		}
		if(j>=a.size()) break; // exhausted all a's. we're done

		if( a[j]->rocid        > rocid         ) continue; // couldn't find rocid in a
		if( a[j]->slot         > slot          ) continue; // couldn't find slot in a
		if( a[j]->channel      > channel       ) continue; // couldn't find channel in a
		
		// Find end of range in a
		uint32_t jend = j+1;
		for(; jend<a.size(); jend++){
			if(a[jend]->rocid        != rocid        ) break;
			if(a[jend]->slot         != slot         ) break;
			if(a[jend]->channel      != channel      ) break;
		}

		// Loop over all combos of both ranges and make associations
		uint32_t jstart = j;
		for(uint32_t ii=istart; ii<iend; ii++){
			for(uint32_t jj=jstart; jj<jend; jj++){
				func(a[jj], b[ii]);
			}
		}
		j = jend;

		if( i>=b.size() ) break;
		if( j>=a.size() ) break;
	}
}

//----------------------------
// MatchPulseF
//----------------------------
template<class T, class U, typename F>
void MatchPulseF(vector<T*> &a, vector<U*> &b, F func)
{
	/// Template routine to loop over two vectors of pointers to
	/// objects derived from DDAQAddress. This will match any hits
	/// coming from the same DAQ module, channel with the same pulse
	/// number.
	///
	/// Note that this assumes the input vectors have been sorted by
	/// rocid, then slot in ascending order.

	// Bail early if nothing to link
	if(b.empty()) return;
	if(a.empty()) return;

	for(uint32_t i=0, j=0; i<b.size(); ){

		uint32_t rocid        = b[i]->rocid;
		uint32_t slot         = b[i]->slot;
		uint32_t channel      = b[i]->channel;
		uint32_t pulse_number = b[i]->pulse_number;

		// Find start and end of range in b
		uint32_t istart = i;
		uint32_t iend   = i+1; // index of first element outside of ROI
		for(; iend<b.size(); iend++){
			if(b[iend]->rocid         != rocid        ) break;
			if(b[iend]->slot          != slot         ) break;
			if(b[iend]->channel       != channel      ) break;
			if(b[iend]->pulse_number  != pulse_number ) break;
		}
		i = iend; // setup for next iteration
		
		// Find start of range in a
		for(; j<a.size(); j++){
			if( a[j]->rocid > rocid ) break;
			if( a[j]->rocid == rocid ){
				if( a[j]->slot > slot ) break;
				if( a[j]->slot == slot ){
					if( a[j]->channel > channel ) break;
					if( a[j]->channel == channel ){
						if( a[j]->pulse_number >= pulse_number ) break;
					}
				}
			}
		}
		if(j>=a.size()) break; // exhausted all a's. we're done

		if( a[j]->rocid        > rocid         ) continue; // couldn't find rocid in a
		if( a[j]->slot         > slot          ) continue; // couldn't find slot in a
		if( a[j]->channel      > channel       ) continue; // couldn't find channel in a
		if( a[j]->pulse_number > pulse_number  ) continue; // couldn't find pulse_number in a
		
		// Find end of range in a
		uint32_t jend = j+1;
		for(; jend<a.size(); jend++){
			if(a[jend]->rocid        != rocid        ) break;
			if(a[jend]->slot         != slot         ) break;
			if(a[jend]->channel      != channel      ) break;
			if(a[jend]->pulse_number != pulse_number ) break;
		}

		// Loop over all combos of both ranges and make associations
		uint32_t jstart = j;
		for(uint32_t ii=istart; ii<iend; ii++){
			for(uint32_t jj=jstart; jj<jend; jj++){
				func(a[jj], b[ii]);
			}
		}
		j = jend;

		if( i>=b.size() ) break;
		if( j>=a.size() ) break;
	}
}

#endif // _LinkAssociations_

