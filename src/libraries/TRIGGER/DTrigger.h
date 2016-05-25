#ifndef _DTrigger_
#define _DTrigger_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>

class DTrigger : public jana::JObject
{
	public:
		JOBJECT_PUBLIC(DTrigger);

		//GETTERS
		uint32_t Get_L1TriggerBits(void) const{return dL1TriggerBits;}
		uint32_t Get_L1FrontPanelTriggerBits(void) const{return dL1FrontPanelTriggerBits;}
		bool Get_IsPhysicsEvent(void) const;

		//SETTERS
		void Set_L1TriggerBits(uint32_t locL1TriggerBits){dL1TriggerBits = locL1TriggerBits;}
		void Set_L1FrontPanelTriggerBits(uint32_t locL1FrontPanelTriggerBits){dL1FrontPanelTriggerBits = locL1FrontPanelTriggerBits;}

		// the second argument to AddString is printf style format
		void toStrings(vector<pair<string,string> >& items) const
		{
			AddString(items, "dL1TriggerBits", "%ld", dL1TriggerBits);
			AddString(items, "dL1FrontPanelTriggerBits", "%ld", dL1FrontPanelTriggerBits);
		}

	private:
		//NOTE: If is EPICS/SYNC/etc. event, both L1 values will be 0
		uint32_t dL1TriggerBits;
		uint32_t dL1FrontPanelTriggerBits;
};

inline bool DTrigger::Get_IsPhysicsEvent(void) const
{
	//Both L1 = 0: EPICS/SYNC/etc. //L1 = 8: PS
	return ((dL1FrontPanelTriggerBits == 0) && (dL1TriggerBits != 0) && (dL1TriggerBits != 8));
}

#endif // _DTrigger_
