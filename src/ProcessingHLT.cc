#include "UserCode/NtuplerAnalyzer/interface/ProcessingHLT.h"
#include "UserCode/NtuplerAnalyzer/interface/recordFuncs.h"
//#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

//using namespace std;

/*
 * should add this stuff:
 * utils::passTriggerPatterns
 *
 * from llvv repository
 */


/*
 * Find all trigger objects which match our pathName
 */


// select only our_hlt_trigger_objects
int Processing_selectHLTobjects(
	std::vector<pat::TriggerObjectStandAlone>& trig_objs,                  // input:  trigger objects in the event
	const edm::TriggerNames& trigNames,                               // input:  names for trigger bits in the event
	std::vector<pat::TriggerObjectStandAlone>& our_hlt_trigger_objects,    // output: trigger objects matching target HLT
	std::string& targetHLT // target HLT string-pattern (like "HLT_IsoMuon24_v" or "HLT_IsoMuon_v4")
	)

{
for (size_t i = 0; i < trig_objs.size(); i++)
	{
	pat::TriggerObjectStandAlone& obj = trig_objs[i];
	obj.unpackPathNames(trigNames);

	bool is_our_hlt = false;
	for (unsigned h = 0; h < obj.pathNames().size(); ++h)
		{
		// HLT_PFJet40_v
		is_our_hlt |= (obj.pathNames()[h].find(targetHLT) != std::string::npos);
		}

	if (is_our_hlt)
		our_hlt_trigger_objects.push_back(obj);
	}

return 0;
}



