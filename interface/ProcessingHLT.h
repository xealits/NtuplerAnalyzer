#ifndef PROCESSINGHLT_H
#define PROCESSINGHLT_H

// for trigger matching:
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
// compiled with this:
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "UserCode/NtuplerAnalyzer/interface/common_definitions.h"

int Processing_selectHLTobjects(
	std::vector<pat::TriggerObjectStandAlone>& trig_objs,                  // input:  trigger objects in the event
	const edm::TriggerNames& trigNames,                               // input:  names for trigger bits in the event
	std::vector<pat::TriggerObjectStandAlone>& our_hlt_trigger_objects,    // output: trigger objects matching target HLT
	std::string& targetHLT // target HLT string-pattern (like "HLT_IsoMuon24_v" or "HLT_IsoMuon_v4")
	);

/*
 * ROOT error
 * not found
 * Processing_selectHLTobjects(std::vector<pat::TriggerObjectStandAlone, std::allocator<pat::TriggerObjectStandAlone> >&, edm::TriggerNames const&, std::vector<pat::TriggerObjectStandAlone, std::allocator<pat::TriggerObjectStandAlone> >&, std::string)
 *
 * can't solve it for 4 hours
 * some CMSSW feature probably
 *
 * don't see any difference between ProcessingHLT and ProcessingElectrons files
 *
 * now:
 * ProcessingHLT.cc:21:5: error: redefinition of 'int Processing_selectHLTobjects...
 * ProcessingHLT.h:28:5: note: 'int Processing_selectHLTobjects(std::vector<pat::TriggerObjectStandAlone>&, const edm::TriggerNames&, std::vector<pat::TriggerObjectStandAlone>&, std::string&)' previously defined here
 *
 * CMSSW, scram, awesomeness
 *
 * there is no magic in programming, just bad error reports
 * 4 hours searching for inconsistent build process of scram
 * -- the error was in declared std::string and defined std::string& input
 * C++ same function names are overloaded and no warning is reported..
 * no scram documentation leaves you hanging on whether your knowledge is complete on how scram and CMSSW buiuld process works
 *
 * good thing: apparently, I got (hopefully) all details of scram's operation
 * yay!
 */

#endif /* PROCESSINGHLT_H */

