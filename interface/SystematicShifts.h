#ifndef SYSTEMATICSHIFTS_H
#define SYSTEMATICSHIFTS_H

#include <map>
using namespace std;

/*
 * since cmssw puts everything into 1 shared lib
 * I have to make these declarations extern in the header
 * and move the definitions into a separate src file,
 * which will correspond to the part of shared lib for sysshifts
 */

typedef enum {SYS_NOMINAL,
	SYS_PU_UP, SYS_PU_DOWN,
	SYS_TOP_PT,
	SYS_JES_UP, SYS_JES_DOWN,
	SYS_JER_UP, SYS_JER_DOWN,
	SYS_BTAG_UP, SYS_BTAG_DOWN,
	SYS_MET_UNCLUSTERED_EN_UP, SYS_MET_UNCLUSTERED_EN_DOWN,
	SYS_TAUID_SF_UP, SYS_TAUID_SF_DOWN} systematic_shift;

static systematic_shift jetSystematics[] = {SYS_NOMINAL,
	SYS_JES_UP,
	SYS_JES_DOWN,
	SYS_JER_UP,
	SYS_JER_DOWN};

/* deprecated
 * TODO: all these systematic enums are deprecated now, remove the old functions which use these as input (in processing jets etc)
static systematic_shift weightSystematics[] = {SYS_NOMINAL,
	SYS_PU_UP, SYS_PU_DOWN,
	SYS_TOP_PT,
	SYS_BTAG_UP, SYS_BTAG_DOWN,
	SYS_TAUID_SF_UP, SYS_TAUID_SF_DOWN};

static systematic_shift btagSystematics[] = {SYS_NOMINAL,
	SYS_BTAG_UP,
	SYS_BTAG_DOWN};

static systematic_shift allSystematics[] = {SYS_NOMINAL,
	SYS_PU_UP, SYS_PU_DOWN,
	SYS_TOP_PT,
	SYS_JES_UP, SYS_JES_DOWN,
	SYS_JER_UP, SYS_JER_DOWN,
	SYS_BTAG_UP, SYS_BTAG_DOWN,
	SYS_MET_UNCLUSTERED_EN_UP, SYS_MET_UNCLUSTERED_EN_DOWN,
	SYS_TAUID_SF_UP, SYS_TAUID_SF_DOWN};

// C++11 feature:
static map<systematic_shift, const char*> systematic_shift_names = {{SYS_NOMINAL, "NOMINAL"},
	{SYS_PU_UP, "PU_UP"}, {SYS_PU_DOWN, "PU_DOWN"},
	{SYS_TOP_PT, "TOP_PT"},
	{SYS_JES_UP,   "JES_UP"}, {SYS_JES_DOWN, "JES_DOWN"},
	{SYS_JER_UP,   "JER_UP"}, {SYS_JER_DOWN, "JER_DOWN"},
	{SYS_BTAG_UP,   "BTAG_UP"}, {SYS_BTAG_DOWN, "BTAG_DOWN"},
	{SYS_MET_UNCLUSTERED_EN_UP,   "MET_UNCLUSTERED_EN_UP"}, {SYS_MET_UNCLUSTERED_EN_DOWN, "MET_UNCLUSTERED_EN_DOWN"},
	{SYS_TAUID_SF_UP, "SYS_TAUID_SF_UP"}, {SYS_TAUID_SF_DOWN, "SYS_TAUID_SF_DOWN"}
};

static map<systematic_shift, const string> btag_sys_points = {{SYS_NOMINAL, "central"},
	{SYS_BTAG_UP, "up"},
	{SYS_BTAG_DOWN, "down"}};
*/

/*
const char * systematic_shift_name[] = {
    "NOMINAL",
    "PU_UP",
    "PU_DOWN",
};
*/

#endif /* SYSTEMATICSHIFTS_H */

