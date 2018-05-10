#include <RoccoR.h>
//#include "RoccoR.cc"

RoccoR  rc("rcdata.2016.v3");

namespace roccor_wrapper
        {

	extern "C" double wrapper_kScaleDT(int Q, double pt, double eta, double phi, int s, int m) {
	    return rc.kScaleDT(Q, pt, eta, phi, s, m);
	}

	extern "C" int wrapper_kScaleDT_double(double *res, int Q, double pt, double eta, double phi, int s, int m) {
	    *res = rc.kScaleDT(Q, pt, eta, phi, s, m);
	    return 55;
	}

	//for MC, if matched gen-level muon (genPt) is available, use this function
	extern "C" double wrapper_kScaleFromGenMC(int Q, double pt, double eta, double phi, int nl, double genPt, double u1, int s, int m) {
	    //return rc.kScaleFromGenMC(Q, pt, eta, phi, nl, genPt, u1, s=0, m=0);
	    //return rc.kScaleFromGenMC(Q, pt, eta, phi, nl, genPt, u1, s=0, m=0);
	    return rc.kScaleFromGenMC(Q, pt, eta, phi, nl, genPt, u1, s, m);
	}

	//if not, then:

	extern "C" double wrapper_kScaleAndSmearMC(int Q, double pt, double eta, double phi, int nl, double u1, double u2, int s, int m) {
	    //return rc.kScaleAndSmearMC(Q, pt, eta, phi, nl, u1, u2, s=0, m=0);
	    return rc.kScaleAndSmearMC(Q, pt, eta, phi, nl, u1, u2, s, m);
	}

	}

