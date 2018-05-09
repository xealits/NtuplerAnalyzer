#include "RoccoR.cc"

RoccoR  rc("rcdata.2016.v3");

extern "C" double wrapper_kScaleDT(int Q, double pt, double eta, double phi, int s, int m) {
    return rc.kScaleDT(Q, pt, eta, phi, s, m);
}

extern "C" int wrapper_kScaleDT_double(double *res, int Q, double pt, double eta, double phi, int s, int m) {
    *res = rc.kScaleDT(Q, pt, eta, phi, s, m);
    return 55;
}


