#include "TVector3.h"


double v_angle(
	Float_t b1eta, Float_t b1phi,
	Float_t b2eta, Float_t b2phi
	)
        {
	TVector3 b1, b2;
	b1.SetPtEtaPhi(1, b1eta, b1phi);
	b2.SetPtEtaPhi(1, b2eta, b2phi);
	b1.SetMag(1);
	b2.SetMag(1);

	return b1.Angle(b2);
        }

double b_angle(
	Float_t b1x, Float_t b1y, Float_t b1z,
	Float_t b2x, Float_t b2y, Float_t b2z
	)
        {
	TVector3 b1, b2;
	b1.SetXYZ(b1x, b1y, b1z);
	b2.SetXYZ(b2x, b2y, b2z);
	return b1.Angle(b2);
        }

double b_anglesum(
	Float_t b1x, Float_t b1y, Float_t b1z,
	Float_t b2x, Float_t b2y, Float_t b2z,
	Float_t b3x, Float_t b3y, Float_t b3z
	)
        {
	TVector3 b1,b2,b3;

	b1.SetXYZ(b1x, b1y, b1z);
	b2.SetXYZ(b2x, b2y, b2z);
	b3.SetXYZ(b3x, b3y, b3z);

	return b1.Angle(b2) + b2.Angle(b3) + b3.Angle(b1);
        }

double b_angle_track(
	Float_t b1x, Float_t b1y, Float_t b1z,
	Float_t t1eta, Float_t t1phi
	)
	{
	TVector3 b, t;
	b.SetXYZ(b1x, b1y, b1z);

	t.SetPtEtaPhi(1, t1eta, t1phi);
	t.SetMag(1);

	return pow(t.Angle(b) - 1.57, 2); // with respect to 90degrees perpendicular
	}



// just for control
// this value is already stored in ntuple
double b_dist(
	Float_t b1x, Float_t b1y, Float_t b1z,
	Float_t b2x, Float_t b2y, Float_t b2z
	)
        {
	TVector3 b1, b2;
	b1.SetXYZ(b1x, b1y, b1z);
	b2.SetXYZ(b2x, b2y, b2z);
	double b_dist = (b1 - b2).Mag();
	return b_dist;
	}

double b_SV(
	Float_t b1x, Float_t b1y, Float_t b1z,
	Float_t b2x, Float_t b2y, Float_t b2z,
	Float_t v1eta, Float_t v1phi,
	Float_t v2eta, Float_t v2phi
	)
        {
	TVector3 b1, b2;
	b1.SetXYZ(b1x, b1y, b1z);
	b2.SetXYZ(b2x, b2y, b2z);
	double b_dist = (b1 - b2).Mag();

	TVector3 t1, t2;
	t1.SetPtEtaPhi(1, v1eta, v1phi);
	t2.SetPtEtaPhi(1, v2eta, v2phi);
	t1.SetMag(1);
	t2.SetMag(1);

	return b_dist / sin(t1.Angle(t2));
	}


TVector3 b_SV_1_xyz1(
	Float_t b1x, Float_t b1y, Float_t b1z,
	Float_t b2x, Float_t b2y, Float_t b2z,
	Float_t v1eta, Float_t v1phi,
	Float_t v2eta, Float_t v2phi
	)
        {
	TVector3 b1, b2;
	b1.SetXYZ(b1x, b1y, b1z);
	b2.SetXYZ(b2x, b2y, b2z);
	double b_dist = (b1 - b2).Mag();

	TVector3 t1, t2;
	t1.SetPtEtaPhi(1, v1eta, v1phi);
	t2.SetPtEtaPhi(1, v2eta, v2phi);
	t1.SetMag(1);
	t2.SetMag(1);

	double flight_length = b_dist / sin(t1.Angle(t2));

	t1.SetMag(flight_length);
	t2.SetMag(flight_length);

	return b1 + t1;
	}

TVector3 b_SV_1_xyz2(
	Float_t b1x, Float_t b1y, Float_t b1z,
	Float_t b2x, Float_t b2y, Float_t b2z,
	Float_t v1eta, Float_t v1phi,
	Float_t v2eta, Float_t v2phi
	)
        {
	TVector3 b1, b2;
	b1.SetXYZ(b1x, b1y, b1z);
	b2.SetXYZ(b2x, b2y, b2z);
	double b_dist = (b1 - b2).Mag();

	TVector3 t1, t2;
	t1.SetPtEtaPhi(1, v1eta, v1phi);
	t2.SetPtEtaPhi(1, v2eta, v2phi);
	t1.SetMag(1);
	t2.SetMag(1);

	double flight_length = b_dist / sin(t1.Angle(t2));

	t1.SetMag(flight_length);
	t2.SetMag(flight_length);

	return b2 + t2;
	}

TVector3 b_SV_2_xyz1(
	Float_t b1x, Float_t b1y, Float_t b1z,
	Float_t b2x, Float_t b2y, Float_t b2z,
	Float_t v1eta, Float_t v1phi,
	Float_t v2eta, Float_t v2phi
	)
        {
	TVector3 b1, b2;
	b1.SetXYZ(b1x, b1y, b1z);
	b2.SetXYZ(b2x, b2y, b2z);
	double b_dist = (b1 - b2).Mag();

	TVector3 t1, t2;
	t1.SetPtEtaPhi(1, v1eta, v1phi);
	t2.SetPtEtaPhi(1, v2eta, v2phi);
	t1.SetMag(1);
	t2.SetMag(1);

	double flight_length = b_dist / sin(t1.Angle(t2));

	t1.SetMag(flight_length);
	t2.SetMag(flight_length);

	return b1 - t1;
	}

TVector3 b_SV_2_xyz2(
	Float_t b1x, Float_t b1y, Float_t b1z,
	Float_t b2x, Float_t b2y, Float_t b2z,
	Float_t v1eta, Float_t v1phi,
	Float_t v2eta, Float_t v2phi
	)
        {
	TVector3 b1, b2;
	b1.SetXYZ(b1x, b1y, b1z);
	b2.SetXYZ(b2x, b2y, b2z);
	double b_dist = (b1 - b2).Mag();

	TVector3 t1, t2;
	t1.SetPtEtaPhi(1, v1eta, v1phi);
	t2.SetPtEtaPhi(1, v2eta, v2phi);
	t1.SetMag(1);
	t2.SetMag(1);

	double flight_length = b_dist / sin(t1.Angle(t2));

	t1.SetMag(flight_length);
	t2.SetMag(flight_length);

	return b2 - t2;
	}

double b_SV_1_dist(
	Float_t b1x, Float_t b1y, Float_t b1z,
	Float_t b2x, Float_t b2y, Float_t b2z,
	Float_t v1eta, Float_t v1phi,
	Float_t v2eta, Float_t v2phi
	)
        {
	TVector3 b1, b2;
	b1.SetXYZ(b1x, b1y, b1z);
	b2.SetXYZ(b2x, b2y, b2z);
	double b_dist = (b1 - b2).Mag();

	TVector3 t1, t2;
	t1.SetPtEtaPhi(1, v1eta, v1phi);
	t2.SetPtEtaPhi(1, v2eta, v2phi);
	t1.SetMag(1);
	t2.SetMag(1);

	double flight_length = b_dist / sin(t1.Angle(t2));

	t1.SetMag(flight_length);
	t2.SetMag(flight_length);

	return ((b1 + t1) - (b2 + t2)).Mag();
	}

double b_SV_2_dist(
	Float_t b1x, Float_t b1y, Float_t b1z,
	Float_t b2x, Float_t b2y, Float_t b2z,
	Float_t v1eta, Float_t v1phi,
	Float_t v2eta, Float_t v2phi
	)
        {
	TVector3 b1, b2;
	b1.SetXYZ(b1x, b1y, b1z);
	b2.SetXYZ(b2x, b2y, b2z);
	double b_dist = (b1 - b2).Mag();

	TVector3 t1, t2;
	t1.SetPtEtaPhi(1, v1eta, v1phi);
	t2.SetPtEtaPhi(1, v2eta, v2phi);
	t1.SetMag(1);
	t2.SetMag(1);

	double flight_length = b_dist / sin(t1.Angle(t2));

	t1.SetMag(flight_length);
	t2.SetMag(flight_length);

	return ((b1 - t1) - (b2 - t2)).Mag();
	}






