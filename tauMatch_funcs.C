#include "TRandom3.h"
#include "TVector3.h"
#include "TVector2.h"
#include "TMatrixD.h"
#include "TMath.h"
#include "Math/SVector.h"
//#include "Math/SVector2.h"
#include "Math/SMatrix.h"


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

double b_Eta(
	Float_t b1x, Float_t b1y, Float_t b1z
	)
        {
	TVector3 b1;
	b1.SetXYZ(b1x, b1y, b1z);
	return b1.Eta();
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


double b_plane_unitary_volume(
	Float_t b1x, Float_t b1y, Float_t b1z,
	Float_t b2x, Float_t b2y, Float_t b2z,
	Float_t b3x, Float_t b3y, Float_t b3z
	)
        {
	TVector3 b1,b2,b3;

	b1.SetXYZ(b1x, b1y, b1z);
	b2.SetXYZ(b2x, b2y, b2z);
	b3.SetXYZ(b3x, b3y, b3z);

	b1.SetMag(1);
	b2.SetMag(1);
	b3.SetMag(1);

	TMatrixD m(3,3);
	m[0][0] = b1[0];
	m[0][1] = b1[1];
	m[0][2] = b1[2];
	m[1][0] = b2[0];
	m[1][1] = b2[1];
	m[1][2] = b2[2];
	m[2][0] = b3[0];
	m[2][1] = b3[1];
	m[2][2] = b3[2];

	return m.Determinant();
        }

double b_plane_angle_sum(
	Float_t b1x, Float_t b1y, Float_t b1z,
	Float_t b2x, Float_t b2y, Float_t b2z,
	Float_t b3x, Float_t b3y, Float_t b3z,
	Float_t v1eta, Float_t v1phi
	)
        {
	TVector3 b1,b2,b3, tau;

	b1.SetXYZ(b1x, b1y, b1z);
	b2.SetXYZ(b2x, b2y, b2z);
	b3.SetXYZ(b3x, b3y, b3z);
	tau.SetPtEtaPhi(1, v1eta, v1phi);

	//b1.SetMag(1);
	//b2.SetMag(1);
	//b3.SetMag(1);

	// find vector products of b _in_direction_ of tau
	// sum of them will give b-plane vector

	TVector3 b12 = b1.Cross(b2);
	TVector3 b23 = b2.Cross(b3);
	TVector3 b31 = b3.Cross(b1);

	if (b12 * tau < 0) b12 *= -1;
	if (b23 * tau < 0) b23 *= -1;
	if (b31 * tau < 0) b31 *= -1;

	return b31.Angle(b12) + b31.Angle(b23) + b12.Angle(b23);
        }

double b_plane_Eta(
	Float_t b1x, Float_t b1y, Float_t b1z,
	Float_t b2x, Float_t b2y, Float_t b2z,
	Float_t b3x, Float_t b3y, Float_t b3z,
	Float_t v1eta, Float_t v1phi
	)
        {
	TVector3 b1,b2,b3, tau;

	b1.SetXYZ(b1x, b1y, b1z);
	b2.SetXYZ(b2x, b2y, b2z);
	b3.SetXYZ(b3x, b3y, b3z);
	tau.SetPtEtaPhi(1, v1eta, v1phi);

	b1.SetMag(1);
	b2.SetMag(1);
	b3.SetMag(1);

	// find vector products of b _in_direction_ of tau
	// sum of them will give b-plane vector

	TVector3 b12 = b1.Cross(b2);
	TVector3 b23 = b2.Cross(b3);
	TVector3 b31 = b3.Cross(b1);

	if (b12 * tau < 0) b12 *= -1;
	if (b23 * tau < 0) b23 *= -1;
	if (b31 * tau < 0) b31 *= -1;

	TVector3 b_plane = b12 + b23 + b31;

	return b_plane.Eta();
        }

double b_plane_Phi(
	Float_t b1x, Float_t b1y, Float_t b1z,
	Float_t b2x, Float_t b2y, Float_t b2z,
	Float_t b3x, Float_t b3y, Float_t b3z,
	Float_t v1eta, Float_t v1phi
	)
        {
	TVector3 b1,b2,b3, tau;

	b1.SetXYZ(b1x, b1y, b1z);
	b2.SetXYZ(b2x, b2y, b2z);
	b3.SetXYZ(b3x, b3y, b3z);
	tau.SetPtEtaPhi(1, v1eta, v1phi);

	b1.SetMag(1);
	b2.SetMag(1);
	b3.SetMag(1);

	// find vector products of b _in_direction_ of tau
	// sum of them will give b-plane vector

	TVector3 b12 = b1.Cross(b2);
	TVector3 b23 = b2.Cross(b3);
	TVector3 b31 = b3.Cross(b1);

	if (b12 * tau < 0) b12 *= -1;
	if (b23 * tau < 0) b23 *= -1;
	if (b31 * tau < 0) b31 *= -1;

	TVector3 b_plane = b12 + b23 + b31;

	return b_plane.Phi();
        }




double b_plane_angle_tau(
	Float_t b1x, Float_t b1y, Float_t b1z,
	Float_t b2x, Float_t b2y, Float_t b2z,
	Float_t b3x, Float_t b3y, Float_t b3z,
	Float_t v1eta, Float_t v1phi
	)
        {
	TVector3 b1,b2,b3, tau;

	b1.SetXYZ(b1x, b1y, b1z);
	b2.SetXYZ(b2x, b2y, b2z);
	b3.SetXYZ(b3x, b3y, b3z);
	tau.SetPtEtaPhi(1, v1eta, v1phi);

	//b1.SetMag(1);
	//b2.SetMag(1);
	//b3.SetMag(1);

	// find vector products of b _in_direction_ of tau
	// sum of them will give b-plane vector

	TVector3 b12 = b1.Cross(b2);
	TVector3 b23 = b2.Cross(b3);
	TVector3 b31 = b3.Cross(b1);

	if (b12 * tau < 0) b12 *= -1;
	if (b23 * tau < 0) b23 *= -1;
	if (b31 * tau < 0) b31 *= -1;

	TVector3 b_plane = b12 + b23 + b31;

	return b_plane.Angle(tau);
        }

double b_plane_dPhi_tau(
	Float_t b1x, Float_t b1y, Float_t b1z,
	Float_t b2x, Float_t b2y, Float_t b2z,
	Float_t b3x, Float_t b3y, Float_t b3z,
	Float_t v1eta, Float_t v1phi
	)
        {
	TVector3 b1,b2,b3, tau;

	b1.SetXYZ(b1x, b1y, b1z);
	b2.SetXYZ(b2x, b2y, b2z);
	b3.SetXYZ(b3x, b3y, b3z);
	tau.SetPtEtaPhi(1, v1eta, v1phi);

	//b1.SetMag(1);
	//b2.SetMag(1);
	//b3.SetMag(1);

	// find vector products of b _in_direction_ of tau
	// sum of them will give b-plane vector

	TVector3 b12 = b1.Cross(b2);
	TVector3 b23 = b2.Cross(b3);
	TVector3 b31 = b3.Cross(b1);

	if (b12 * tau < 0) b12 *= -1;
	if (b23 * tau < 0) b23 *= -1;
	if (b31 * tau < 0) b31 *= -1;

	TVector3 b_plane = b12 + b23 + b31;

	return (b_plane - tau).Phi();
        }

double b_plane_dEta_tau(
	Float_t b1x, Float_t b1y, Float_t b1z,
	Float_t b2x, Float_t b2y, Float_t b2z,
	Float_t b3x, Float_t b3y, Float_t b3z,
	Float_t v1eta, Float_t v1phi
	)
        {
	TVector3 b1,b2,b3, tau;

	b1.SetXYZ(b1x, b1y, b1z);
	b2.SetXYZ(b2x, b2y, b2z);
	b3.SetXYZ(b3x, b3y, b3z);
	tau.SetPtEtaPhi(1, v1eta, v1phi);

	//b1.SetMag(1);
	//b2.SetMag(1);
	//b3.SetMag(1);

	// find vector products of b _in_direction_ of tau
	// sum of them will give b-plane vector

	TVector3 b12 = b1.Cross(b2);
	TVector3 b23 = b2.Cross(b3);
	TVector3 b31 = b3.Cross(b1);

	if (b12 * tau < 0) b12 *= -1;
	if (b23 * tau < 0) b23 *= -1;
	if (b31 * tau < 0) b31 *= -1;

	TVector3 b_plane = b12 + b23 + b31;

	return (b_plane - tau).Eta();
        }




// it seems tau-lep has smaller angles than lep-jets
double track_anglesum(
	Float_t v1eta, Float_t v1phi,
	Float_t v2eta, Float_t v2phi,
	Float_t v3eta, Float_t v3phi
	)
        {

	TVector3 t1, t2, t3;
	t1.SetPtEtaPhi(1, v1eta, v1phi);
	t2.SetPtEtaPhi(1, v2eta, v2phi);
	t3.SetPtEtaPhi(1, v3eta, v3phi);

	t1.SetMag(1);
	t2.SetMag(1);
	t3.SetMag(1);

	return t1.Angle(t2) + t2.Angle(t3) + t3.Angle(t1);
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


/*
 * b-s are not perpendicular to tracks
 * check if they lay in the plane of track
 */

// check the perpendicularity of perpendicular part of b
double b_perp_angle(
	Float_t b1x, Float_t b1y, Float_t b1z,
	Float_t v1eta, Float_t v1phi
	)
        {
	TVector3 b_vec;
	b_vec.SetXYZ(b1x, b1y, b1z);

	TVector3 t1;
	t1.SetPtEtaPhi(1, v1eta, v1phi);
	t1.SetMag(1);

	TVector3 b_long = t1 * (b_vec.Dot(t1));
	TVector3 b_perp = b_vec - b_long;

	return b_perp.Angle(t1);
	}

double b_tau_angle(
	Float_t b1x, Float_t b1y, Float_t b1z,
	Float_t v1eta, Float_t v1phi
	)
        {
	TVector3 b_vec, tau;
	b_vec.SetXYZ(b1x, b1y, b1z);
	tau.SetPtEtaPhi(1, v1eta, v1phi);

	return b_vec.Angle(tau);
	}

double b_tracksum_angle(
	Float_t b1x, Float_t b1y, Float_t b1z,
	Float_t v1pt, Float_t v1eta, Float_t v1phi,
	Float_t v2pt, Float_t v2eta, Float_t v2phi,
	Float_t v3pt, Float_t v3eta, Float_t v3phi
	)
        {
	TVector3 b_vec;
	b_vec.SetXYZ(b1x, b1y, b1z);

	TVector3 t1, t2, t3;
	t1.SetPtEtaPhi(v1pt, v1eta, v1phi);
	//t1.SetMag(1);
	t2.SetPtEtaPhi(v2pt, v2eta, v2phi);
	//t2.SetMag(1);
	t3.SetPtEtaPhi(v3pt, v3eta, v3phi);
	//t3.SetMag(1);

	TVector3 t_sum = t1 + t2 + t3;
	t_sum.SetMag(1);

	return b_vec.Angle(t_sum);
	}


double b_tracksum1_angle(
	Float_t b1x, Float_t b1y, Float_t b1z,
	Float_t v1eta, Float_t v1phi,
	Float_t v2eta, Float_t v2phi,
	Float_t v3eta, Float_t v3phi
	)
        {
	TVector3 b_vec;
	b_vec.SetXYZ(b1x, b1y, b1z);

	TVector3 t1, t2, t3;
	t1.SetPtEtaPhi(1, v1eta, v1phi);
	t1.SetMag(1);
	t2.SetPtEtaPhi(1, v2eta, v2phi);
	t2.SetMag(1);
	t3.SetPtEtaPhi(1, v3eta, v3phi);
	t3.SetMag(1);

	TVector3 t_sum = t1 + t2 + t3;
	t_sum.SetMag(1);

	return b_vec.Angle(t_sum);
	}

double b_plane_angle1(
	Float_t b1x, Float_t b1y, Float_t b1z,
	Float_t v1eta, Float_t v1phi,
	Float_t v2eta, Float_t v2phi,
	Float_t v3eta, Float_t v3phi
	)
        {
	TVector3 b_vec;
	b_vec.SetXYZ(b1x, b1y, b1z);

	TVector3 t1, t2, t3;
	t1.SetPtEtaPhi(1, v1eta, v1phi);
	t1.SetMag(1);
	t2.SetPtEtaPhi(1, v2eta, v2phi);
	t2.SetMag(1);
	t3.SetPtEtaPhi(1, v3eta, v3phi);
	t3.SetMag(1);

	TVector3 t_sum = t1 + t2 + t3;
	t_sum.SetMag(1);

	TVector3 b_long = t_sum * (b_vec.Dot(t_sum));
	TVector3 b_perp = b_vec - b_long;

	TVector3 t1_long = t_sum * (t1.Dot(t_sum));
	TVector3 t1_perp = t1 - t1_long;
	//TVector3 t2_long = t_sum * (t2.Dot(t_sum));
	//TVector3 t2_perp = t2 - t2_long;
	//TVector3 t3_long = t_sum * (t3.Dot(t_sum));
	//TVector3 t3_perp = t3 - t3_long;

	return b_perp.Angle(t1_perp);
	}

double b_plane_angle2(
	Float_t b1x, Float_t b1y, Float_t b1z,
	Float_t v1eta, Float_t v1phi,
	Float_t v2eta, Float_t v2phi,
	Float_t v3eta, Float_t v3phi
	)
        {
	TVector3 b_vec;
	b_vec.SetXYZ(b1x, b1y, b1z);

	TVector3 t1, t2, t3;
	t1.SetPtEtaPhi(1, v1eta, v1phi);
	t1.SetMag(1);
	t2.SetPtEtaPhi(1, v2eta, v2phi);
	t2.SetMag(1);
	t3.SetPtEtaPhi(1, v3eta, v3phi);
	t3.SetMag(1);

	TVector3 t_sum = t1 + t2 + t3;
	t_sum.SetMag(1);

	TVector3 b_long = t_sum * (b_vec.Dot(t_sum));
	TVector3 b_perp = b_vec - b_long;

	//TVector3 t1_long = t_sum * (t1.Dot(t_sum));
	//TVector3 t1_perp = t1 - t1_long;
	TVector3 t2_long = t_sum * (t2.Dot(t_sum));
	TVector3 t2_perp = t2 - t2_long;
	//TVector3 t3_long = t_sum * (t3.Dot(t_sum));
	//TVector3 t3_perp = t3 - t3_long;

	return b_perp.Angle(t2_perp);
	}

double b_plane_angle3(
	Float_t b1x, Float_t b1y, Float_t b1z,
	Float_t v1eta, Float_t v1phi,
	Float_t v2eta, Float_t v2phi,
	Float_t v3eta, Float_t v3phi
	)
        {
	TVector3 b_vec;
	b_vec.SetXYZ(b1x, b1y, b1z);

	TVector3 t1, t2, t3;
	t1.SetPtEtaPhi(1, v1eta, v1phi);
	t1.SetMag(1);
	t2.SetPtEtaPhi(1, v2eta, v2phi);
	t2.SetMag(1);
	t3.SetPtEtaPhi(1, v3eta, v3phi);
	t3.SetMag(1);

	TVector3 t_sum = t1 + t2 + t3;
	t_sum.SetMag(1);

	TVector3 b_long = t_sum * (b_vec.Dot(t_sum));
	TVector3 b_perp = b_vec - b_long;

	//TVector3 t1_long = t_sum * (t1.Dot(t_sum));
	//TVector3 t1_perp = t1 - t1_long;
	//TVector3 t2_long = t_sum * (t2.Dot(t_sum));
	//TVector3 t2_perp = t2 - t2_long;
	TVector3 t3_long = t_sum * (t3.Dot(t_sum));
	TVector3 t3_perp = t3 - t3_long;

	return b_perp.Angle(t3_perp);
	}


// testing second version of b projection
// correcting Z component of b, which has largest reconstruction error
// other components are not changed
// --- worse than projections
double corZ_angle_trk(
	Float_t b1x, Float_t b1y, Float_t b1z,
	Float_t v1eta, Float_t v1phi,
	Float_t v2eta, Float_t v2phi
	)
        {
	TVector3 b_vec, trk, tau;
	b_vec.SetXYZ(-b1x, -b1y, -b1z);
	trk.SetPtEtaPhi(1, v1eta, v1phi);
	tau.SetPtEtaPhi(1, v2eta, v2phi);

	tau.SetMag(1);

	// the transverse plane is built from tau vector
	// correct b
	// return angle between b and transverse part of track

	// corrected bz component
	double bz_cor = - (b_vec.x()*tau.x() + b_vec.y()*tau.y()) / tau.z();
	b_vec.SetZ(bz_cor);

	// find perpendicular of track
	TVector3 trk_long = tau * (trk.Dot(tau));
	TVector3 trk_perp = trk - trk_long;

	return b_vec.Angle(trk_perp);
	//return trk_perp.Angle(tau);
	}

// let's try anyway
double corZ_average_SV(
	Float_t b1x, Float_t b1y, Float_t b1z,
	Float_t b2x, Float_t b2y, Float_t b2z,
	Float_t b3x, Float_t b3y, Float_t b3z,
	Float_t v1pt, Float_t v1eta, Float_t v1phi,
	Float_t v2pt, Float_t v2eta, Float_t v2phi,
	Float_t v3pt, Float_t v3eta, Float_t v3phi
	)
        {
	Float_t tracker_error = 0.2; // approximately systematic error on positions

	TVector3 b_vec1, b_vec2, b_vec3;
	b_vec1.SetXYZ(-b1x, -b1y, -b1z);
	b_vec2.SetXYZ(-b2x, -b2y, -b2z);
	b_vec3.SetXYZ(-b3x, -b3y, -b3z);

	TVector3 t1, t2, t3;
	t1.SetPtEtaPhi(v1pt, v1eta, v1phi);
	t1.SetMag(1);
	t2.SetPtEtaPhi(v2pt, v2eta, v2phi);
	t2.SetMag(1);
	t3.SetPtEtaPhi(v3pt, v3eta, v3phi);
	t3.SetMag(1);

	TVector3 t_sum = t1 + t2 + t3;
	t_sum.SetMag(1);

	// after establishing direction of tau
	// tracks are only geometrical lines
	t1.SetMag(1);
	t2.SetMag(1);
	t3.SetMag(1);

	// instead of projections correct Z
	//TVector3 b_long1 = t_sum * (b_vec1.Dot(t_sum));
	//TVector3 b_perp1 = b_vec1 - b_long1;
	//TVector3 b_long2 = t_sum * (b_vec2.Dot(t_sum));
	//TVector3 b_perp2 = b_vec2 - b_long2;
	//TVector3 b_long3 = t_sum * (b_vec3.Dot(t_sum));
	//TVector3 b_perp3 = b_vec3 - b_long3;
	double bz_cor = - (b_vec1.x()*t_sum.x() + b_vec1.y()*t_sum.y()) / t_sum.z();
	b_vec1.SetZ(bz_cor);
	bz_cor = - (b_vec2.x()*t_sum.x() + b_vec2.y()*t_sum.y()) / t_sum.z();
	b_vec2.SetZ(bz_cor);
	bz_cor = - (b_vec3.x()*t_sum.x() + b_vec3.y()*t_sum.y()) / t_sum.z();
	b_vec3.SetZ(bz_cor);

	// in the perp plane find b long to tracks

	// perpendicular parts of tracks
	TVector3 t1_long = t_sum * (t1.Dot(t_sum));
	TVector3 t1_perp = t1 - t1_long;
	TVector3 t2_long = t_sum * (t2.Dot(t_sum));
	TVector3 t2_perp = t2 - t2_long;
	TVector3 t3_long = t_sum * (t3.Dot(t_sum));
	TVector3 t3_perp = t3 - t3_long;

	t1_perp.SetMag(1);
	t2_perp.SetMag(1);
	t3_perp.SetMag(1);

	//TVector3 perp_b_long1 = t1_perp * (b_perp1.Dot(t1_perp));
	////TVector3 perp_b_perp1 = b_perp1 - perp_b_long1;
	//TVector3 perp_b_long2 = t2_perp * (b_perp2.Dot(t2_perp));
	////TVector3 perp_b_perp2 = b_perp2 - perp_b_long2;
	//TVector3 perp_b_long3 = t3_perp * (b_perp3.Dot(t3_perp));
	////TVector3 perp_b_perp3 = b_perp3 - perp_b_long3;

	// just using corrected b vecs
	TVector3 perp_b_long1 = t1_perp * (b_vec1.Dot(t1_perp));
	TVector3 perp_b_long2 = t2_perp * (b_vec2.Dot(t2_perp));
	TVector3 perp_b_long3 = t3_perp * (b_vec3.Dot(t3_perp));


	// the best point calculation
	TVector3 dV = t1 - t2;
	TVector3 dB = perp_b_long1 - perp_b_long2;
	double x12 = dV.Dot(dB) / dV.Mag2();
	dV = t2 - t3;
	dB = perp_b_long2 - perp_b_long3;
	double x23 = dV.Dot(dB) / dV.Mag2();
	dV = t3 - t1;
	dB = perp_b_long3 - perp_b_long1;
	double x31 = dV.Dot(dB) / dV.Mag2();

	// and systematic error of tracker
	double syst12 = tracker_error / t1.Angle(t2); // technically / Sin (or Tan), but Sin = Angle with these angles
	double syst23 = tracker_error / t2.Angle(t3); // technically / Sin (or Tan), but Sin = Angle with these angles
	double syst31 = tracker_error / t3.Angle(t1); // technically / Sin (or Tan), but Sin = Angle with these angles
	//double syst = pow(syst12, 2) + pow(syst23, 2) + pow(syst31, 2);
	double syst12_weight = 1/syst12;
	double syst23_weight = 1/syst23;
	double syst31_weight = 1/syst31;

	//double x_average = (x12 + x23 + x31) / 3;
	//double x_deviation = pow(x12 - x_average, 2) + pow(x23 - x_average, 2) + pow(x31 - x_average, 2);
	//double x_dev_syst = x_deviation + syst;

	// weighted average with tracker errors
	double x_average = (x12*syst12_weight + x23*syst23_weight + x31*syst31_weight) / (syst12_weight + syst23_weight + syst31_weight);
	double x_deviation = (syst12_weight*pow(x12 - x_average, 2) + syst23_weight*pow(x23 - x_average, 2) + syst31_weight*pow(x31 - x_average, 2))/(2*(syst12_weight + syst23_weight + syst31_weight)/3);

	return x_average / sqrt(x_deviation);
	}



// testing third version, aiming at optimal direction
// avarage of b-track plane intersections
double b_track_plane_intersections_average(
	Float_t b1x, Float_t b1y, Float_t b1z,
	Float_t b2x, Float_t b2y, Float_t b2z,
	Float_t b3x, Float_t b3y, Float_t b3z,
	Float_t v1pt, Float_t v1eta, Float_t v1phi,
	Float_t v2pt, Float_t v2eta, Float_t v2phi,
	Float_t v3pt, Float_t v3eta, Float_t v3phi,
	Float_t taupt, Float_t taueta, Float_t tauphi
	)
        {
	TVector3 b_vec1, b_vec2, b_vec3;
	b_vec1.SetXYZ(b1x, b1y, b1z); // 100% known that z here has giant error -- need to do something with it
	b_vec2.SetXYZ(b2x, b2y, b2z);
	b_vec3.SetXYZ(b3x, b3y, b3z);

	TVector3 t1, t2, t3;
	t1.SetPtEtaPhi(v1pt, v1eta, v1phi);
	t2.SetPtEtaPhi(v2pt, v2eta, v2phi);
	t3.SetPtEtaPhi(v3pt, v3eta, v3phi);

	//TVector3 tau = t1+t2+t3;
	TVector3 tau;
	tau.SetPtEtaPhi(taupt, taueta, tauphi);

	TVector3 bpl1 = b_vec1.Cross(t1);
	TVector3 bpl2 = b_vec2.Cross(t2);
	TVector3 bpl3 = b_vec3.Cross(t3);

	TVector3 intr12 = bpl1.Cross(bpl2);
	TVector3 intr23 = bpl2.Cross(bpl3);
	TVector3 intr31 = bpl3.Cross(bpl1);

	// point them at tau
	if (intr12 * tau < 0) intr12 *= -1;
	if (intr23 * tau < 0) intr23 *= -1;
	if (intr31 * tau < 0) intr31 *= -1;

	//intr12.SetMag(1);
	//intr23.SetMag(1);
	//intr31.SetMag(1);

	TVector3 average = intr12 + intr23 + intr31;

	//return intr12.Angle(intr23) + intr23.Angle(intr31) + intr31.Angle(intr12);
	return average.Angle(tau);
	}


// check angle between perpendicular b and track
double b_track_plane_intersections_average_perp_angle(
	int which, //
	Float_t b1x, Float_t b1y, Float_t b1z,
	Float_t b2x, Float_t b2y, Float_t b2z,
	Float_t b3x, Float_t b3y, Float_t b3z,
	Float_t v1pt, Float_t v1eta, Float_t v1phi,
	Float_t v2pt, Float_t v2eta, Float_t v2phi,
	Float_t v3pt, Float_t v3eta, Float_t v3phi,
	Float_t taupt, Float_t taueta, Float_t tauphi
	)
        {
	TVector3 b_vec1, b_vec2, b_vec3;
	b_vec1.SetXYZ(b1x, b1y, b1z); // 100% known that z here has giant error -- need to do something with it
	b_vec2.SetXYZ(b2x, b2y, b2z);
	b_vec3.SetXYZ(b3x, b3y, b3z);

	TVector3 t1, t2, t3;
	t1.SetPtEtaPhi(v1pt, v1eta, v1phi);
	t2.SetPtEtaPhi(v2pt, v2eta, v2phi);
	t3.SetPtEtaPhi(v3pt, v3eta, v3phi);

	//TVector3 tau = t1+t2+t3;
	TVector3 tau;
	tau.SetPtEtaPhi(taupt, taueta, tauphi);

	TVector3 bpl1 = b_vec1.Cross(t1);
	TVector3 bpl2 = b_vec2.Cross(t2);
	TVector3 bpl3 = b_vec3.Cross(t3);

	TVector3 intr12 = bpl1.Cross(bpl2);
	TVector3 intr23 = bpl2.Cross(bpl3);
	TVector3 intr31 = bpl3.Cross(bpl1);

	// point them at tau
	if (intr12 * tau < 0) intr12 *= -1;
	if (intr23 * tau < 0) intr23 *= -1;
	if (intr31 * tau < 0) intr31 *= -1;

	//intr12.SetMag(1);
	//intr23.SetMag(1);
	//intr31.SetMag(1);

	TVector3 average = intr12 + intr23 + intr31;
	average.SetMag(1);

	// find perpendicular b-s
	TVector3 b_long, b_perp, t_long, t_perp;

	switch (which)
		{
		case 0:
			b_long = average * (b_vec1.Dot(average));
			b_perp = b_vec1 - b_long;
			t_long = average * (t1.Dot(average));
			t_perp = t1 - t_long;
			return b_perp.Angle(t_perp);
			//return b_perp.Angle(average);
			//return t_perp.Angle(average);
			break;
		case 1:
			b_long = average * (b_vec2.Dot(average));
			b_perp = b_vec2 - b_long;
			t_long = average * (t2.Dot(average));
			t_perp = t2 - t_long;
			return b_perp.Angle(t_perp);
			//return average.Angle(t_perp);
			//return average.Angle(b_perp);
			break;
		case 2:
			b_long = average * (b_vec3.Dot(average));
			b_perp = b_vec3 - b_long;
			t_long = average * (t3.Dot(average));
			t_perp = t3 - t_long;
			return b_perp.Angle(t_perp);
			break;
		case 3:
			// for tests
			b_long = average * (b_vec3.Dot(average));
			b_perp = b_vec3 - b_long;
			t_long = average * (t3.Dot(average));
			t_perp = t3 - t_long;
			return b_perp.Angle(t_perp);
			break;
		default:
			return average.Angle(tau);
		}
	return -1;
	}
// strange result
// doesn't break with typos


double max_track_bis_angle(
	int which,
	Float_t v1pt, Float_t v1eta, Float_t v1phi,
	Float_t v2pt, Float_t v2eta, Float_t v2phi,
	Float_t v3pt, Float_t v3eta, Float_t v3phi
	)
	{
	TVector3 t1, t2, t3;
	t1.SetPtEtaPhi(v1pt, v1eta, v1phi);
	t2.SetPtEtaPhi(v2pt, v2eta, v2phi);
	t3.SetPtEtaPhi(v3pt, v3eta, v3phi);

	if (which == 0) // geometrical bis
		{
		t1.SetMag(1);
		t2.SetMag(1);
		t3.SetMag(1);
		}

	TVector3 bis = t1+t2+t3;

	double an1 = t1.Angle(bis);
	double an2 = t2.Angle(bis);
	double an3 = t3.Angle(bis);

	double max_an = (an1 > an2? an1 : an2);
	if (an3>max_an) max_an = an3;

	return max_an;
	}


TRandom3 *r3 = new TRandom3();

// optimal directions
// ok, it's late, let's try simple SV with this plane direction
double optimal_directions_intersections_raw_SV(
	int which, // for tests
	Float_t b1x, Float_t b1y, Float_t b1z,
	Float_t b2x, Float_t b2y, Float_t b2z,
	Float_t b3x, Float_t b3y, Float_t b3z,
	Float_t v1pt, Float_t v1eta, Float_t v1phi,
	Float_t v2pt, Float_t v2eta, Float_t v2phi,
	Float_t v3pt, Float_t v3eta, Float_t v3phi,
	Float_t taupt, Float_t taueta, Float_t tauphi
	)
        {
	Float_t tracker_error = 0.002; // approximately systematic error on positions
	// it will cancel out with weights

	TVector3 b_vec1, b_vec2, b_vec3;
	b_vec1.SetXYZ(-b1x, -b1y, -b1z); // 100% known that z here has giant error -- need to do something with it
	b_vec2.SetXYZ(-b2x, -b2y, -b2z);
	b_vec3.SetXYZ(-b3x, -b3y, -b3z);

	TVector3 t1, t2, t3;
	t1.SetPtEtaPhi(v1pt, v1eta, v1phi);
	t2.SetPtEtaPhi(v2pt, v2eta, v2phi);
	t3.SetPtEtaPhi(v3pt, v3eta, v3phi);

	// weighted bis direction -- used in simple b SV (Friday result)
	TVector3 t_sum = t1 + t2 + t3;
	t_sum.SetMag(1);

	// after establishing direction of tau
	// tracks are only geometrical lines
	t1.SetMag(1);
	t2.SetMag(1);
	t3.SetMag(1);

	//TVector3 tau = t1+t2+t3;
	TVector3 tau;
	tau.SetPtEtaPhi(taupt, taueta, tauphi);

	// find the "optimal direction"
	// -- direction of minimal angles betwee tracks and b-s in perpendicular plane
	// it should maximize sum of angles between tracks and b-s in transverse plane
	// and minimize the angles between intersections
	double max_angle_sum = 0;
	double min_intersection_angle_sum = 999;
	// randomly change Z coordinate of b-s to find the maximum
	TVector3 max_average;

	double orig_b1_z = b_vec1.Z();
	double orig_b2_z = b_vec2.Z();
	double orig_b3_z = b_vec3.Z();
	double max_b1_z = orig_b1_z;
	double max_b2_z = orig_b2_z;
	double max_b3_z = orig_b3_z;

	for (unsigned int i = 0; i<1000; i++)
		{
		//b_vec1.SetZ(2 * orig_b1_z * r3->Uniform());
		//b_vec2.SetZ(2 * orig_b2_z * r3->Uniform());
		//b_vec3.SetZ(2 * orig_b3_z * r3->Uniform());
		//   Z +- 1.5Z = -0.5Z -- + 2.5Z = r(3Z) - 0.5Z
		b_vec1.SetZ(3 * orig_b1_z * r3->Uniform() - 0.5*orig_b1_z);
		b_vec2.SetZ(3 * orig_b2_z * r3->Uniform() - 0.5*orig_b2_z);
		b_vec3.SetZ(3 * orig_b3_z * r3->Uniform() - 0.5*orig_b3_z);

		TVector3 bpl1 = b_vec1.Cross(t1);
		TVector3 bpl2 = b_vec2.Cross(t2);
		TVector3 bpl3 = b_vec3.Cross(t3);

		TVector3 intr12 = bpl1.Cross(bpl2);
		TVector3 intr23 = bpl2.Cross(bpl3);
		TVector3 intr31 = bpl3.Cross(bpl1);

		// point them at tau
		if (intr12 * tau < 0) intr12 *= -1;
		if (intr23 * tau < 0) intr23 *= -1;
		if (intr31 * tau < 0) intr31 *= -1;

		// to discard effect from randomly big Z affecting stuff
		// set directions to 1
		intr12.SetMag(1);
		intr23.SetMag(1);
		intr31.SetMag(1);

		double intersection_angle_sum = intr12.Angle(intr23) + intr23.Angle(intr31) + intr31.Angle(intr12);

		// so this is the optimal direction
		TVector3 average = intr12 + intr23 + intr31;
		average.SetMag(1);
		// got plane direction

		// and to optimal direction
		// find perpendicular b-s
		TVector3 b_long1 = average * (b_vec1.Dot(average));
		TVector3 b_perp1 = b_vec1 - b_long1;
		TVector3 b_long2 = average * (b_vec2.Dot(average));
		TVector3 b_perp2 = b_vec2 - b_long2;
		TVector3 b_long3 = average * (b_vec3.Dot(average));
		TVector3 b_perp3 = b_vec3 - b_long3;

		// perpendicular parts of tracks
		TVector3 t1_long = average * (t1.Dot(average));
		TVector3 t1_perp = t1 - t1_long;
		TVector3 t2_long = average * (t2.Dot(average));
		TVector3 t2_perp = t2 - t2_long;
		TVector3 t3_long = average * (t3.Dot(average));
		TVector3 t3_perp = t3 - t3_long;

		double angle_sum = b_perp1.Angle(t1_perp) + b_perp2.Angle(t2_perp) + b_perp3.Angle(t3_perp);
		if (angle_sum > max_angle_sum && intersection_angle_sum < min_intersection_angle_sum)
			{
			min_intersection_angle_sum = intersection_angle_sum;
			max_angle_sum = angle_sum;
			max_average = average;
			max_b1_z = b_vec1.Z();
			max_b2_z = b_vec2.Z();
			max_b3_z = b_vec3.Z();
			}
		}
	b_vec1.SetZ(max_b1_z);
	b_vec2.SetZ(max_b2_z);
	b_vec3.SetZ(max_b3_z);

	// and to optimal direction
	// find perpendicular b-s
	TVector3 b_long1 = max_average * (b_vec1.Dot(max_average));
	TVector3 b_perp1 = b_vec1 - b_long1;
	TVector3 b_long2 = max_average * (b_vec2.Dot(max_average));
	TVector3 b_perp2 = b_vec2 - b_long2;
	TVector3 b_long3 = max_average * (b_vec3.Dot(max_average));
	TVector3 b_perp3 = b_vec3 - b_long3;

	// perpendicular parts of tracks
	TVector3 t1_long = max_average * (t1.Dot(max_average));
	TVector3 t1_perp = t1 - t1_long;
	TVector3 t2_long = max_average * (t2.Dot(max_average));
	TVector3 t2_perp = t2 - t2_long;
	TVector3 t3_long = max_average * (t3.Dot(max_average));
	TVector3 t3_perp = t3 - t3_long;

	// project found b-s to perp tracks
	// in principle it should not be needed, since the direction is found to fit them together well
	// but let's try to get to simple SV best result
	t1_perp.SetMag(1);
	t2_perp.SetMag(1);
	t3_perp.SetMag(1);

	TVector3 b_long_perp1 = t1_perp * (b_perp1.Dot(t1_perp));
	TVector3 b_long_perp2 = t2_perp * (b_perp2.Dot(t2_perp));
	TVector3 b_long_perp3 = t3_perp * (b_perp3.Dot(t3_perp));



	// perpendiculars to bis direction, for reference
	// find perpendicular b-s
	TVector3 b_bis_long1 = t_sum * (b_vec1.Dot(t_sum));
	TVector3 b_bis_perp1 = b_vec1 - b_bis_long1;
	TVector3 b_bis_long2 = t_sum * (b_vec2.Dot(t_sum));
	TVector3 b_bis_perp2 = b_vec2 - b_bis_long2;
	TVector3 b_bis_long3 = t_sum * (b_vec3.Dot(t_sum));
	TVector3 b_bis_perp3 = b_vec3 - b_bis_long3;

	// perpendicular parts of tracks
	TVector3 t1_bis_long = t_sum * (t1.Dot(t_sum));
	TVector3 t1_bis_perp = t1 - t1_bis_long;
	TVector3 t2_bis_long = t_sum * (t2.Dot(t_sum));
	TVector3 t2_bis_perp = t2 - t2_bis_long;
	TVector3 t3_bis_long = t_sum * (t3.Dot(t_sum));
	TVector3 t3_bis_perp = t3 - t3_bis_long;

	// in the perp plane find b long to tracks
	// -- nope, no additional correction to b-s

	// the best point calculation
	// with just transverse b-s
	//TVector3 dV = t1 - t2;
	//TVector3 dB = b_perp1 - b_perp2;
	//double x12 = dV.Dot(dB) / dV.Mag2();
	//dV = t2 - t3;
	//dB = b_perp2 - b_perp3;
	//double x23 = dV.Dot(dB) / dV.Mag2();
	//dV = t3 - t1;
	//dB = b_perp3 - b_perp1;
	//double x31 = dV.Dot(dB) / dV.Mag2();

	// the best point calculation with projected b-s
	TVector3 dV = t1 - t2;
	TVector3 dB = b_long_perp1 - b_long_perp2;
	double x12 = - dV.Dot(dB) / dV.Mag2();
	dV = t2 - t3;
	dB = b_long_perp2 - b_long_perp3;
	double x23 = - dV.Dot(dB) / dV.Mag2();
	dV = t3 - t1;
	dB = b_long_perp3 - b_long_perp1;
	double x31 = - dV.Dot(dB) / dV.Mag2();

	TVector3 bp12 = b_long_perp1 + x12 * t1;
	TVector3 bp21 = b_long_perp2 + x12 * t2;
	TVector3 bp_1 = 0.5*(bp12 + bp21);

	TVector3 bp23 = b_long_perp2 + x23 * t2;
	TVector3 bp32 = b_long_perp3 + x23 * t3;
	TVector3 bp_2 = 0.5*(bp23 + bp32);

	TVector3 bp31 = b_long_perp3 + x31 * t3;
	TVector3 bp13 = b_long_perp1 + x31 * t1;
	TVector3 bp_3 = 0.5*(bp31 + bp13);

	TVector3 bp_average = 0.3333*(bp_1 + bp_2 + bp_3);
	TVector3 bp_dev1 = bp_1 - bp_average;
	TVector3 bp_dev2 = bp_2 - bp_average;
	TVector3 bp_dev3 = bp_3 - bp_average;

	// and systematic error of tracker
	double syst12 = tracker_error / t1.Angle(t2); // technically / Sin (or Tan), but Sin = Angle with these angles
	double syst23 = tracker_error / t2.Angle(t3); // technically / Sin (or Tan), but Sin = Angle with these angles
	double syst31 = tracker_error / t3.Angle(t1); // technically / Sin (or Tan), but Sin = Angle with these angles
	//double syst = pow(syst12, 2) + pow(syst23, 2) + pow(syst31, 2);
	double syst12_weight = 1/syst12;
	double syst23_weight = 1/syst23;
	double syst31_weight = 1/syst31;

	// not weighted averages
	double x_average = (x12 + x23 + x31) / 3;
	double x_deviation = (pow(x12 - x_average, 2) + pow(x23 - x_average, 2) + pow(x31 - x_average, 2)) * 0.3333;
	//double x_dev_syst = x_deviation + syst;

	// weighted average with tracker errors
	//double x_average = (x12*syst12_weight + x23*syst23_weight + x31*syst31_weight) / (syst12_weight + syst23_weight + syst31_weight);
	//double x_deviation = (syst12_weight*pow(x12 - x_average, 2) + syst23_weight*pow(x23 - x_average, 2) + syst31_weight*pow(x31 - x_average, 2))/(2*(syst12_weight + syst23_weight + syst31_weight)/3);

	switch (which)
		{
		// test optimal direction
		case 1:
			return max_average.Eta();
		case 2:
			return max_average.Angle(tau);

		// test on optimalness of optimal direction
		case 3: // and same angles w.r. to bis direction
			return b_bis_perp1.Angle(t1_bis_perp) + b_bis_perp2.Angle(t2_bis_perp) + b_bis_perp3.Angle(t3_bis_perp);
		case 4: // angles between b and tracks in perp plane of the direction
			return b_perp1.Angle(t1_perp) + b_perp2.Angle(t2_perp) + b_perp3.Angle(t3_perp);
		case 5:
			return max_angle_sum;
		case 6:
			return min_intersection_angle_sum;

		// best point tests
		case 7:
			return (bp12 - bp21).Mag();
		case 8:
			return (bp23 - bp32).Mag();
		case 9:
			return (bp31 - bp13).Mag();
		// distances between b-s of tracks
		case 10:
			return (b_long_perp1 - b_long_perp2).Mag();
		case 11:
			return (b_long_perp2 - b_long_perp3).Mag();
		case 12:
			return (b_long_perp3 - b_long_perp1).Mag();

		// flight Len, SV convergence, deviations
		case 13: // flight Len
			if (bp_average.Dot(tau) > 0)
				return bp_average.Mag();
			else
				return - bp_average.Mag();
		case 14: // convergence
			return sqrt(pow((bp12 - bp21).Mag(), 2) + pow((bp23 - bp32).Mag(), 2) + pow((bp31 - bp13).Mag(), 2));
		case 15: // deviations
			return sqrt(pow(bp_dev1.Mag(), 2) + pow(bp_dev2.Mag(), 2) + pow(bp_dev3.Mag(), 2));
		case 16: // flight Sign with deviation
			if (bp_average.Dot(tau) > 0)
				return bp_average.Mag() / sqrt(pow(bp_dev1.Mag(), 2) + pow(bp_dev2.Mag(), 2) + pow(bp_dev3.Mag(), 2));
			else
				return - bp_average.Mag() / sqrt(pow(bp_dev1.Mag(), 2) + pow(bp_dev2.Mag(), 2) + pow(bp_dev3.Mag(), 2));
		}

	return x_average / sqrt(x_deviation);
	}





// simple SV
// finds best point for each pair
// with the following simple algorithm

double b_simple_bestpoint(
	Float_t b1x, Float_t b1y, Float_t b1z,
	Float_t b2x, Float_t b2y, Float_t b2z,
	Float_t v1eta, Float_t v1phi,
	Float_t v2eta, Float_t v2phi,
	Float_t v3eta, Float_t v3phi
	)
        {
	TVector3 b_vec1, b_vec2;
	b_vec1.SetXYZ(b1x, b1y, b1z);
	b_vec2.SetXYZ(b2x, b2y, b2z);

	TVector3 t1, t2, t3;
	t1.SetPtEtaPhi(1, v1eta, v1phi);
	t1.SetMag(1);
	t2.SetPtEtaPhi(1, v2eta, v2phi);
	t2.SetMag(1);
	t3.SetPtEtaPhi(1, v3eta, v3phi);
	t3.SetMag(1);

	// sum of track vectors -- roughly direction of tau
	TVector3 t_sum = t1 + t2 + t3;
	t_sum.SetMag(1);

	// find perpendicular b-s
	TVector3 b_long1 = t_sum * (b_vec1.Dot(t_sum));
	TVector3 b_perp1 = b_vec1 - b_long1;

	TVector3 b_long2 = t_sum * (b_vec2.Dot(t_sum));
	TVector3 b_perp2 = b_vec2 - b_long2;

	// in the perp plane project b perp on tracks

	// perpendicular parts of tracks
	TVector3 t1_long = t_sum * (t1.Dot(t_sum));
	TVector3 t1_perp = t1 - t1_long;
	TVector3 t2_long = t_sum * (t2.Dot(t_sum));
	TVector3 t2_perp = t2 - t2_long;
	//TVector3 t3_long = t_sum * (t3.Dot(t_sum));
	//TVector3 t3_perp = t3 - t3_long;

	b_long1 = t1_perp * (b_perp1.Dot(t1_perp));
	b_perp1 = b_perp1 - b_long1;

	b_long2 = t2_perp * (b_perp2.Dot(t2_perp));
	b_perp2 = b_perp2 - b_long2;

	// the best point calculation
	TVector3 dV = t1 - t2;
	TVector3 dB = b_perp1 - b_perp2;
	return dV.Dot(dB) / dV.Mag2();
	}


double b_simple_SV(
	Float_t b1x, Float_t b1y, Float_t b1z,
	Float_t b2x, Float_t b2y, Float_t b2z,
	Float_t b3x, Float_t b3y, Float_t b3z,
	Float_t v1eta, Float_t v1phi,
	Float_t v2eta, Float_t v2phi,
	Float_t v3eta, Float_t v3phi
	)
        {
	TVector3 b_vec1, b_vec2, b_vec3;
	b_vec1.SetXYZ(b1x, b1y, b1z);
	b_vec2.SetXYZ(b2x, b2y, b2z);
	b_vec3.SetXYZ(b3x, b3y, b3z);

	TVector3 t1, t2, t3;
	t1.SetPtEtaPhi(1, v1eta, v1phi);
	t1.SetMag(1);
	t2.SetPtEtaPhi(1, v2eta, v2phi);
	t2.SetMag(1);
	t3.SetPtEtaPhi(1, v3eta, v3phi);
	t3.SetMag(1);

	TVector3 t_sum = t1 + t2 + t3;
	t_sum.SetMag(1);

	TVector3 b_long1 = t_sum * (b_vec1.Dot(t_sum));
	TVector3 b_perp1 = b_vec1 - b_long1;

	TVector3 b_long2 = t_sum * (b_vec2.Dot(t_sum));
	TVector3 b_perp2 = b_vec2 - b_long2;

	TVector3 b_long3 = t_sum * (b_vec3.Dot(t_sum));
	TVector3 b_perp3 = b_vec3 - b_long3;

	// in the perp plane project b on tracks

	// perpendicular parts of tracks
	TVector3 t1_long = t_sum * (t1.Dot(t_sum));
	TVector3 t1_perp = t1 - t1_long;
	TVector3 t2_long = t_sum * (t2.Dot(t_sum));
	TVector3 t2_perp = t2 - t2_long;
	TVector3 t3_long = t_sum * (t3.Dot(t_sum));
	TVector3 t3_perp = t3 - t3_long;

	t1_perp.SetMag(1);
	t2_perp.SetMag(1);
	t3_perp.SetMag(1);

	b_long1 = t1_perp * (b_perp1.Dot(t1_perp));
	//b_perp1 = b_perp1 - b_long1;

	b_long2 = t2_perp * (b_perp2.Dot(t2_perp));
	//b_perp2 = b_perp2 - b_long2;

	b_long3 = t3_perp * (b_perp3.Dot(t3_perp));
	//b_perp3 = b_perp3 - b_long3;


	// the best point calculation
	TVector3 dV = t1 - t2;
	TVector3 dB = b_long1 - b_long2;
	double x12 = dV.Dot(dB) / dV.Mag2();
	dV = t2 - t3;
	dB = b_long2 - b_long3;
	double x23 = dV.Dot(dB) / dV.Mag2();
	dV = t3 - t1;
	dB = b_long3 - b_long1;
	double x31 = dV.Dot(dB) / dV.Mag2();

	double x_average = (x12 + x23 + x31) / 3;
	double x_deviation = sqrt(pow(x12 - x_average, 2) + pow(x23 - x_average, 2) + pow(x31 - x_average, 2));

	return x_average / x_deviation;
	}

// FINAL SV calc
// add pt-s to tracks
double b_simple_SV_pt(
	Float_t b1x, Float_t b1y, Float_t b1z,
	Float_t b2x, Float_t b2y, Float_t b2z,
	Float_t b3x, Float_t b3y, Float_t b3z,
	Float_t v1pt, Float_t v1eta, Float_t v1phi,
	Float_t v2pt, Float_t v2eta, Float_t v2phi,
	Float_t v3pt, Float_t v3eta, Float_t v3phi
	)
        {
	Float_t tracker_error = 0.2; // approximately systematic error on positions
	// 0.0002 is not noticable
	// 0.002 has slight effect
	// cannot pass it as function parameter in root/cint/tformula!!!!!!

	TVector3 b_vec1, b_vec2, b_vec3;
	b_vec1.SetXYZ(b1x, b1y, b1z);
	b_vec2.SetXYZ(b2x, b2y, b2z);
	b_vec3.SetXYZ(b3x, b3y, b3z);

	TVector3 t1, t2, t3;
	t1.SetPtEtaPhi(v1pt, v1eta, v1phi);
	t2.SetPtEtaPhi(v2pt, v2eta, v2phi);
	t3.SetPtEtaPhi(v3pt, v3eta, v3phi);

	TVector3 t_sum = t1 + t2 + t3;
	t_sum.SetMag(1);

	// after establishing direction of tau
	// tracks are only geometrical lines
	t1.SetMag(1);
	t2.SetMag(1);
	t3.SetMag(1);

	TVector3 b_long1 = t_sum * (b_vec1.Dot(t_sum));
	TVector3 b_perp1 = b_vec1 - b_long1;

	TVector3 b_long2 = t_sum * (b_vec2.Dot(t_sum));
	TVector3 b_perp2 = b_vec2 - b_long2;

	TVector3 b_long3 = t_sum * (b_vec3.Dot(t_sum));
	TVector3 b_perp3 = b_vec3 - b_long3;

	// in the perp plane find b long to tracks

	// perpendicular parts of tracks
	TVector3 t1_long = t_sum * (t1.Dot(t_sum));
	TVector3 t1_perp = t1 - t1_long;
	TVector3 t2_long = t_sum * (t2.Dot(t_sum));
	TVector3 t2_perp = t2 - t2_long;
	TVector3 t3_long = t_sum * (t3.Dot(t_sum));
	TVector3 t3_perp = t3 - t3_long;

	t1_perp.SetMag(1);
	t2_perp.SetMag(1);
	t3_perp.SetMag(1);

	TVector3 perp_b_long1 = t1_perp * (b_perp1.Dot(t1_perp));
	//TVector3 perp_b_perp1 = b_perp1 - perp_b_long1;

	TVector3 perp_b_long2 = t2_perp * (b_perp2.Dot(t2_perp));
	//TVector3 perp_b_perp2 = b_perp2 - perp_b_long2;

	TVector3 perp_b_long3 = t3_perp * (b_perp3.Dot(t3_perp));
	//TVector3 perp_b_perp3 = b_perp3 - perp_b_long3;


	// the best point calculation
	TVector3 dV = t1 - t2;
	TVector3 dB = perp_b_long1 - perp_b_long2;
	double x12 = dV.Dot(dB) / dV.Mag2();
	dV = t2 - t3;
	dB = perp_b_long2 - perp_b_long3;
	double x23 = dV.Dot(dB) / dV.Mag2();
	dV = t3 - t1;
	dB = perp_b_long3 - perp_b_long1;
	double x31 = dV.Dot(dB) / dV.Mag2();

	// and systematic error of tracker
	double syst12 = tracker_error / t1.Angle(t2); // technically / Sin (or Tan), but Sin = Angle with these angles
	double syst23 = tracker_error / t2.Angle(t3); // technically / Sin (or Tan), but Sin = Angle with these angles
	double syst31 = tracker_error / t3.Angle(t1); // technically / Sin (or Tan), but Sin = Angle with these angles
	//double syst = pow(syst12, 2) + pow(syst23, 2) + pow(syst31, 2);
	double syst12_weight = 1/syst12;
	double syst23_weight = 1/syst23;
	double syst31_weight = 1/syst31;

	//double x_average = (x12 + x23 + x31) / 3;
	//double x_deviation = pow(x12 - x_average, 2) + pow(x23 - x_average, 2) + pow(x31 - x_average, 2);
	//double x_dev_syst = x_deviation + syst;

	// weighted average with tracker errors
	double x_average = (x12*syst12_weight + x23*syst23_weight + x31*syst31_weight) / (syst12_weight + syst23_weight + syst31_weight);
	double x_deviation = (syst12_weight*pow(x12 - x_average, 2) + syst23_weight*pow(x23 - x_average, 2) + syst31_weight*pow(x31 - x_average, 2))/(2*(syst12_weight + syst23_weight + syst31_weight)/3);

	return x_average / sqrt(x_deviation);
	}

// penalty deviation from direction
double b_simple_SV_pt_direction_penalty(
	int which, // for tests etc
	Float_t b1x, Float_t b1y, Float_t b1z,
	Float_t b2x, Float_t b2y, Float_t b2z,
	Float_t b3x, Float_t b3y, Float_t b3z,
	Float_t v1pt, Float_t v1eta, Float_t v1phi,
	Float_t v2pt, Float_t v2eta, Float_t v2phi,
	Float_t v3pt, Float_t v3eta, Float_t v3phi
	)
        {
	Float_t tracker_error = 0.2; // approximately systematic error on positions
	// 0.0002 is not noticable
	// 0.002 has slight effect
	// cannot pass it as function parameter in root/cint/tformula!!!!!!

	TVector3 b_vec1, b_vec2, b_vec3;
	// B DIRECTION NEGATIVE
	b_vec1.SetXYZ(-b1x, -b1y, -b1z);
	b_vec2.SetXYZ(-b2x, -b2y, -b2z);
	b_vec3.SetXYZ(-b3x, -b3y, -b3z);

	TVector3 t1, t2, t3;
	t1.SetPtEtaPhi(v1pt, v1eta, v1phi);
	t2.SetPtEtaPhi(v2pt, v2eta, v2phi);
	t3.SetPtEtaPhi(v3pt, v3eta, v3phi);

	TVector3 t_sum = t1 + t2 + t3;
	t_sum.SetMag(1);

	// after establishing direction of tau
	// tracks are only geometrical lines
	t1.SetMag(1);
	t2.SetMag(1);
	t3.SetMag(1);

	TVector3 b_long1 = t_sum * (b_vec1.Dot(t_sum));
	TVector3 b_perp1 = b_vec1 - b_long1;

	TVector3 b_long2 = t_sum * (b_vec2.Dot(t_sum));
	TVector3 b_perp2 = b_vec2 - b_long2;

	TVector3 b_long3 = t_sum * (b_vec3.Dot(t_sum));
	TVector3 b_perp3 = b_vec3 - b_long3;

	// in the perp plane find b long to tracks

	// perpendicular parts of tracks
	TVector3 t1_long = t_sum * (t1.Dot(t_sum));
	TVector3 t1_perp = t1 - t1_long;
	TVector3 t2_long = t_sum * (t2.Dot(t_sum));
	TVector3 t2_perp = t2 - t2_long;
	TVector3 t3_long = t_sum * (t3.Dot(t_sum));
	TVector3 t3_perp = t3 - t3_long;

	t1_perp.SetMag(1);
	t2_perp.SetMag(1);
	t3_perp.SetMag(1);

	TVector3 perp_b_long1 = t1_perp * (b_perp1.Dot(t1_perp));
	//TVector3 perp_b_perp1 = b_perp1 - perp_b_long1;

	TVector3 perp_b_long2 = t2_perp * (b_perp2.Dot(t2_perp));
	//TVector3 perp_b_perp2 = b_perp2 - perp_b_long2;

	TVector3 perp_b_long3 = t3_perp * (b_perp3.Dot(t3_perp));
	//TVector3 perp_b_perp3 = b_perp3 - perp_b_long3;


	// the best point calculation
	// approxiamtion, left for reference
	TVector3 dV = t1 - t2;
	TVector3 dB = perp_b_long1 - perp_b_long2;
	double x12 = dV.Dot(dB) / dV.Mag2();
	dV = t2 - t3;
	dB = perp_b_long2 - perp_b_long3;
	double x23 = dV.Dot(dB) / dV.Mag2();
	dV = t3 - t1;
	dB = perp_b_long3 - perp_b_long1;
	double x31 = dV.Dot(dB) / dV.Mag2();

	double pair_error1 = (perp_b_long1 + x12*t1 - (perp_b_long2 + x12*t2)).Mag();
	double pair_error2 = (perp_b_long2 + x23*t2 - (perp_b_long3 + x23*t3)).Mag();
	double pair_error3 = (perp_b_long3 + x31*t3 - (perp_b_long1 + x31*t1)).Mag();

	// and systematic error of tracker
	double syst12 = tracker_error / t1.Angle(t2); // technically / Sin (or Tan), but Sin = Angle with these angles
	double syst23 = tracker_error / t2.Angle(t3); // technically / Sin (or Tan), but Sin = Angle with these angles
	double syst31 = tracker_error / t3.Angle(t1); // technically / Sin (or Tan), but Sin = Angle with these angles
	//double syst12 = sqrt(pow(pair_error1, 2) + pow(tracker_error / t1.Angle(t2), 2)); // technically / Sin (or Tan), but Sin = Angle with these angles
	//double syst23 = sqrt(pow(pair_error2, 2) + pow(tracker_error / t2.Angle(t3), 2)); // technically / Sin (or Tan), but Sin = Angle with these angles
	//double syst31 = sqrt(pow(pair_error3, 2) + pow(tracker_error / t3.Angle(t1), 2)); // technically / Sin (or Tan), but Sin = Angle with these angles
	//double syst = pow(syst12, 2) + pow(syst23, 2) + pow(syst31, 2);
	double syst12_weight = 1/syst12;
	double syst23_weight = 1/syst23;
	double syst31_weight = 1/syst31;

	//double x_average = (x12 + x23 + x31) / 3;
	//double x_deviation = pow(x12 - x_average, 2) + pow(x23 - x_average, 2) + pow(x31 - x_average, 2);
	//double x_dev_syst = x_deviation + syst;

	// weighted average with tracker errors (old, left for reference)
	double x_average = (x12*syst12_weight + x23*syst23_weight + x31*syst31_weight) / (syst12_weight + syst23_weight + syst31_weight);
	double x_deviation = (syst12_weight*pow(x12 - x_average, 2) + syst23_weight*pow(x23 - x_average, 2) + syst31_weight*pow(x31 - x_average, 2))/(2*(syst12_weight + syst23_weight + syst31_weight)/3);

	// best point with full solution of linear equations -- it's not exactly possible! these are 3D vectors, they might not intersect!


	// use mean best points
	// and penalty convergence
	TVector3 bestpoint_12 = 0.5 * (perp_b_long1 + x12*t1 + perp_b_long2 + x12*t2);
	TVector3 bestpoint_23 = 0.5 * (perp_b_long2 + x23*t2 + perp_b_long3 + x23*t3);
	TVector3 bestpoint_31 = 0.5 * (perp_b_long3 + x31*t3 + perp_b_long1 + x31*t1);

	// average bp weighted with convergance
	// deviations are raw to penalty convergence
	syst12_weight = 1/pair_error1;
	syst23_weight = 1/pair_error2;
	syst31_weight = 1/pair_error3;
	TVector3 bestpoint_average_weighted = (bestpoint_12*syst12_weight + bestpoint_23*syst23_weight + bestpoint_31*syst31_weight) * (1/(syst12_weight + syst23_weight + syst31_weight));
	double bav_w = bestpoint_average_weighted.Mag();
	TVector3 bestpoint_dev12_w = bestpoint_12 - bestpoint_average_weighted;
	TVector3 bestpoint_dev23_w = bestpoint_23 - bestpoint_average_weighted;
	TVector3 bestpoint_dev31_w = bestpoint_31 - bestpoint_average_weighted;
	double bpdev_w = (pow(bestpoint_dev12_w.Mag(), 2) + pow(bestpoint_dev23_w.Mag(), 2) + pow(bestpoint_dev31_w.Mag(), 2))*(0.6666);

	// let's try penalty divertion
	bestpoint_average_weighted.SetMag(1);
	TVector3 bpdev_w1_long = bestpoint_average_weighted * (bestpoint_dev12_w.Dot(bestpoint_average_weighted)); // dev on projection
	TVector3 bpdev_w1_perp = bestpoint_dev12_w - bpdev_w1_long;
	TVector3 bpdev_w2_long = bestpoint_average_weighted * (bestpoint_dev23_w.Dot(bestpoint_average_weighted));
	TVector3 bpdev_w2_perp = bestpoint_dev23_w - bpdev_w2_long;
	TVector3 bpdev_w3_long = bestpoint_average_weighted * (bestpoint_dev31_w.Dot(bestpoint_average_weighted));
	TVector3 bpdev_w3_perp = bestpoint_dev31_w - bpdev_w3_long;

	double bpdev_w_long = (pow(bpdev_w3_long.Mag(), 2) + pow(bpdev_w2_long.Mag(), 2) + pow(bpdev_w1_long.Mag(), 2))*(0.6666);
	double bpdev_w_perp = (pow(bpdev_w3_perp.Mag(), 2) + pow(bpdev_w2_perp.Mag(), 2) + pow(bpdev_w1_perp.Mag(), 2))*(0.6666);


	//TVector3 bestpoint_12 = perp_b_long1 + x12*t1;
	//TVector3 bestpoint_23 = perp_b_long2 + x23*t2;
	//TVector3 bestpoint_31 = perp_b_long3 + x31*t3;
	TVector3 bestpoint_average = (bestpoint_12 + bestpoint_23 + bestpoint_31) * (0.3333);
	TVector3 bestpoint_dev12 = bestpoint_12 - bestpoint_average;
	TVector3 bestpoint_dev23 = bestpoint_23 - bestpoint_average;
	TVector3 bestpoint_dev31 = bestpoint_31 - bestpoint_average;

	// now project them on t_sum, find perp parts, find deviation along perp average
	// t_sum MIGHT BE WRONG
	TVector3 bpav_long = t_sum * (bestpoint_average.Dot(t_sum)); // 1) best point along the direction
	double bpav_proj = bpav_long.Mag();
	TVector3 bpav_perp = bestpoint_average - bpav_long;          // 2) best point perpendicular to direction, divertion from direction
	double bpav_divert = bpav_perp.Mag();
	// now deviations are needed
	bpav_perp.SetMag(1);

	TVector3 bpdev1_long = t_sum * (bestpoint_dev12.Dot(t_sum)); // dev on projection
	TVector3 bpdev1_perp = bpav_perp * (bpav_perp.Dot(bestpoint_dev12 - bpdev1_long)); // dev on divertion
	TVector3 bpdev2_long = t_sum * (bestpoint_dev23.Dot(t_sum));
	TVector3 bpdev2_perp = bpav_perp * (bpav_perp.Dot(bestpoint_dev23 - bpdev2_long));
	TVector3 bpdev3_long = t_sum * (bestpoint_dev31.Dot(t_sum));
	TVector3 bpdev3_perp = bpav_perp * (bpav_perp.Dot(bestpoint_dev31 - bpdev3_long));

	// dev-s of projection and divertion
	//double bpdev_proj = bp
	//double bpdev_proj = (syst12_weight*pow(bpdev1_long.Mag(), 2) + syst23_weight*pow(bpdev2_long.Mag(), 2) + syst31_weight*pow(bpdev3_long.Mag(), 2))/(2*(syst12_weight + syst23_weight + syst31_weight)/3);
	//double bpdev_divert = (syst12_weight*pow(bpdev1_perp.Mag(), 2) + syst23_weight*pow(bpdev2_perp.Mag(), 2) + syst31_weight*pow(bpdev3_perp.Mag(), 2))/(2*(syst12_weight + syst23_weight + syst31_weight)/3);
	// no weights
	double bpdev_proj   = (pow(bpdev1_long.Mag(), 2) + pow(bpdev2_long.Mag(), 2) + pow(bpdev3_long.Mag(), 2))*(0.6666);
	double bpdev_divert = (pow(bpdev1_perp.Mag(), 2) + pow(bpdev2_perp.Mag(), 2) + pow(bpdev3_perp.Mag(), 2))*(0.6666);

	switch (which)
		{
		// test finding paired best points
		case 1:
			// test best point 12
			return pair_error1;
		case 2:
			// test best point 23
			return pair_error2;
		case 3:
			// test best point 31
			return pair_error3;

		// test average best point etc
		case 4:
			return bpav_proj;
		case 5:
			return bpav_divert;
		case 6:
			return bpdev_proj;
		case 7:
			return bpdev_divert;

		// old calculation results, for reference
		case 8:
			return x_average;
		case 9:
			return x_deviation;
		case 10:
			if (x_deviation > 0)
				return x_average / sqrt(x_deviation);
			else
				return - x_average / sqrt(x_deviation);

		// different sign parameters
		case 11: // no penalties
			return (bpav_proj/bpdev_proj);
		case 12: // penalty divertion
			return (bpav_proj/bpdev_proj) / (1 + bpav_divert/bpdev_divert);
		case 13: // penalty convergence
			return (bpav_proj/bpdev_proj) / (1 + pair_error1 + pair_error2 + pair_error3);
		// deviations and penalty values
		case 14:
			return bpav_divert/bpdev_divert;
		case 15:
			return pair_error1 + pair_error2 + pair_error3;

		// convergance-weighted average, but raw dev to penalty convergence
		case 16:
			return bav_w;
		case 17:
			return bpdev_w;
		case 18:
			return bav_w / bpdev_w;
		case 19:
			return bpdev_w_long;
		case 20:
			return bpdev_w_perp;
		case 21:
			return (bav_w / bpdev_w_long) / (1 + bpdev_w_perp);
		}

	//return x_average / sqrt(x_deviation);
	return (bpav_proj/bpdev_proj) / (1 + bpav_divert/bpdev_divert + pair_error1 + pair_error2 + pair_error3);
	}












// SV in xy plane, ignoring z
//

// shouldn't be rperp
double xy_b_track_angle(
	Float_t b1x, Float_t b1y, Float_t b1z,
	Float_t v1eta, Float_t v1phi
	)
        {
	TVector3 b1, tau;

	b1.SetXYZ(b1x, b1y, 0);
	tau.SetPtEtaPhi(1, v1eta, v1phi);
	tau.SetZ(0);

	return b1.Angle(tau);
	}


// degenerative ROOT and its' crap packages...

// x
// y
//  1 2
struct vec2 {
	double x;
	double y;
	};

struct mat2 {
	double x1;
	double x2;
	double y1;
	double y2;
	};

struct vec2 dot2(struct mat2 mat, struct vec2 vec)
	//struct mat2 tracks_inv_m = {.x1 = tracks[0][0], .x2 = tracks[1][0], .y1 = tracks[0][1], .y2 = tracks[1][1]};
	{
	return {vec.x*mat.x1 + vec.y*mat.x2, vec.x*mat.y1 + vec.y*mat.y2};
	}

//struct vec2 xy_b_track_pair_bestpoint(
double xy_b_delirium_Phi(
	Float_t b1x, Float_t b1y
	)
        {
	TVector3 b1;
	b1.SetXYZ(b1x, b1y, 0);
	return b1.Phi();
	}

double xy_b_delirium_tauPhi(
	Float_t taueta, Float_t tauphi
	)
        {
	TVector3 b1, v2;
	b1.SetPtEtaPhi(1, taueta, tauphi);
	b1.SetZ(0);
	b1.SetMag(1);
	v2.SetPtEtaPhi(1, 0, tauphi);
	return b1.Angle(v2);
	}


//struct vec2 xy_b_track_pair_bestpoint(
//	int which,
//	Float_t b1x, Float_t b1y, Float_t b1z,
//	Float_t b2x, Float_t b2y, Float_t b2z,
//	Float_t v1eta, Float_t v1phi,
//	Float_t v2eta, Float_t v2phi
//	)
//        {
//	struct vec2 out = {0,0};
//	return out;
//	}

struct vec2 xy_b_track_pair_bestpoint(
//double xy_b_track_pair_bestpoint_test(
	int which,
	Float_t b1x, Float_t b1y, Float_t b1z,
	Float_t b2x, Float_t b2y, Float_t b2z,
	Float_t v1eta, Float_t v1phi,
	Float_t v2eta, Float_t v2phi
	//Float_t tauphi
	)
        {
	TVector3 b1, b2, t1, t2;

	// b is inverted in ntupler
	b1.SetXYZ(-b1x, -b1y, 0);
	b2.SetXYZ(-b2x, -b2y, 0);

	t1.SetPtEtaPhi(1, v1eta, v1phi);
	t1.SetZ(0);
	t2.SetPtEtaPhi(1, v2eta, v2phi);
	t2.SetZ(0);

	t1.SetMag(1);
	t2.SetMag(1);

	TVector3 t_sum = t1 + t2;
	t_sum.SetMag(1);

	//TVector3 tau;
	//tau.SetPtEtaPhi(1, 0, tauphi);
	//tau.SetMag(1);

	//TVectorD* t1_xy = new TVectorD(2);
	//TVectorD* t2_xy = new TVectorD(2);
	//t1_xy->SetXY(t1.X(), t1.Y());
	//t2_xy->SetXY(t2.X(), t2.Y());
	//ROOT::Math::SVector2 t1_xy(t1.X(), t1.Y());
	//ROOT::Math::SVector2 t2_xy(t2.X(), t2.Y());
	//SVector2 t1_xy(t1.X(), t1.Y());
	//SVector2 t2_xy(t2.X(), t2.Y());
	//SVector<double,2> t1_xy(t1.X(), t1.Y());
	//SVector<double,2> t2_xy(t2.X(), t2.Y());
	struct vec2 t1_xy = {.x = t1.X(), .y =  t1.Y()};
	struct vec2 t2_xy = {.x = t2.X(), .y =  t2.Y()};

	TMatrixD tracks_orig(2,2);
	tracks_orig[0][0] =  t1.X();
	tracks_orig[0][1] =  t1.Y();
	tracks_orig[1][0] = -t2.X();
	tracks_orig[1][1] = -t2.Y();
	//SMatrix22 tracks;
	TMatrixD tracks(2,2);
	tracks[0][0] =  t1.X();
	tracks[0][1] =  t1.Y();
	tracks[1][0] = -t2.X();
	tracks[1][1] = -t2.Y();
	// not sure about the row/coloumn order here
	// t1 multiplies by x
	// -t2 by y
	tracks.Invert();
	// the only thing needed from this crap

	//TVectorD *bs_xy = new TVectorD(2);
	//bs_xy->SetXY(b2.x() - b1.x(), b2.y() - b1.y());
	//ROOT::Math::SVector2 bs_xy(b2.x() - b1.x(), b2.y() - b1.y());
	//SVector2 bs_xy(b2.x() - b1.x(), b2.y() - b1.y());
	//auto solution = tracks * (bs_xy);

	struct vec2 bs_xy = {.x = b2.x() - b1.x(), .y = b2.y() - b1.y()};
	//struct mat2 tracks_inv_m = {.x1 = tracks[0][0], .x2 = tracks[0][1], .y1 = tracks[1][0], .y2 = tracks[1][1]};
	struct mat2 tracks_inv_m = {.x1 = tracks[0][0], .x2 = tracks[1][0], .y1 = tracks[0][1], .y2 = tracks[1][1]};

	struct vec2 solution = dot2(tracks_inv_m, bs_xy);

	//auto best_point = t1_xy*solution.X() + b1;
	// should = to
	//auto best_point2 = t2_xy*solution.y() + b2;

	struct vec2 best_point1, best_point2;

	// tests
	// t1 - x parameter solution
	// t2 - y parameter solution
	best_point1.x = t1_xy.x * solution.x + b1.x();
	best_point1.y = t1_xy.y * solution.x + b1.y();
	best_point2.x = t2_xy.x * solution.y + b2.x();
	best_point2.y = t2_xy.y * solution.y + b2.y();

	TVector3 t1_long, t1_perp, t2_long, t2_perp;

	/*
	switch (which)
		{
		case -12:
			return b1.X();
		case -13:
			return b1.Y();
		case -14:
			return b1.Z();
		// check solution
		case 1:
			//return (solution.x * t1[0] - solution.y * t2[0]) - (b2.x() - b1.x());
			return (solution.x * t1.x() - solution.y * t2.x()) - (b2.x() - b1.x());
		case 2:
			//return (solution.x * t1[1] - solution.y * t2[1]) - (b2.y() - b1.y());
			return (solution.x * t1.y() - solution.y * t2.y()) - (b2.y() - b1.y());
		case 3:
			return (tracks_orig * tracks)[0][0];
			break;
		case 4:
			return (tracks_orig * tracks)[0][1];
			break;
		case 5:
			return (tracks_orig * tracks)[1][0];
			break;
		case 6:
			return (tracks_orig * tracks)[1][1];
			break;
		case 7:
			return abs(best_point1.x - best_point2.x) + abs(best_point1.y - best_point2.y);
			break;
		case 8:
			return solution.x;
		case 9:
			return solution.y;
		case 10:
			// track perp part is along b or against b
			//t1_long = t_sum * (t1.Dot(t_sum));
			//t1_perp = t1 - t1_long;
			t1_long = tau * (t1.Dot(tau));
			t1_perp = t1 - t1_long;
			return t1_perp * b1;

		case 11:
			//t2_long = t_sum * (t2.Dot(t_sum));
			//t2_perp = t2 - t2_long;
			t2_long = tau * (t2.Dot(tau));
			t2_perp = t2 - t2_long;
			return t2_perp * b2;

		case 12:
			// project t1 on t2
			t2.SetMag(1);
			t1_long = t2 * (t1.Dot(t2));
			t1_perp = t1 - t1_long;
			return t1_perp * b1;
		case 13:
			// project t1 on t2
			t1.SetMag(1);
			t2_long = t1 * (t2.Dot(t1));
			t2_perp = t2 - t2_long;
			return t2_perp * b2;

		}
	return 0;
	}
	*/

	struct vec2 best_point;
	if (which > 0)
		{
		best_point.x = t1_xy.x * solution.x + b1.x();
		best_point.y = t1_xy.y * solution.x + b1.y();
		}
	else
		{
		best_point.x = t2_xy.x * solution.y + b2.x();
		best_point.y = t2_xy.y * solution.y + b2.y();
		}
	//struct vec2 best_point = (which > 0 ? struct vec2 {.x = t1_xy.x * solution.x + b1.x(), .y= t1_xy.y * solution.y + b1.y()} : struct vec2 {.x = t2_xy.x * solution.x + b2.x(), .y= t2_xy.y * solution.y + b2.y()});
	//struct vec2 best_point2 = {.x = t2_xy.x * solution.x + b2.x(), .y= t2_xy.y * solution.y + b2.y()};
	return best_point;
	}


double xy_b_track_pair_bestpoint_tests(
	int which,
	Float_t b1x, Float_t b1y, Float_t b1z,
	Float_t b2x, Float_t b2y, Float_t b2z,
	Float_t v1eta, Float_t v1phi,
	Float_t v2eta, Float_t v2phi
	)
	{
	TVector3 b1, b2, t1, t2;

	b1.SetXYZ(b1x, b1y, 0);
	t1.SetPtEtaPhi(1, v1eta, v1phi);
	t1.SetZ(0);

	b2.SetXYZ(b2x, b2y, 0);
	t2.SetPtEtaPhi(1, v2eta, v2phi);
	t2.SetZ(0);

	t1.SetMag(1);
	t2.SetMag(1);

	//TVectorD* t1_xy = new TVectorD(2);
	//TVectorD* t2_xy = new TVectorD(2);
	//t1_xy->SetXY(t1.X(), t1.Y());
	//t2_xy->SetXY(t2.X(), t2.Y());
	//ROOT::Math::SVector2 t1_xy(t1.X(), t1.Y());
	//ROOT::Math::SVector2 t2_xy(t2.X(), t2.Y());
	//SVector2 t1_xy(t1.X(), t1.Y());
	//SVector2 t2_xy(t2.X(), t2.Y());
	//SVector<double,2> t1_xy(t1.X(), t1.Y());
	//SVector<double,2> t2_xy(t2.X(), t2.Y());
	struct vec2 t1_xy = {.x = t1.X(), .y =  t1.Y()};
	struct vec2 t2_xy = {.x = t2.X(), .y =  t2.Y()};

	//SMatrix22 tracks;
	TMatrixD tracks_orig(2,2);
	tracks_orig[0][0] = t1[0];
	tracks_orig[0][1] = t1[1];
	tracks_orig[1][0] = -t2[0];
	tracks_orig[1][1] = -t2[1];
	TMatrixD tracks(2,2);
	tracks[0][0] = t1[0];
	tracks[0][1] = t1[1];
	tracks[1][0] = -t2[0];
	tracks[1][1] = -t2[1];
	// not sure about the row/coloumn order here
	// t1 multiplies by x
	// -t2 by y
	tracks.Invert();
	// the only thing needed from this crap

	//TVectorD *bs_xy = new TVectorD(2);
	//bs_xy->SetXY(b2.x() - b1.x(), b2.y() - b1.y());
	//ROOT::Math::SVector2 bs_xy(b2.x() - b1.x(), b2.y() - b1.y());
	//SVector2 bs_xy(b2.x() - b1.x(), b2.y() - b1.y());
	//auto solution = tracks * (bs_xy);

	struct vec2 bs_xy = {.x = b2.x() - b1.x(), .y = b2.y() - b1.y()};
	//struct mat2 tracks_inv_m = {.x1 = tracks[0][0], .x2 = tracks[0][1], .y1 = tracks[1][0], .y2 = tracks[1][1]};
	struct mat2 tracks_inv_m = {.x1 = tracks[0][0], .x2 = tracks[1][0], .y1 = tracks[0][1], .y2 = tracks[1][1]};

	struct vec2 solution = dot2(tracks_inv_m, bs_xy);

	//auto best_point = t1_xy*solution.X() + b1;
	// should = to
	//auto best_point2 = t2_xy*solution.y() + b2;

	//struct vec2 best_point = (which > 0 ? {.x = t1_xy.x * solution.x + b1.x(), .y= t1_xy.y * solution.y + b1.y()} : {.x = t2_xy.x * solution.x + b2.x(), .y= t2_xy.y * solution.y + b2.y()});
	//struct vec2 best_point2 = {.x = t2_xy.x * solution.x + b2.x(), .y= t2_xy.y * solution.y + b2.y()};
	//return best_point;

	switch (which)
		{
		// check solution
		case 1:
			return (solution.x * t1[0] - solution.y * t2[0]) - (b2.x() - b1.x());
			break;
		case 2:
			return (solution.x * t1[1] - solution.y * t2[1]) - (b2.y() - b1.y());
			break;
		case 3:
			return (tracks_orig * tracks)[0][0];
			break;
		case 4:
			return (tracks_orig * tracks)[0][1];
			break;
		case 5:
			return (tracks_orig * tracks)[1][0];
			break;
		case 6:
			return (tracks_orig * tracks)[1][1];
			break;
		}
	return 0;
	}

// 
double xy_plane_b_tau_directions(
	int which,
	Float_t b1x, Float_t b1y, Float_t b1z,
	Float_t b2x, Float_t b2y, Float_t b2z,
	Float_t b3x, Float_t b3y, Float_t b3z,
	Float_t v1eta, Float_t v1phi,
	Float_t v2eta, Float_t v2phi,
	Float_t v3eta, Float_t v3phi
	)
        {
	TVector3 b1, b2, b3;
	b1.SetXYZ(b1x, b1y, 0);
	b2.SetXYZ(b2x, b2y, 0);
	b3.SetXYZ(b3x, b3y, 0);
	TVector3 t1, t2, t3;
	t1.SetPtEtaPhi(1, 0, v1phi); // maybe just eta = 0?
	t2.SetPtEtaPhi(1, 0, v2phi);
	t3.SetPtEtaPhi(1, 0, v3phi);
	TVector3 t_sum = t1 + t2 + t3;
	t_sum.SetMag(1);

	switch(which)
		{
		case 1:
			return t_sum.Dot(b1);
		case 2:
			return t_sum.Dot(b2);
		case 3:
			return t_sum.Dot(b3);
		}

	return 0;
	}

// 
//TVectorD xy_b_track_pair_bestpoint(
//double xy_b_track_pair_bestpoint(
double xy_average_SV(
	int which,
	Float_t b1x, Float_t b1y, Float_t b1z,
	Float_t b2x, Float_t b2y, Float_t b2z,
	Float_t b3x, Float_t b3y, Float_t b3z,
	Float_t v1eta, Float_t v1phi,
	Float_t v2eta, Float_t v2phi,
	Float_t v3eta, Float_t v3phi
	)
        {
	struct vec2 bp12 = xy_b_track_pair_bestpoint( 1,
		b1x, b1y, b1z,
		b2x, b2y, b2z,
		v1eta, v1phi,
		v2eta, v2phi);
	struct vec2 bp23 = xy_b_track_pair_bestpoint( 1,
		b2x, b2y, b2z,
		b3x, b3y, b3z,
		v2eta, v2phi,
		v3eta, v3phi);
	struct vec2 bp31 = xy_b_track_pair_bestpoint( 1,
		b3x, b3y, b3z,
		b1x, b1y, b1z,
		v3eta, v3phi,
		v1eta, v1phi);

	double av_x = (bp12.x + bp23.x + bp31.x) / 3;
	double av_y = (bp12.y + bp23.y + bp31.y) / 3;

	double dev_x = sqrt((pow(bp12.x - av_x, 2) + pow(bp23.x - av_x, 2) + pow(bp31.x - av_x, 2)) / 3);
	double dev_y = sqrt((pow(bp12.y - av_y, 2) + pow(bp23.y - av_y, 2) + pow(bp31.y - av_y, 2)) / 3);

	TVector3 t1, t2, t3, t_sum;
	double tau_dir = 0;
	switch (which)
		{
		case -3:
			// proj on tau scaled by deviations
			t1.SetPtEtaPhi(1, 0, v1phi); // maybe just eta = 0?
			t2.SetPtEtaPhi(1, 0, v2phi);
			t3.SetPtEtaPhi(1, 0, v3phi);
			//t1.SetZ(0);
			//t2.SetZ(0);
			//t3.SetZ(0);
			t_sum = t1 + t2 + t3;
			t_sum.SetMag(1);

			tau_dir = t_sum.X() * av_x / dev_x + t_sum.Y() * av_y / dev_y;

			return tau_dir;

		case -2:
			// flightLength in xy proj on tau
			t1.SetPtEtaPhi(1, 0, v1phi); // maybe just eta = 0?
			t2.SetPtEtaPhi(1, 0, v2phi);
			t3.SetPtEtaPhi(1, 0, v3phi);
			//t1.SetZ(0);
			//t2.SetZ(0);
			//t3.SetZ(0);
			t_sum = t1 + t2 + t3;
			t_sum.SetMag(1);

			tau_dir = t_sum.X() * av_x + t_sum.Y() * av_y;

			return tau_dir;

		case -1:
			t1.SetPtEtaPhi(1, 0, v1phi); // maybe just eta = 0?
			t2.SetPtEtaPhi(1, 0, v2phi);
			t3.SetPtEtaPhi(1, 0, v3phi);
			//t1.SetZ(0);
			//t2.SetZ(0);
			//t3.SetZ(0);
			t_sum = t1 + t2 + t3;
			t_sum.SetMag(1);

			tau_dir = t_sum.X() * av_x + t_sum.Y() * av_y;

			if (tau_dir > 0)
				return  sqrt(pow(av_x, 2) + pow(av_y, 2));
			else
				return -sqrt(pow(av_x, 2) + pow(av_y, 2));

		case 0:
			t1.SetPtEtaPhi(1, 0, v1phi); // maybe just eta = 0?
			t2.SetPtEtaPhi(1, 0, v2phi);
			t3.SetPtEtaPhi(1, 0, v3phi);
			//t1.SetZ(0);
			//t2.SetZ(0);
			//t3.SetZ(0);
			t_sum = t1 + t2 + t3;
			t_sum.SetMag(1);

			tau_dir = t_sum.X() * av_x + t_sum.Y() * av_y;

			if (tau_dir > 0)
				return sqrt(pow(av_x/dev_x, 2) + pow(av_y/dev_y, 2));
			else
				return -sqrt(pow(av_x/dev_x, 2) + pow(av_y/dev_y, 2));
			//return av_x;
		case 1:
			return av_x;
		case 2:
			return av_y;
		case 3:
			return dev_x;
		case 4:
			return dev_y;
		}

	return 0;
	}





