#include "TVector3.h"
#include "TMatrixD.h"


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
double b_cor_angle_trk(
	Float_t b1x, Float_t b1y, Float_t b1z,
	Float_t v1eta, Float_t v1phi,
	Float_t v2eta, Float_t v2phi
	)
        {
	TVector3 b_vec, trk, tau;
	b_vec.SetXYZ(b1x, b1y, b1z);
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


// ok, it's late, let's try simple SV with this plane direction
double b_track_plane_intersections_raw_SV(
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
	// got plane direction

	// find perpendicular b-s

	TVector3 b_long1 = average * (b_vec1.Dot(average));
	TVector3 b_perp1 = b_vec1 - b_long1;
	TVector3 b_long2 = average * (b_vec2.Dot(average));
	TVector3 b_perp2 = b_vec2 - b_long2;
	TVector3 b_long3 = average * (b_vec3.Dot(average));
	TVector3 b_perp3 = b_vec3 - b_long3;

	// in the perp plane find b long to tracks

	// perpendicular parts of tracks
	TVector3 t1_long = average * (t1.Dot(average));
	TVector3 t1_perp = t1 - t1_long;
	TVector3 t2_long = average * (t2.Dot(average));
	TVector3 t2_perp = t2 - t2_long;
	TVector3 t3_long = average * (t3.Dot(average));
	TVector3 t3_perp = t3 - t3_long;

	// no additional correction to b-s

	// the best point calculation
	TVector3 dV = t1 - t2;
	TVector3 dB = b_perp1 - b_perp2;
	double x12 = dV.Dot(dB) / dV.Mag2();
	dV = t2 - t3;
	dB = b_perp2 - b_perp3;
	double x23 = dV.Dot(dB) / dV.Mag2();
	dV = t3 - t1;
	dB = b_perp3 - b_perp1;
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










