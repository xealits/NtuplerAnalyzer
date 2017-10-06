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

	double x_average = (x12 + x23 + x31) / 3;
	double x_deviation = sqrt(pow(x12 - x_average, 2) + pow(x23 - x_average, 2) + pow(x31 - x_average, 2));

	return x_average / x_deviation;
	}










