#include "stdafx.h"
#include "RtCommon.h"

double GreatCircleDistance(DVec3D& s, DVec3D& e)
{
	double ls2 = s.x * s.x + s.y * s.y + s.z * s.z;
	double le2 = e.x * e.x + e.y * e.y + e.z * e.z;
	double ls = sqrt(ls2);
	double le = sqrt(le2);

	if (fabs(ls - le) > TOLERANCE)
	{
		cout << "S and E are not in the same sphere!" << endl;
		return 0.0;
	}

	double c2 = (s.x - e.x) * (s.x - e.x) + (s.y - e.y) * (s.y - e.y) + (s.z - e.z) * (s.z - e.z);

	double d = acos((ls2 + le2 - c2) / (2.0 * ls * le));

	return d;
}

DVec3D llxyz(const GeoCoord& sv)
{
	DVec3D v;
	v.x = cos(sv.lat) * cos(sv.lon);
	v.y = cos(sv.lat) * sin(sv.lon);
	v.z = sin(sv.lat);

	const double epsilon = 0.000000000000001;
	if (fabs(v.x) < epsilon)
		v.x = 0.0;
	if (fabs(v.y) < epsilon)
		v.y = 0.0;
	if (fabs(v.z) < epsilon)
		v.z = 0.0;

	return v;
}

GeoCoord xyzll(const DVec3D& v0)
{
	GeoCoord sv;

	sv.lat = -1.0;
	sv.lon = -1.0;

	DVec3D v = v0;

	if (fabs(v.z) - 1 < TOLERANCE)
	{
		if (v.z > 1)
			v.z = 1;
		if (v.z < -1)
			v.z = -1;

		sv.lat = asin(v.z);

		if ((sv.lat == M_PI / 2) || (sv.lat == -M_PI / 2))
			sv.lon = 0;
		else
			sv.lon = atan2(v.y, v.x);

		return sv;
	}
	else
	{
		printf("Error: in function xyzll, asin domain error.\n");
		return sv;
	}
}

double Intersec(const GeoCoord& pt, const GeoCoord& pt1, const GeoCoord& pt2)
{
	DVec3D a_ = llxyz(pt);
	DVec3D b_ = llxyz(pt1);
	DVec3D c_ = llxyz(pt2);

	DVec3D t2 = a_.crossproduct(a_, c_);

	if (fabs(t2.magnitude()) < TOLERANCE)
	{
		return 0.0;
	}

	DVec3D t1 = a_.crossproduct(a_, b_);

	double tt = t1.dotproduct(t2) / (t1.magnitude() * t2.magnitude());

	if ((tt > 1) && (tt - 1) < TOLERANCE)
		tt = 1;

	double angle = acos(tt);

	return angle;
}