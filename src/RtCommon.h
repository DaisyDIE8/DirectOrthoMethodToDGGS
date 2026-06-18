#pragma once

#include "stdafx.h"
#include "RtConstants.h"

enum InterpolateType
{
	BILINEAR,
	NEAREST,
};

struct DVec2D
{
	double x;

	double y;

	DVec2D() : x(0), y(0)
	{
	}

	DVec2D(double px, double py) : x(px), y(py)
	{
	}
};

struct IVec2D
{
	int i;

	int j;

	IVec2D() : i(0), j(0)
	{
	}

	IVec2D(int px, int py) : i(px), j(py)
	{
	}
};

struct DVec3D
{
	double x;

	double y;

	double z;

	DVec3D() : x(0.0), y(0.0), z(0.0)
	{
	}

	DVec3D(double px, double py, double pz) : x(px), y(py), z(pz)
	{
	}

	DVec3D crossproduct(const DVec3D& pt1, const DVec3D& pt2)
	{
		DVec3D m;
		m.x = pt1.y * pt2.z - pt2.y * pt1.z;
		m.y = -pt1.x * pt2.z + pt2.x * pt1.z;
		m.z = pt1.x * pt2.y - pt2.x * pt1.y;

		return m;
	}

	double magnitude(void) const
	{
		return sqrt(x * x + y * y + z * z);
	}

	double dotproduct(const DVec3D& pt1)
	{
		return this->x * pt1.x + this->y * pt1.y + this->z * pt1.z;
	}
};

struct GeoCoord
{
	double lat;

	double lon;

	GeoCoord() : lat(0.0), lon(0.0)
	{
	}

	GeoCoord(double plat, double plon) : lat(plat), lon(plon)
	{
	}
};

struct DGridCoord2D
{
	unsigned short idx;

	DVec2D pt;

	DGridCoord2D() : idx(0), pt(DVec2D())
	{
	}

	DGridCoord2D(unsigned short i, const DVec2D& p) : idx(i), pt(p)
	{
	}

	DGridCoord2D(unsigned short i, double x, double y) : idx(i), pt(x, y)
	{
	}
};

struct DGridCoord2Dij
{
	unsigned short idx;

	IVec2D pt;

	DGridCoord2Dij() : idx(0), pt(IVec2D())
	{
	}

	DGridCoord2Dij(unsigned short i, const IVec2D& p) : idx(i), pt(p)
	{
	}

	DGridCoord2Dij(unsigned short i, int x, int y) : idx(i), pt(x, y)
	{
	}
};

struct PreCompGeo
{
	GeoCoord pt;

	double sinLat;

	double sinLon;

	double cosLat;

	double cosLon;
};

struct PreCompInTri
{
	double p0x0;
	double p0y0;
	double p0z0;

	double p0x1;
	double p0y1;
	double p0z1;

	double p0x2;
	double p0y2;
	double p0z2;

	double t00;
	double t01;
	double t02;
};

double GreatCircleDistance(DVec3D& s, DVec3D& e);

DVec3D llxyz(const GeoCoord& g);

GeoCoord xyzll(const DVec3D& v);

double Intersec(const GeoCoord& pt, const GeoCoord& pt1, const GeoCoord& pt2);
