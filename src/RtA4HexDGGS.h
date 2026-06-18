#pragma once
#include "rtriacontahedron.h"

class RtSphHexDGGS :
	public Rtriacontahedron
{
public:
	RtSphHexDGGS(double, unsigned short);

	int GetSphTri(const GeoCoord& pt, const PreCompInTri* intri);

	DGridCoord2D RTForwISlice(const GeoCoord& p);

	GeoCoord RTInvSlice(DGridCoord2D& p);

	DGridCoord2D GetTricoordfromIJ(const DGridCoord2Dij);

	DGridCoord2Dij GetIJfromTricoord(const DGridCoord2D& gp);

	GeoCoord CalGeoFromij(const DGridCoord2Dij& ij);

	DGridCoord2Dij CalijFromGeo(const GeoCoord& geo);

	void SetL();

public:

	unsigned short level_;

	int m_;

	double lengthIJ_;

	double lengthK_;

	double lengthH_;

private:

	double R_;

	GeoCoord geoVerts_[32];

	DVec3D desVerts_[32];

	PreCompInTri ptIn_[60];
};
