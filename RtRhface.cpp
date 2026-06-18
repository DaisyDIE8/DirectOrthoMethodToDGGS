#include "stdafx.h"
#include "RtRhface.h"

RtRhface::RtRhface(void)
{
	id_ = -1;

	for (int i = 0; i < 4; i++)
	{
		verts_id[i] = -1;
		verts_geo[i] = GeoCoord();
		verts_3D[i] = DVec3D();
	}
}

RtRhface::RtRhface(unsigned short id)
{
	id_ = id;

	for (int i = 0; i < 4; i++)
	{
		verts_id[i] = -1;
		verts_geo[i] = GeoCoord();
		verts_3D[i] = DVec3D();
	}
}

RtRhface::~RtRhface(void)
{
}

void RtRhface::SetID(unsigned short id)
{
	id_ = id;
}

void RtRhface::SetVertsID(int id1, int id2, int id3, int id4)
{
	verts_id[0] = id1;
	verts_id[1] = id2;
	verts_id[2] = id3;
	verts_id[3] = id4;
}

void RtRhface::SetVertsgeo(const GeoCoord& v0, const GeoCoord& v1, const GeoCoord& v2, const GeoCoord& v3)
{
	verts_geo[0] = v0;
	verts_geo[1] = v1;
	verts_geo[2] = v2;
	verts_geo[3] = v3;
}

void RtRhface::SetVerts3d(const DVec3D& v0, const DVec3D& v1, const DVec3D& v2, const DVec3D& v3)
{
	verts_3D[0] = v0;
	verts_3D[1] = v1;
	verts_3D[2] = v2;
	verts_3D[3] = v3;
}

int* const RtRhface::GetVertsId()
{
	return verts_id;
}

unsigned short const RtRhface::GetID()
{
	return id_;
}

GeoCoord* const RtRhface::GetVertsGeo()
{
	return verts_geo;
}

void RtRhface::SetNeighFaceID(int id1, int id2, int id3, int id4)
{
	NeighFace_id[0] = id1;
	NeighFace_id[1] = id2;
	NeighFace_id[2] = id3;
	NeighFace_id[3] = id4;
}

int* const RtRhface::GetNeighFaceId()
{
	return NeighFace_id;
}