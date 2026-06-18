#include "stdafx.h"
#include "RtTriface.h"

RtTriface::RtTriface(void)
{
	id_ = -1;

	verts_id[0] = -1;
	verts_id[1] = -1;
	verts_id[2] = -1;

	verts_geo[0] = GeoCoord();
	verts_geo[1] = GeoCoord();
	verts_geo[2] = GeoCoord();

	verts_3D[0] = DVec3D();
	verts_3D[1] = DVec3D();
	verts_3D[2] = DVec3D();
}

RtTriface::RtTriface(unsigned short id)
{
	id_ = id;

	verts_id[0] = -1;
	verts_id[1] = -1;
	verts_id[2] = -1;

	verts_geo[0] = GeoCoord();
	verts_geo[1] = GeoCoord();
	verts_geo[2] = GeoCoord();

	verts_3D[0] = DVec3D();
	verts_3D[1] = DVec3D();
	verts_3D[2] = DVec3D();
}

RtTriface::~RtTriface(void)
{
}

void RtTriface::SetId(unsigned short id)
{
	id_ = id;
}

void RtTriface::SetVertsId(int id1, int id2, int id3)
{
	verts_id[0] = id1;
	verts_id[1] = id2;
	verts_id[2] = id3;
}

void RtTriface::SetVertsgeo(const GeoCoord& v0, const GeoCoord& v1, const GeoCoord& v2)
{
	verts_geo[0] = v0;
	verts_geo[1] = v1;
	verts_geo[2] = v2;
}

void RtTriface::SetVerts3d(const DVec3D& v0, const DVec3D& v1, const DVec3D& v2)
{
	verts_3D[0] = v0;
	verts_3D[1] = v1;
	verts_3D[2] = v2;
}

int const RtTriface::GetID()
{
	return id_;
}

int* const RtTriface::GetVertsID()
{
	return verts_id;
}

GeoCoord* const RtTriface::GetVertsgeo()
{
	return verts_geo;
}