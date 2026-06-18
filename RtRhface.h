#pragma once

#include "RtCommon.h"
#include "RtTriface.h"

class RtRhface
{
public:

	RtRhface(void);

	RtRhface(unsigned short id);

	void SetID(unsigned short id);

	void SetVertsID(int id1, int id2, int id3, int id4);

	void SetVertsgeo(const GeoCoord& v0, const GeoCoord& v1, const GeoCoord& v2, const GeoCoord& v3);

	void SetVerts3d(const DVec3D& v0, const DVec3D& v1, const DVec3D& v2, const DVec3D& v3);

	void SetNeighFaceID(int id1, int id2, int id3, int id4);

	int* const GetVertsId();

	int* const GetNeighFaceId();

	GeoCoord* const GetVertsGeo();

	unsigned short const GetID();

	~RtRhface(void);

private:

	unsigned short id_;

	int verts_id[4];

	GeoCoord verts_geo[4];

	DVec3D verts_3D[4];

	int NeighFace_id[4];
};
