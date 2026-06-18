#pragma once
#include "RtCommon.h"

class RtTriface
{
public:

	RtTriface(void);

	RtTriface(unsigned short id);

	void SetId(unsigned short id);

	void SetVertsId(int id1, int id2, int id3);

	void SetVertsgeo(const GeoCoord& v0, const GeoCoord& v1, const GeoCoord& v2);

	void SetVerts3d(const DVec3D& v0, const DVec3D& v1, const DVec3D& v2);

	int const GetID();

	int* const GetVertsID();

	GeoCoord* const GetVertsgeo();

	~RtTriface(void);

protected:

	unsigned short id_;

	int verts_id[3];

	GeoCoord verts_geo[3];

	DVec3D verts_3D[3];
};
