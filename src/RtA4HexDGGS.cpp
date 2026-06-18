#include "stdafx.h"
#include "RtA4HexDGGS.h"

RtSphHexDGGS::RtSphHexDGGS(double R, unsigned short lev)
	:level_(lev), R_(R), Rtriacontahedron(R* sqrt(M_PI / (3 * sqrt(5))))
{
	m_ = (int)pow(2, level_);
	lengthIJ_ = L_ / m_;
	lengthK_ = b_ / m_;
	lengthH_ = a_ / m_;

	geoVerts_[0].lat = M_PI / 2.0;
	geoVerts_[0].lon = 0;

	for (int i = 1; i <= 5; i++)
	{
		geoVerts_[i].lat = 0.91843818700958013;
		geoVerts_[i].lon = (72 * i - 36) * M_PI / 180.0;

		geoVerts_[i + 5].lat = 26.565051177 * M_PI / 180.0;
		geoVerts_[i + 5].lon = (72 * i - 72) * M_PI / 180.0;

		geoVerts_[i + 10].lat = 0.18871053078293504;
		geoVerts_[i + 10].lon = (72 * i - 36) * M_PI / 180.0;

		geoVerts_[i + 15].lat = -0.18871053078293504;
		geoVerts_[i + 15].lon = (72 * i) * M_PI / 180.0;

		geoVerts_[i + 20].lat = -26.565051177 * M_PI / 180.0;
		geoVerts_[i + 20].lon = (72 * i - 36) * M_PI / 180.0;

		geoVerts_[i + 25].lat = -0.91843818700958013;
		geoVerts_[i + 25].lon = (72 * i) * M_PI / 180.0;
	}

	geoVerts_[31].lat = -M_PI / 2.0;
	geoVerts_[31].lon = 0.0;

	SetL();

	for (int i = 0; i < 60; i++)
	{
		GeoCoord* t = tri_[i].GetVertsgeo();

		DVec3D ptri[3];
		for (int j = 0; j < 3; j++)
			ptri[j] = llxyz(t[j]);

		ptIn_[i].p0x0 = ptri[0].y * ptri[1].z - ptri[1].y * ptri[0].z;
		ptIn_[i].p0y0 = ptri[0].x * ptri[1].z - ptri[1].x * ptri[0].z;
		ptIn_[i].p0z0 = ptri[0].x * ptri[1].y - ptri[1].x * ptri[0].y;

		ptIn_[i].t00 = ptri[2].x * ptIn_[i].p0x0 -
			ptri[2].y * ptIn_[i].p0y0
			+ ptri[2].z * ptIn_[i].p0z0;

		ptIn_[i].p0x1 = ptri[0].y * ptri[2].z - ptri[2].y * ptri[0].z;
		ptIn_[i].p0y1 = ptri[0].x * ptri[2].z - ptri[2].x * ptri[0].z;
		ptIn_[i].p0z1 = ptri[0].x * ptri[2].y - ptri[2].x * ptri[0].y;

		ptIn_[i].t01 = ptri[1].x * ptIn_[i].p0x1 -
			ptri[1].y * ptIn_[i].p0y1
			+ ptri[1].z * ptIn_[i].p0z1;

		ptIn_[i].p0x2 = ptri[1].y * ptri[2].z - ptri[2].y * ptri[1].z;
		ptIn_[i].p0y2 = ptri[1].x * ptri[2].z - ptri[2].x * ptri[1].z;
		ptIn_[i].p0z2 = ptri[1].x * ptri[2].y - ptri[2].x * ptri[1].y;

		ptIn_[i].t02 = ptri[0].x * ptIn_[i].p0x2 -
			ptri[0].y * ptIn_[i].p0y2
			+ ptri[0].z * ptIn_[i].p0z2;
	}
}

int RtSphHexDGGS::GetSphTri(const GeoCoord& pt, const PreCompInTri* pptri)
{
	double p0;

	DVec3D pp = llxyz(pt);

	for (int i = 0; i < 60; i++)
	{
		const PreCompInTri& intri = pptri[i];

		p0 = pp.x * intri.p0x0 - pp.y * intri.p0y0 + pp.z * intri.p0z0;

		if (p0 * intri.t00 < 0.0)
			continue;

		p0 = pp.x * intri.p0x1 - pp.y * intri.p0y1 + pp.z * intri.p0z1;

		if (p0 * intri.t01 < 0.0)
			continue;

		p0 = pp.x * intri.p0x2 - pp.y * intri.p0y2 + pp.z * intri.p0z2;

		if (p0 * intri.t02 < 0.0)
			continue;

		return i;
	}

	return -1;
}

void RtSphHexDGGS::SetL()
{
	for (int i = 0; i < 32; i++)
	{
		desVerts_[i].x = R_ * cos(geoVerts_[i].lat) * cos(geoVerts_[i].lon);
		desVerts_[i].y = R_ * cos(geoVerts_[i].lat) * sin(geoVerts_[i].lon);
		desVerts_[i].z = R_ * sin(geoVerts_[i].lat);
	}

	for (int i = 0; i < 30; i++)
	{
		int* id_set = new int[4];
		id_set = rhomb_[i].GetVertsId();

		rhomb_[i].SetVertsgeo(geoVerts_[id_set[0]], geoVerts_[id_set[1]], geoVerts_[id_set[2]], geoVerts_[id_set[3]]);
		rhomb_[i].SetVerts3d(desVerts_[id_set[0]], desVerts_[id_set[1]], desVerts_[id_set[2]], desVerts_[id_set[3]]);
	}

	for (int i = 0; i < 60; i++)
	{
		int m;
		if (i % 2 == 0)
			m = i / 2;
		else
			m = (i - 1) / 2;

		int* id_set = new int[4];
		id_set = rhomb_[m].GetVertsId();

		if (i % 2 == 0)
		{
			tri_[i].SetVertsgeo(geoVerts_[id_set[0]], geoVerts_[id_set[1]], geoVerts_[id_set[2]]);
			tri_[i].SetVerts3d(desVerts_[id_set[0]], desVerts_[id_set[1]], desVerts_[id_set[2]]);
		}
		else
		{
			tri_[i].SetVertsgeo(geoVerts_[id_set[3]], geoVerts_[id_set[2]], geoVerts_[id_set[1]]);
			tri_[i].SetVerts3d(desVerts_[id_set[3]], desVerts_[id_set[2]], desVerts_[id_set[1]]);
		}
	}
}

DGridCoord2D RtSphHexDGGS::RTForwISlice(const GeoCoord& p)
{
	DGridCoord2D p2d;

	p2d.idx = GetSphTri(p, ptIn_);

	int* tri_id = tri_[p2d.idx].GetVertsID();

	GeoCoord B = geoVerts_[tri_id[0]];
	GeoCoord A = geoVerts_[tri_id[1]];
	GeoCoord C = geoVerts_[tri_id[2]];

	DVec2D B_(0.0, 1.164268433 / 2.0);
	DVec2D A_(-0.7195574638 / 2.0, 0.0);
	DVec2D C_(0.7195574638 / 2.0, 0.0);

	double b = 0.7195574638;
	double a = 1.164268433;

	double r = Intersec(B, A, p);

	if (fabs(r) < TOLERANCE)
	{
		DVec3D B_3d = llxyz(B);
		DVec3D P_3d = llxyz(p);
		double x = GreatCircleDistance(B_3d, P_3d);

		double m = sqrt((1 - cos(x)) / (1 - cos(0.652358139)));
		p2d.pt.x = A_.x * m;
		p2d.pt.y = a * (1 - m) / 2.0;
		return p2d;
	}

	double delta_cos = -cos(r) / 2.0 + M_SQRT3_2 * sin(r) * cos(0.652358139);
	double delta = acos(delta_cos);

	double m1 = 11 - (r + delta) * 15 / M_PI;

	DVec2D D_(0.0, 0.0);
	D_.x = b * (1 - 2 * m1) / 2.0;

	DVec3D B_3d = llxyz(B);
	DVec3D P_3d = llxyz(p);
	double x = GreatCircleDistance(B_3d, P_3d);

	double q = (0.5 + cos(r) * delta_cos) / (sin(r) * sin(delta));
	double m2 = sqrt((1 - cos(x)) / (1 - q));

	p2d.pt.x = D_.x * m2;
	p2d.pt.y = a * (1 - m2) / 2.0;

	p2d.pt.x *= R_;
	p2d.pt.y *= R_;

	return p2d;
}

GeoCoord RtSphHexDGGS::RTInvSlice(DGridCoord2D& p)
{
	DVec2D B_(0.0, 1.164268433 / 2.0);
	DVec2D A_(-0.7195574638 / 2.0, 0.0);
	DVec2D C_(0.7195574638 / 2.0, 0.0);

	int* id = tri_[p.idx].GetVertsID();

	DVec3D B = desVerts_[id[0]];
	DVec3D A = desVerts_[id[1]];
	DVec3D C = desVerts_[id[2]];

	B.x /= R_;  B.y /= R_;	B.z /= R_;
	A.x /= R_;  A.y /= R_;	A.z /= R_;
	C.x /= R_;  C.y /= R_;	C.z /= R_;

	p.pt.x /= R_; p.pt.y /= R_;

	double m;

	DVec2D D_(0.0, 0.0);
	DVec3D D;

	if (fabs(p.pt.x) < TOLERANCE)
	{
		m = 0.5;
	}
	else
	{
		double k = (p.pt.y - 0.5821342165) / p.pt.x;
		double x0 = -0.5821342165 / k;

		m = x0 / 0.7195574638 + 0.5;
		D_.x = x0;
	}

	double Y = (m / 15.0 + 2.0 / 3.0) * M_PI;
	double cos_Y = cos(Y);
	double sin_Y = sqrt(1 - cos_Y * cos_Y);

	double b = 0.729727656;

	double tan_r = (cos_Y + 0.5) / ((M_SQRT3_2 * COS_A) - sin_Y);
	double sin_r = tan_r / sqrt(1 + tan_r * tan_r);
	double cos_r = sin_r / tan_r;

	double sin_ad = sin_r * SIN_A / (sin_Y * cos_r - cos_Y * sin_r);
	double cos_ad = sqrt(1 - sin_ad * sin_ad);

	double l0 = SIN_B * cos_ad - COS_B * sin_ad;

	D.x = (C.x * sin_ad + A.x * l0) / SIN_B;
	D.y = (C.y * sin_ad + A.y * l0) / SIN_B;
	D.z = (C.z * sin_ad + A.z * l0) / SIN_B;

	double m2 = (B_.y - p.pt.y) * (B_.y - p.pt.y) / (B_.y * B_.y);

	double q = B.x * D.x + B.y * D.y + B.z * D.z;

	double cos_x = 1 - (1 - q) * m2;
	double sin_x = sqrt(1 - cos_x * cos_x);

	double bd = GreatCircleDistance(D, B);

	double sin_bd = sin(bd);
	double cos_bd = sqrt(1 - sin_bd * sin_bd);

	DVec3D p_3d;
	double q0 = 1.0 / sin_bd;

	double sin_bdx = sin_bd * cos_x - cos_bd * sin_x;

	p_3d.x = q0 * (sin_x * D.x + sin_bdx * B.x);
	p_3d.y = q0 * (sin_x * D.y + sin_bdx * B.y);
	p_3d.z = q0 * (sin_x * D.z + sin_bdx * B.z);

	GeoCoord P_geo = xyzll(p_3d);

	return P_geo;
}

DGridCoord2D RtSphHexDGGS::GetTricoordfromIJ(const DGridCoord2Dij input)
{
	if (input.idx < 0 || input.idx > 29)
	{
		cout << "the face ID is wrong!" << endl;
		exit(1);
	}

	int IMax = (int)pow(2, level_) - 1;
	int JMax = (int)pow(2, level_) - 1;
	int IMin = 0;
	int JMin = 0;
	switch (input.idx)
	{
	case 5:case 6:case 7:case 8:case 9:
	case 20:case 21:case 22:case 23:case 24:
	{
		JMin += 1;
		JMax += 1;
		break;
	}
	default:break;
	}

	if (input.pt.i<IMin || input.pt.i>IMax || input.pt.j<JMin || input.pt.j>JMax)
	{
		cout << "the I or J coordinate is wrong!" << endl;
		exit(1);
	}

	DGridCoord2D transCoord;
	transCoord.pt.x = (input.pt.i + input.pt.j) * lengthK_ / 2 - b_ / 2;
	transCoord.pt.y = (input.pt.j - input.pt.i) * lengthH_ / 2;

	transCoord.idx = input.idx * 2;
	if (transCoord.pt.y < 0)
	{
		transCoord.idx += 1;
		transCoord.pt.x = -transCoord.pt.x;
		transCoord.pt.y = -transCoord.pt.y;
	}
	return transCoord;
}

DGridCoord2Dij RtSphHexDGGS::GetIJfromTricoord(const DGridCoord2D& gpIn)
{
	DGridCoord2Dij resul;

	DGridCoord2D gp = gpIn;
	if (gp.idx % 2)
	{
		gp.pt.x = -gp.pt.x;
		gp.pt.y = -gp.pt.y;
	}
	gp.idx /= 2;
	resul.idx = gp.idx;

	gp.pt.x += b_ / 2;

	gp.pt.x = (gp.pt.x - gp.pt.y / TAN_HALF_ANG) / 2 / COS_HALF_ANG;
	gp.pt.y = gp.pt.x + gp.pt.y / SIN_HALF_ANG;
	if (gp.pt.x < 0 || gp.pt.y < 0)
	{
		cout << "there is something wrong with corI or corJ!" << endl;
		exit(1);
	}

	gp.pt.x /= lengthIJ_;
	gp.pt.y /= lengthIJ_;

	resul.pt.i = (int)gp.pt.x;
	resul.pt.j = (int)gp.pt.y;

	double resI = gp.pt.x - resul.pt.i;
	double resJ = gp.pt.y - resul.pt.j;

	if (resI < M_5)
	{
		if (resI < M_3)
		{
			if (resJ > resI * M_5 + M_5)
				resul.pt.j += 1;
		}
		else
		{
			if (resJ > 1 - resI)
				resul.pt.j += 1;
			if (resJ > 1 - resI && resJ < 2 * resI)
				resul.pt.i += 1;
		}
	}
	else if (resI < M2_3)
	{
		if (resJ > 1 - resI)
			resul.pt.j += 1;
		if (resJ > 1 - resI || resJ < 2 * resI - 1)
			resul.pt.i += 1;
	}
	else
	{
		resul.pt.i += 1;
		if (resJ > resI * M_5)
			resul.pt.j += 1;
	}

	switch (resul.idx)
	{
	case 5:case 6:case 7:case 8:case 9:
	{
		if (resul.pt.j == 0)
		{
			if (resul.pt.i < m_)
			{
				resul.idx = rhomb_[resul.idx].GetNeighFaceId()[0];
				resul.pt.j = resul.pt.i;
				resul.pt.i = 0;
			}
			else
			{
				resul.idx = rhomb_[rhomb_[resul.idx].GetNeighFaceId()[3]].GetNeighFaceId()[3];
				resul.pt.i = 0;
				resul.pt.j = m_;
			}
		}
		else if (resul.pt.i == m_)
		{
			resul.idx = rhomb_[resul.idx].GetNeighFaceId()[3];
			resul.pt.i = m_ - resul.pt.j;
			resul.pt.j = 0;
		}
		break;
	}

	case 20:case 21:case 22:case 23:case 24:
	{
		if (resul.pt.j == 0)
		{
			if (resul.pt.i < m_)
			{
				resul.idx = rhomb_[resul.idx].GetNeighFaceId()[0];
				resul.pt.j = resul.pt.i;
				resul.pt.i = 0;
			}
			else
			{
				resul.idx = rhomb_[rhomb_[resul.idx].GetNeighFaceId()[0]].GetNeighFaceId()[2];
				resul.pt.i = 0;
				resul.pt.j = m_;
			}
		}
		else if (resul.pt.i == m_)
		{
			resul.idx = rhomb_[resul.idx].GetNeighFaceId()[3];
			resul.pt.i = m_ - resul.pt.j;
			resul.pt.j = 0;
		}
		break;
	}

	case 11:case 13:case 15:case 17:case 19:
	{
		if (resul.pt.i == m_)
		{
			if (resul.pt.j == 0)
			{
				resul.idx = rhomb_[rhomb_[resul.idx].GetNeighFaceId()[3]].GetNeighFaceId()[3];
				resul.pt.i = 0;
				resul.pt.j = m_;
			}
			else
			{
				resul.idx = rhomb_[resul.idx].GetNeighFaceId()[3];
				resul.pt.i = m_ - resul.pt.j;
				resul.pt.j = 0;
			}
		}
		else if (resul.pt.j == m_)
		{
			resul.idx = rhomb_[resul.idx].GetNeighFaceId()[2];
			resul.pt.j = m_ - resul.pt.i;
			resul.pt.i = 0;
		}
		break;
	}

	case 10:case 12:case 14:case 16:case 18:
	{
		if (resul.pt.j == m_)
		{
			if (resul.pt.i == 0)
			{
				resul.idx = rhomb_[rhomb_[resul.idx].GetNeighFaceId()[2]].GetNeighFaceId()[2];
			}
			else
			{
				resul.idx = rhomb_[resul.idx].GetNeighFaceId()[2];
				resul.pt.j = m_ - resul.pt.i;
				resul.pt.i = 0;
			}
		}
		else if (resul.pt.i == m_)
		{
			resul.idx = rhomb_[resul.idx].GetNeighFaceId()[3];
			resul.pt.i = resul.pt.j;
			resul.pt.j = m_;
		}
		break;
	}

	case 0:case 1:case 2:case 3:case 4:
	{
		if (resul.pt.j == m_)
		{
			if (resul.pt.i == 0)
			{
				resul.idx = 30;
			}
			else
			{
				resul.idx = rhomb_[resul.idx].GetNeighFaceId()[2];
				resul.pt.j = m_ - resul.pt.i;
				resul.pt.i = 0;
			}
		}
		else if (resul.pt.i == m_)
		{
			resul.idx = rhomb_[rhomb_[resul.idx].GetNeighFaceId()[2]].GetNeighFaceId()[0];
			resul.pt.i = resul.pt.j;
			resul.pt.j = m_;
		}
		break;
	}

	case 25:case 26:case 27:case 28:case 29:
	{
		if (resul.pt.i == m_)
		{
			if (resul.pt.j == 0)
			{
				resul.idx = 31;
			}
			else
			{
				resul.idx = rhomb_[resul.idx].GetNeighFaceId()[3];
				resul.pt.i = m_ - resul.pt.j;
				resul.pt.j = 0;
			}
		}
		else if (resul.pt.j == m_)
		{
			resul.idx = rhomb_[resul.idx].GetNeighFaceId()[2];
			resul.pt.j = m_ - resul.pt.i;
			resul.pt.i = 0;
		}
		break;
	}

	default:cout << "something wrong with faceID! " << endl;
		break;
	}
	return resul;
}

GeoCoord RtSphHexDGGS::CalGeoFromij(const DGridCoord2Dij& input)
{
	if (input.idx == 30)
	{
		return GeoCoord(M_PI / 2.0, 0.0);
	}
	if (input.idx == 31)
	{
		return GeoCoord(-M_PI / 2.0, 0.0);
	}

	DGridCoord2D transCoord = GetTricoordfromIJ(input);

	return RTInvSlice(transCoord);
}

DGridCoord2Dij RtSphHexDGGS::CalijFromGeo(const GeoCoord& geo)
{
	DGridCoord2D gp = RTForwISlice(geo);
	return GetIJfromTricoord(gp);
}