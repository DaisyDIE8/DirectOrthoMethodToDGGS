#pragma once
#include "stdafx.h"
#include "gdal_priv.h"

const double epsilon = 0.000000000000001;

struct ground_point_struct
{
	bool operator== (const ground_point_struct& gp)
	{
		if (fabs(gp.x - x) > 1e-8 || fabs(gp.y - y) > 1e-8 || fabs(gp.z - z) > 1e-3)
			return false;
		else
			return true;
	}

	double x;
	double y;
	double z;
};

struct image_point_struct
{
	bool operator== (const image_point_struct& im)
	{
		if (fabs(im.sample - sample) > 1e-6 || fabs(im.line - line) > 1e-6)
			return false;
		else
			return true;
	}

	double line;
	double sample;
};

typedef struct
{
	int col;
	int row;
	double gridX;
	double gridY;
	double Xmin;
	double Xmax;
	double Ymin;
	double Ymax;
	double Zmin;
	double Zmax;
}DEM_GRD;

extern bool  InterpolateZ(const double& EX, const double& NY, double* HZ, const DEM_GRD& deminfo, float* dem);

extern bool StringContains(const std::string& str, const std::string& sub_str);

extern std::string StringReplace(const std::string& str, const std::string& old_str,
	const std::string& new_str);

extern std::string StringGetFileName(const std::string& path);