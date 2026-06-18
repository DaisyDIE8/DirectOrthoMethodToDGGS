#include <iostream>
#include <vector>
#include <list>

#include "gdal.h"
#include "gdal_priv.h"
#include "rpcsensormodel.h"
#include "tifrpcdataset.h"

TifRPCDataSet::TifRPCDataSet()
{
	_image_filename = "";
	_rpc_filename = "";
	_image_id = "";
	_image_samples = 0;
	_image_lines = 0;
	_mean_terrain_height = 0.0;
	_lev0gsd = 0.0;

	_rfmType = eRPB;

	ground_point_struct gp;
	gp.x = gp.y = gp.z = 0;
	for (int i = 0; i < 4; i++)
		_footprints[i] = gp;

	_sm = new RPCSensorModel;
}

TifRPCDataSet::TifRPCDataSet(const TifRPCDataSet& rhs)
{
	_image_filename = rhs._image_filename;
	_rpc_filename = rhs._rpc_filename;
	_image_id = rhs._image_id;
	_image_samples = rhs._image_samples;
	_image_lines = rhs._image_lines;
	_mean_terrain_height = rhs._mean_terrain_height;
	_lev0gsd = rhs._lev0gsd;

	memcpy(_footprints, rhs._footprints, sizeof(ground_point_struct) * 4);

	RPCSensorModel* sm = rhs.sensorModel();
	_sm = new RPCSensorModel(*sm);
}

TifRPCDataSet::~TifRPCDataSet()
{
	if (!_sm)
	{
		delete _sm;
		_sm = nullptr;
	}
}

bool TifRPCDataSet::loadImage(const std::string tiffile)
{
	GDALDataset* poTifDS;
	poTifDS = (GDALDataset*)GDALOpen(tiffile.c_str(), GA_ReadOnly);
	if (poTifDS == NULL) {
		GDALClose(poTifDS);
		std::cerr << "Open raster file to read failed: " << tiffile;
		return false;
	}

	_image_filename = tiffile;
	_image_id = StringGetFileName(_image_filename);
	_image_samples = static_cast<size_t>(poTifDS->GetRasterXSize());
	_image_lines = static_cast<size_t>(poTifDS->GetRasterYSize());
	GDALClose(poTifDS);

	std::string rpcFile, rpbFile;

	size_t idx = tiffile.find_last_of('.');
	std::string tifBasePath = tiffile.substr(0, idx);
	rpcFile = tifBasePath + "_rpc.txt";
	rpbFile = tifBasePath + ".rpb";

	if (_sm->readRPB(rpbFile)) {
		_rpc_filename = rpbFile;
		_rfmType = eRPB;
	}
	else if (_sm->readRPC(rpcFile)) {
		_rpc_filename = rpcFile;
		_rfmType = eRPC;
	}
	else
		return false;

	return computeFootPrints();
}

bool TifRPCDataSet::readRPB(const std::string rpbFile)
{
	if (_sm->readRPB(rpbFile)) {
		_rpc_filename = rpbFile;
		_rfmType = eRPB;
	}
	else
		return false;

	_mean_terrain_height = _sm->heightOffset();
	return computeFootPrints();
}

bool TifRPCDataSet::computeFootPrints()
{
	ground_point_struct gp;
	gp.z = _mean_terrain_height;

	ground_point_struct footprints[4];

	if (!_sm->Conv_SL_To_Map(0, 0, _mean_terrain_height, &gp.x, &gp.y))
		return false;
	footprints[0] = gp;

	if (!_sm->Conv_SL_To_Map(_image_samples, 0, _mean_terrain_height, &gp.x, &gp.y))
		return false;
	footprints[1] = gp;

	if (!_sm->Conv_SL_To_Map(_image_samples, _image_lines, _mean_terrain_height, &gp.x, &gp.y))
		return false;
	footprints[2] = gp;

	if (!_sm->Conv_SL_To_Map(0, _image_lines, _mean_terrain_height, &gp.x, &gp.y))
		return false;
	footprints[3] = gp;

	double minLon = 9.0e10, minLat = 9.0e10, maxLon = -9.0e10, maxLat = -9.0e10;

	double avgLon = 0, avgLat = 0;
	for (int i = 0; i < 4; ++i)
	{
		avgLon += footprints[i].x;
		avgLat += footprints[i].y;

		if (footprints[i].x < minLon)
			minLon = footprints[i].x;
		if (footprints[i].x > maxLon)
			maxLon = footprints[i].x;

		if (footprints[i].y < minLat)
			minLat = footprints[i].y;
		if (footprints[i].y > maxLat)
			maxLat = footprints[i].y;
	}

	avgLon /= 4;
	avgLat /= 4;

	_lev0gsd = (maxLon - minLon) * (maxLat - minLat) / (_image_samples * _image_lines);
	_lev0gsd = sqrt(_lev0gsd);

	for (int i = 0; i < 4; ++i)
	{
		if (footprints[i].x < avgLon && footprints[i].y > avgLat)
			_footprints[0] = footprints[i];

		if (footprints[i].x > avgLon && footprints[i].y > avgLat)
			_footprints[1] = footprints[i];

		if (footprints[i].x > avgLon && footprints[i].y < avgLat)
			_footprints[2] = footprints[i];

		if (footprints[i].x < avgLon && footprints[i].y < avgLat)
			_footprints[3] = footprints[i];
	}
	return true;
}

bool TifRPCDataSet::vGrndToImg(const ground_point_struct& gp, image_point_struct& out)
{
	return _sm->Conv_Map_To_SL(gp.x, gp.y, gp.z, &out.sample, &out.line);
}

bool TifRPCDataSet::vImgToGrnd(const image_point_struct& img_pt_in, ground_point_struct& gp_out)
{
	return _sm->Conv_SL_To_Map(img_pt_in.sample, img_pt_in.line, gp_out.z, &gp_out.x, &gp_out.y);
}