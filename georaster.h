#pragma once

#include <vector>
#include <string>

#include "gdal.h"
#include "gdal_priv.h"
#include "common.h"

class CGeoRaster
{
public:
	CGeoRaster();
	CGeoRaster(const CGeoRaster& dom);
	CGeoRaster(const std::string filename);
	virtual ~CGeoRaster();

	bool readBlock(const int ulRow, const int ulCol, const int ww, const int hh, void* data, GDALDataType datatype, const int band = 1);
	bool readMetaData(const std::string filename);

	bool loadInMemory(const int band = 1);

	void setInvalidValue(int val) { _nodata_value = val; };
	void setDemRaster(bool isdem) { _is_dem = isdem; };
	int  isDemRaster() { return _is_dem; }

private:
	void initialize();

public:
	std::string supportFileName() const { return _support_filename; };
	std::string imageFileName() const { return _image_filename; };
	std::string imageID() const { return _image_id; };

	size_t  width() const { return _width; };
	size_t  height() const { return _height; };
	size_t  bands() const { return _bands; };
	GDALDataType  dataType() const { return _datatype; };
	int  invalidValue() const { return _nodata_value; };

	double* minGray() { return _min_gray; };
	double* maxGray() { return _max_gray; };
	const double* geoTransform() { return _geotransform; }

	void setSupportFileName(const std::string filename) { _support_filename = filename; }
	void setImageFileName(const std::string filename) { _image_filename = filename; }
	void setImageID(std::string imageid) { _image_id = imageid; };

	void setWidth(const size_t& ww) { _width = ww; }
	void setHeight(const size_t& hh) { _height = hh; }
	void setBands(const size_t& bands) { _bands = bands; }
	void setDataType(const GDALDataType& type) { _datatype = type; }
	void setData(void* data) { _data = data; _loaded_in_memory = true; }
	void setGeoTransform(const double geotrans[6]) { memcpy(_geotransform, geotrans, 6 * sizeof(double)); };

	void upperLeftCorner(double& ulx, double& uly) const { ulx = _footprints[0].x; uly = _footprints[0].y; };
	void gsd(double& gsdx, double& gsdy) const { gsdx = _geotransform[1]; gsdy = _geotransform[5]; };
	ground_point_struct* footPrints() { return _footprints; };
	ground_point_struct* footPrintsWGS84() { return _footprints_latlon; };

	OGRSpatialReference referenceCS() const { return _refCS; };
	void setReferenceCS(const OGRSpatialReference cs) { _refCS = cs; };

	void ground2Image(const double& ex, const double& ny, double& samp, double& line);
	void image2Ground(const double& samp, const double& line, double& ex, double& ny);

	bool isPointInside(const double& lon, const double& lat);

	bool computeFootprints();

protected:
	std::string _support_filename;
	std::string _image_filename;
	std::string _image_id;

	int  _nodata_value;
	int  _is_dem;

	size_t  _width;
	size_t  _height;
	size_t  _bands;
	GDALDataType  _datatype;

	double* _min_gray;
	double* _max_gray;

	bool    _loaded_in_memory;
	void* _data;

	double  _geotransform[6];

	ground_point_struct    _footprints[4];
	ground_point_struct    _footprints_latlon[4];

	OGRSpatialReference _refCS;
};
