#include <sstream>
#include <iostream>
#include <assert.h>

#include "georaster.h"

CGeoRaster::CGeoRaster()
{
	initialize();
}

CGeoRaster::CGeoRaster(const std::string filename)
{
	initialize();
	readMetaData(filename);
}

void  CGeoRaster::initialize()
{
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");

	_support_filename = "";
	_image_filename = "";
	_image_id = "";

	_nodata_value = -32768;
	_is_dem = false;

	_width = 0;
	_height = 0;
	_bands = 1;
	_datatype = GDT_Byte;

	_min_gray = NULL;
	_max_gray = NULL;

	_data = NULL;
	_loaded_in_memory = false;

	ground_point_struct gp;
	gp.x = gp.y = gp.z = 0;
	for (int i = 0; i < 4; i++) {
		_footprints[i] = gp;
		_footprints_latlon[i] = gp;
	}

	memset(_geotransform, 0, sizeof(double) * 6);

	_refCS.SetWellKnownGeogCS("WGS84");
	if (GDAL_VERSION_MAJOR >= 3)
		_refCS.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
}

CGeoRaster::CGeoRaster(const CGeoRaster& dom)
{
	if (this != &dom) {
		_support_filename = dom._support_filename;
		_image_filename = dom._image_filename;

		_width = dom._width;
		_height = dom._height;
		_bands = dom._bands;
		_datatype = dom._datatype;

		if (_min_gray) {
			delete[]_min_gray;
			_min_gray = NULL;
		}

		if (_max_gray) {
			delete[]_max_gray;
			_max_gray = NULL;
		}

		_min_gray = new double[_bands];
		_max_gray = new double[_bands];
		memcpy(_min_gray, dom._min_gray, sizeof(double) * _bands);
		memcpy(_max_gray, dom._max_gray, sizeof(double) * _bands);

		memcpy(_footprints, dom._footprints, sizeof(ground_point_struct) * 4);
		memcpy(_footprints_latlon, dom._footprints_latlon, sizeof(ground_point_struct) * 4);

		memcpy(_geotransform, dom._geotransform, sizeof(double) * 6);

		_refCS = dom._refCS;
	}
}

CGeoRaster::~CGeoRaster()
{
	if (_min_gray) {
		delete[]_min_gray;
		_min_gray = NULL;
	}

	if (_max_gray) {
		delete[]_max_gray;
		_max_gray = NULL;
	}

	if (_data != NULL) {
		CPLFree(_data);
		_data = NULL;
	}
}

bool CGeoRaster::isPointInside(const double& lon, const double& lat)
{
	double samp, line;
	ground2Image(lon, lat, samp, line);
	return samp > 0 && samp < _width && line > 0 && line < _height;
}

bool CGeoRaster::loadInMemory(const int band)
{
	GDALDataset* poDomDS;
	poDomDS = (GDALDataset*)GDALOpen(_image_filename.c_str(), GA_ReadOnly);
	if (poDomDS == NULL) {
		GDALClose(poDomDS);
		std::cerr << "Open raster file to read failed: " << _image_filename;
		return false;
	}

	if (_data != NULL) {
		CPLFree(_data);
		_data = NULL;
	}

	GDALRasterBand* pSrcBand;
	pSrcBand = poDomDS->GetRasterBand(band);
	int iDepth = GDALGetDataTypeSize(_datatype) / 8;

	_data = CPLMalloc(_width * _height * iDepth);
	if (_data == NULL) {
		std::cerr << "Read raster file failed, no sufficient memory. " << _image_filename;
		GDALClose(poDomDS);
		return false;
	}

	if (CE_Failure == pSrcBand->RasterIO(GF_Read, 0, 0, _width, _height, _data, _width, _height, _datatype, 0, 0)) {
		std::cerr << "Read raster file failed. " << _image_filename;
		GDALClose((GDALDatasetH)poDomDS);
		return false;
	}

	GDALClose(poDomDS);

	_loaded_in_memory = true;

	return true;
}

void CGeoRaster::ground2Image(const double& ex, const double& ny, double& samp, double& line)
{
	samp = _geotransform[5] * (ex - _geotransform[0]) - _geotransform[2] * (ny - _geotransform[3]);
	line = _geotransform[1] * (ny - _geotransform[3]) - _geotransform[4] * (ex - _geotransform[0]);

	samp /= _geotransform[1] * _geotransform[5] - _geotransform[2] * _geotransform[4];
	line /= _geotransform[1] * _geotransform[5] - _geotransform[2] * _geotransform[4];
}

void CGeoRaster::image2Ground(const double& samp, const double& line, double& ex, double& ny)
{
	ex = _geotransform[0] + _geotransform[1] * samp + _geotransform[2] * line;
	ny = _geotransform[3] + _geotransform[4] * samp + _geotransform[5] * line;
}

bool CGeoRaster::readMetaData(const std::string filename)
{
	GDALDataset* poDomDS;
	poDomDS = (GDALDataset*)GDALOpen(filename.c_str(), GA_ReadOnly);
	if (poDomDS == NULL) {
		GDALClose(poDomDS);
		std::cerr << "Open raster file to read failed: " << filename;
		return false;
	}

	_width = static_cast<size_t>(poDomDS->GetRasterXSize());
	_height = static_cast<size_t>(poDomDS->GetRasterYSize());
	_bands = static_cast<size_t>(poDomDS->GetRasterCount());
	_datatype = poDomDS->GetRasterBand(1)->GetRasterDataType();

	poDomDS->GetGeoTransform(_geotransform);

	image2Ground(0, 0, _footprints[0].x, _footprints[0].y);
	image2Ground(_width - 1, 0, _footprints[1].x, _footprints[1].y);
	image2Ground(_width - 1, _height - 1, _footprints[2].x, _footprints[2].y);
	image2Ground(0, _height - 1, _footprints[3].x, _footprints[3].y);

	if (_min_gray) {
		delete[]_min_gray;
		_min_gray = NULL;
	}

	if (_max_gray) {
		delete[]_max_gray;
		_max_gray = NULL;
	}

	_min_gray = new double[_bands];
	_max_gray = new double[_bands];

	int bGotMin, bGotMax;
	double adfMinMax[2];
	GDALRasterBand* pSrcBand;
	for (int i = 0; i < _bands; ++i) {
		pSrcBand = poDomDS->GetRasterBand(i + 1);
		_min_gray[i] = pSrcBand->GetMinimum(&bGotMin);
		_max_gray[i] = pSrcBand->GetMaximum(&bGotMax);
		if (!(bGotMax && bGotMin)) {
			GDALComputeRasterMinMax(pSrcBand, true, adfMinMax);
			_min_gray[i] = adfMinMax[0];
			_max_gray[i] = adfMinMax[1];
		}
	}

	char* pszProjection = NULL;
	if (poDomDS->GetProjectionRef() != NULL) {
		pszProjection = (char*)poDomDS->GetProjectionRef();
		if (_refCS.importFromWkt(&pszProjection) == CE_None) {
			CPLDebug("gdalsrsinfo", "got SRS from GDAL ");
		}
	}
	else {
		delete[] _max_gray;
		delete[] _min_gray;
		GDALClose(poDomDS);
		std::cerr << "The projection of geo raster was not defined.  " << filename;
		return false;
	}

	GDALClose(poDomDS);

	OGRSpatialReference    oLatLong;
	if (GDAL_VERSION_MAJOR >= 3)
		oLatLong.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);

	oLatLong.SetWellKnownGeogCS("WGS84");

	OGRCoordinateTransformation* poTransDem2LatLon = OGRCreateCoordinateTransformation(&_refCS, &oLatLong);
	assert(poTransDem2LatLon != nullptr);

	memcpy(_footprints_latlon, _footprints, sizeof(ground_point_struct) * 4);
	for (int i = 0; i < 4; ++i)
		poTransDem2LatLon->Transform(1, &_footprints_latlon[i].x, &_footprints_latlon[i].y);

	OGRCoordinateTransformation::DestroyCT(poTransDem2LatLon);

	_image_id = StringGetFileName(filename);
	_image_filename = filename;
	return true;
}

bool CGeoRaster::computeFootprints()
{
	if (_width == 0 || _height == 0)
		return false;

	if (_geotransform[1] < 1e-7 || std::abs(_geotransform[5]) < 1e-7)
		return false;

	OGRSpatialReference    oLatLong;
	if (GDAL_VERSION_MAJOR >= 3)
		oLatLong.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);

	oLatLong.SetWellKnownGeogCS("WGS84");

	OGRCoordinateTransformation* poTransDem2LatLon = OGRCreateCoordinateTransformation(&_refCS, &oLatLong);
	if (poTransDem2LatLon == nullptr)
		return false;

	image2Ground(0, 0, _footprints[0].x, _footprints[0].y);
	image2Ground(_width - 1, 0, _footprints[1].x, _footprints[1].y);
	image2Ground(_width - 1, _height - 1, _footprints[2].x, _footprints[2].y);
	image2Ground(0, _height - 1, _footprints[3].x, _footprints[3].y);

	memcpy(_footprints_latlon, _footprints, sizeof(ground_point_struct) * 4);
	for (int i = 0; i < 4; ++i)
		poTransDem2LatLon->Transform(1, &_footprints_latlon[i].x, &_footprints_latlon[i].y);

	OGRCoordinateTransformation::DestroyCT(poTransDem2LatLon);
	return true;
}

bool CGeoRaster::readBlock(const int ulRow, const int ulCol, const int ww, const int hh, void* data, GDALDataType datatype, const int band)
{
	GDALDataset* poDomDS;
	poDomDS = (GDALDataset*)GDALOpen(_image_filename.c_str(), GA_ReadOnly);
	if (poDomDS == NULL) {
		std::cerr << "Read Geo Raster Failed." << _image_filename;
		GDALClose(poDomDS);
		return false;
	}

	GDALRasterBand* poBand;
	poBand = poDomDS->GetRasterBand(band);

	if (poBand->RasterIO(GF_Read, ulCol, ulRow, ww, hh, data, ww, hh, datatype, 0, 0) == CE_Failure) {
		std::cerr << "Read Geo Raster File Failed." << _image_filename;
		GDALClose(poDomDS);
		return false;
	}

	GDALClose(poDomDS);
	return true;
}