#pragma once

#include <vector>
#include <string>
#include <memory>

#include "common.h"

enum RFMType
{
	eRPB = 0,
	eRPC
};

class RPCSensorModel;
class TifRPCDataSet
{
public:
	TifRPCDataSet();
	TifRPCDataSet(const TifRPCDataSet& ds);
	virtual ~TifRPCDataSet();

	bool  loadImage(const std::string tiffile);
	bool  readRPB(const std::string rpbFile);

	void  setImageFileName(const std::string filename) { _image_filename = filename; }
	void  setRpcFile(const std::string filename) { _rpc_filename = filename; }
	void  setMeanTerrainHeight(const double hgt) { _mean_terrain_height = hgt; }

	std::string  imageFileName() const { return _image_filename; }
	std::string  imageId() const { return _image_id; }
	std::string  rpcFileName() const { return _rpc_filename; }
	size_t  imageSamples()const { return _image_samples; }
	size_t  imageLines() const { return _image_lines; }
	double  meanTerrainHeight()const { return _mean_terrain_height; }
	double  level0GSD() const { return _lev0gsd; }
	bool    isRPB() const { return _rfmType == eRPB ? true : false; }
	RPCSensorModel* sensorModel()const { return _sm; }

	const ground_point_struct* footPrints() const { return _footprints; };

public:
	virtual bool vGrndToImg(const ground_point_struct& gp, image_point_struct& out);

	virtual bool vImgToGrnd(const image_point_struct& img_pt_in, ground_point_struct& gp_out);

	bool  computeFootPrints();

protected:
	std::string _image_filename;
	std::string _image_id;
	std::string _rpc_filename;

	size_t    _image_lines;
	size_t    _image_samples;

	double _mean_terrain_height;
	double _lev0gsd;

	ground_point_struct    _footprints[4];

	RFMType  _rfmType;
	RPCSensorModel* _sm;
};
