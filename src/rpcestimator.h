#pragma once

#include "common.h"
#include <vector>
#include <string>
#include <memory>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include"RtA4HexDGGS.h"

using namespace Eigen;

/*
 * @brief	 Extended RPC generation function.
 *
 * @param inputImage: The path of the original image.
 * @param rpbFile: The path of the original rpb file.
 * @param rpbexFile: The path of the output extended rpbex file.
 *
 * @return	 0 for success and 1 for failure.
 */
bool rpcmake(std::string, std::string, std::string);

class TifRPCDataSet;
class RpcEstimator
{
public:
	RpcEstimator(std::shared_ptr<TifRPCDataSet> sm, double minz, double maxz) {
		_sm = sm; _min_z = minz; _max_z = maxz;
		_NumLine.resize(20);
		_DenLine.resize(20);
		_NumSamp.resize(20);
		_DenSamp.resize(20);
	};

	RpcEstimator(std::shared_ptr<TifRPCDataSet> sm) {
		_sm = sm; _min_z = 0; _max_z = 1000;
		_NumLine.resize(20);
		_DenLine.resize(20);
		_NumSamp.resize(20);
		_DenSamp.resize(20);
	};

	RpcEstimator() {
		_sm = nullptr; _min_z = 0; _max_z = 1000;
		_NumLine.resize(20);
		_DenLine.resize(20);
		_NumSamp.resize(20);
		_DenSamp.resize(20);
	};

	void solve();

	bool  solveLSM(const std::vector<ground_point_struct>& vecgps,
		const std::vector<image_point_struct>& vecips);

	void  setSensorModel(std::shared_ptr<TifRPCDataSet> sm) { _sm = sm; }
	void  setZRange(const double zmin, const double zmax) { _min_z = zmin; _max_z = zmax; }

	void  generateGrids(const int grid,
		const int lyrs,
		std::vector<ground_point_struct>& vecgps,
		std::vector<image_point_struct>& vecips);

	bool  toRPB(std::string rpcfilename);

	double rmsError() { return _rms_error; };
	double maxError() { return _max_error; };

protected:
	void	normalize(const std::vector<ground_point_struct>& vecgps,
		const std::vector<image_point_struct>& vecips);

	double  polyTrans(double nv, double nu, double nw, double coe[20]);

	void    evalErrors();

	void updateMatrix(const std::vector<double>& f,
		const std::vector<double>& x,
		const std::vector<double>& y,
		const std::vector<double>& z,
		MatrixXd& M);

	void updateWeights(const std::vector<double>& x,
		const std::vector<double>& y,
		const std::vector<double>& z,
		const VectorXd& coeff,
		VectorXd& W);

	void computeCoefficients(const std::vector<double>& f,
		const std::vector<double>& x,
		const std::vector<double>& y,
		const std::vector<double>& z,
		VectorXd& outCoeffs);

public:
	double  lineOffset() { return _line_offset; }
	double  sampOffset() { return _samp_offset; }
	double  longOffset() { return _lon_offset; }
	double  latOffset() { return _lat_offset; }
	double  heightOffset() { return _hgt_offset; }

	double  lineScale() { return _line_scale; }
	double  sampScale() { return _samp_scale; }
	double  longScale() { return _lon_scale; }
	double  latScale() { return _lat_scale; }
	double  heightScale() { return _hgt_scale; }

	const std::vector<double>& lineNumCoefficients() { return _NumLine; }
	const std::vector<double>& lineDenCoefficients() { return _DenLine; }
	const std::vector<double>& sampNumCoefficients() { return _NumSamp; }
	const std::vector<double>& sampDenCoefficients() { return _DenSamp; }

protected:
	std::vector<double> _latNorm;
	std::vector<double> _lonNorm;
	std::vector<double> _hgtNorm;
	std::vector<double> _sampNorm;
	std::vector<double> _lineNorm;

	std::vector<double> _NumLine;
	std::vector<double> _DenLine;
	std::vector<double> _NumSamp;
	std::vector<double> _DenSamp;

	std::shared_ptr<TifRPCDataSet> _sm;

	double  _min_z;
	double  _max_z;

	double  _max_error;
	double  _rms_error;

	double	_line_offset;
	double	_samp_offset;
	double  _lon_offset;
	double  _lat_offset;
	double  _hgt_offset;

	double	_line_scale;
	double	_samp_scale;
	double  _lon_scale;
	double  _lat_scale;
	double  _hgt_scale;
};
