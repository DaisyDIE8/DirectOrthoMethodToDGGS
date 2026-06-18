#include "rpcestimator.h"
#include "tifrpcdataset.h"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>

void   RpcEstimator::normalize(
	const std::vector<ground_point_struct>& vecgps,
	const std::vector<image_point_struct>& vecips
) {
	assert(vecips.size() == vecgps.size());

	int numberOfPoints = vecips.size();

	_latNorm.clear();
	_lonNorm.clear();
	_hgtNorm.clear();
	_sampNorm.clear();
	_lineNorm.clear();

	_latNorm.reserve(numberOfPoints);
	_lonNorm.reserve(numberOfPoints);
	_hgtNorm.reserve(numberOfPoints);
	_sampNorm.reserve(numberOfPoints);
	_lineNorm.reserve(numberOfPoints);

	double accLat = 0.;
	double accLon = 0.;
	double accAlt = 0.;

	for (auto i = 0; i < numberOfPoints; i++) {
		accLon += vecgps[i].x;
		accLat += vecgps[i].y;
		accAlt += vecgps[i].z;
	}

	_lon_offset = accLon / numberOfPoints;
	_lat_offset = accLat / numberOfPoints;
	_hgt_offset = accAlt / numberOfPoints;

	double minc = std::numeric_limits<double>::max();
	double minl = std::numeric_limits<double>::max();
	double maxc = std::numeric_limits<double>::min();
	double maxl = std::numeric_limits<double>::min();

	for (auto i = 0; i < numberOfPoints; i++) {
		minc = std::min(vecips[i].sample, minc);
		maxc = std::max(vecips[i].sample, maxc);

		minl = std::min(vecips[i].line, minl);
		maxl = std::max(vecips[i].line, maxl);
	}

	_samp_offset = (minc + maxc - 1) / 2.0;
	_line_offset = (minl + maxl - 1) / 2.0;

	auto height = std::abs(minl - maxl) + 1;
	auto width = std::abs(minc - maxc) + 1;

	_samp_scale = width / 2.0;
	_line_scale = height / 2.0;

	double maxDeltaLon = std::numeric_limits<double>::min();
	double maxDeltaLat = std::numeric_limits<double>::min();
	double maxDeltaAlt = std::numeric_limits<double>::min();
	for (auto i = 0; i < numberOfPoints; i++) {
		auto deltaLon = vecgps[i].x - _lon_offset;
		auto deltaLat = vecgps[i].y - _lat_offset;
		auto deltaAlt = vecgps[i].z - _hgt_offset;

		maxDeltaLon = std::max(maxDeltaLon, std::abs(deltaLon));
		maxDeltaLat = std::max(maxDeltaLat, std::abs(deltaLat));
		maxDeltaAlt = std::max(maxDeltaAlt, std::abs(deltaAlt));
	}

	if (maxDeltaLat < 1.0)
		maxDeltaLat = 1.0;

	if (maxDeltaLon < 1.0)
		maxDeltaLon = 1.0;

	if (maxDeltaAlt < 1.0)
		maxDeltaAlt = 1.0;

	_lat_scale = maxDeltaLat;
	_lon_scale = maxDeltaLon;
	_hgt_scale = maxDeltaAlt;

	for (auto i = 0; i < numberOfPoints; i++) {
		auto deltaLon = vecgps[i].x - _lon_offset;
		auto deltaLat = vecgps[i].y - _lat_offset;
		auto deltaAlt = vecgps[i].z - _hgt_offset;

		_sampNorm.push_back((vecips[i].sample - _samp_offset) / _samp_scale);
		_lineNorm.push_back((vecips[i].line - _line_offset) / _line_scale);
		_lonNorm.push_back(deltaLon / _lon_scale);
		_latNorm.push_back(deltaLat / _lat_scale);
		_hgtNorm.push_back(deltaAlt / _hgt_scale);
	}
}

bool  RpcEstimator::solveLSM(
	const std::vector<ground_point_struct>& vecgps,
	const std::vector<image_point_struct>& vecips)
{
	normalize(vecgps, vecips);

	solve();

	return true;
}

void RpcEstimator::updateMatrix(const std::vector<double>& f,
	const std::vector<double>& x,
	const std::vector<double>& y,
	const std::vector<double>& z,
	MatrixXd& M)
{
	for (unsigned int i = 0; i < f.size(); i++)
	{
		M(i, 0) = 1;
		M(i, 1) = x[i];
		M(i, 2) = y[i];
		M(i, 3) = z[i];
		M(i, 4) = x[i] * y[i];
		M(i, 5) = x[i] * z[i];
		M(i, 6) = y[i] * z[i];
		M(i, 7) = x[i] * x[i];
		M(i, 8) = y[i] * y[i];
		M(i, 9) = z[i] * z[i];
		M(i, 10) = x[i] * y[i] * z[i];
		M(i, 11) = x[i] * x[i] * x[i];
		M(i, 12) = x[i] * y[i] * y[i];
		M(i, 13) = x[i] * z[i] * z[i];
		M(i, 14) = x[i] * x[i] * y[i];
		M(i, 15) = y[i] * y[i] * y[i];
		M(i, 16) = y[i] * z[i] * z[i];
		M(i, 17) = x[i] * x[i] * z[i];
		M(i, 18) = y[i] * y[i] * z[i];
		M(i, 19) = z[i] * z[i] * z[i];
		M(i, 20) = -f[i] * x[i];
		M(i, 21) = -f[i] * y[i];
		M(i, 22) = -f[i] * z[i];
		M(i, 23) = -f[i] * x[i] * y[i];
		M(i, 24) = -f[i] * x[i] * z[i];
		M(i, 25) = -f[i] * y[i] * z[i];
		M(i, 26) = -f[i] * x[i] * x[i];
		M(i, 27) = -f[i] * y[i] * y[i];
		M(i, 28) = -f[i] * z[i] * z[i];
		M(i, 29) = -f[i] * x[i] * y[i] * z[i];
		M(i, 30) = -f[i] * x[i] * x[i] * x[i];
		M(i, 31) = -f[i] * x[i] * y[i] * y[i];
		M(i, 32) = -f[i] * x[i] * z[i] * z[i];
		M(i, 33) = -f[i] * x[i] * x[i] * y[i];
		M(i, 34) = -f[i] * y[i] * y[i] * y[i];
		M(i, 35) = -f[i] * y[i] * z[i] * z[i];
		M(i, 36) = -f[i] * x[i] * x[i] * z[i];
		M(i, 37) = -f[i] * y[i] * y[i] * z[i];
		M(i, 38) = -f[i] * z[i] * z[i] * z[i];
	}
}

void RpcEstimator::updateWeights(const std::vector<double>& x,
	const std::vector<double>& y,
	const std::vector<double>& z,
	const VectorXd& coeff,
	VectorXd& W)
{
	std::vector<double> row(coeff.size());

	for (unsigned int i = 0; i < x.size(); i++)
	{
		W(i) = 0.;

		row[0] = 1;
		row[1] = x[i];
		row[2] = y[i];
		row[3] = z[i];
		row[4] = x[i] * y[i];
		row[5] = x[i] * z[i];
		row[6] = y[i] * z[i];
		row[7] = x[i] * x[i];
		row[8] = y[i] * y[i];
		row[9] = z[i] * z[i];
		row[10] = x[i] * y[i] * z[i];
		row[11] = x[i] * x[i] * x[i];
		row[12] = x[i] * y[i] * y[i];
		row[13] = x[i] * z[i] * z[i];
		row[14] = x[i] * x[i] * y[i];
		row[15] = y[i] * y[i] * y[i];
		row[16] = y[i] * z[i] * z[i];
		row[17] = x[i] * x[i] * z[i];
		row[18] = y[i] * y[i] * z[i];
		row[19] = z[i] * z[i] * z[i];

		for (unsigned int j = 0; j < row.size(); j++)
		{
			W(i) += row[j] * coeff(j);
		}

		if (W(i) > epsilon)
		{
			W(i) = 1.0 / W(i);
		}
	}
}

void RpcEstimator::computeCoefficients(const std::vector<double>& f,
	const std::vector<double>& x,
	const std::vector<double>& y,
	const std::vector<double>& z,
	VectorXd& outCoeffs)
{
	size_t numberOfPoints = f.size();

	MatrixXd M(numberOfPoints, 39);
	VectorXd r(f.size());
	for (auto i = 0; i < f.size(); i++)
		r(i) = f[i];

	VectorXd weights(numberOfPoints);
	for (auto i = 0; i < numberOfPoints; i++)
		weights(i) = 1.0;

	double res = std::numeric_limits<double>::max();
	for (int i = 0; i < 10; i++)
	{
		if (res < epsilon)
		{
			break;
		}

		Eigen::SparseMatrix<double> spW(numberOfPoints, numberOfPoints);
		for (auto j = 0; j < numberOfPoints; j++) {
			spW.coeffRef(j, j) = weights(j) * weights(j);
		}

		updateMatrix(f, x, y, z, M);

		Eigen::MatrixXd  A = M.transpose() * spW * M;
		Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);

		Eigen::MatrixXd V = svd.matrixV(), U = svd.matrixU();
		Eigen::VectorXd diag = svd.singularValues();

		for (unsigned int idx = 0; idx < diag.size(); idx++)
		{
			if (diag[idx] > epsilon)
			{
				diag[idx] = 1.0 / diag[idx];
			}
			else
			{
				diag[idx] = 0.0;
			}
		}

		Eigen::MatrixXd invA = V * diag.asDiagonal() * U.transpose();

		outCoeffs = invA * M.transpose() * spW * r;

		Eigen::VectorXd denominator = outCoeffs.segment(19, 20);
		denominator(0) = 1.;

		updateWeights(x, y, z, denominator, weights);

		Eigen::VectorXd  residual = M.transpose() * spW * (M * outCoeffs - r);

		res = residual.norm();
	}
}

void RpcEstimator::solve()
{
	auto numberOfPoints = _latNorm.size();

	if (numberOfPoints < 20) {
		std::cerr << "At least 20 points are required to estimate the parameters of a RPC model. Only " << numberOfPoints << " points were given." << std::endl;
		return;
	}

	Eigen::VectorXd colCoeffs;
	computeCoefficients(_sampNorm, _lonNorm, _latNorm, _hgtNorm, colCoeffs);

	Eigen::VectorXd lineCoeffs;
	computeCoefficients(_lineNorm, _lonNorm, _latNorm, _hgtNorm, lineCoeffs);

	for (auto i = 0; i < 20; i++) {
		_NumLine[i] = lineCoeffs(i);
		_NumSamp[i] = colCoeffs(i);
	}

	_DenLine[0] = 1.0;
	_DenSamp[0] = 1.0;
	for (auto i = 0; i < 19; i++) {
		_DenLine[i + 1] = lineCoeffs(i + 20);
		_DenSamp[i + 1] = colCoeffs(i + 20);
	}

	evalErrors();
}

void  RpcEstimator::generateGrids(const int grid,
	const int lyrs,
	std::vector<ground_point_struct>& vecgps,
	std::vector<image_point_struct>& vecips) {
	size_t ww = _sm->imageSamples();
	size_t hh = _sm->imageLines();

	int cols = static_cast<int>(ceil(ww * 1.0 / grid));
	int rows = static_cast<int>(ceil(hh * 1.0 / grid));

	for (auto i = 0; i <= rows; i++) {
		image_point_struct im;
		im.line = i * grid;
		if (im.line >= (hh - 1))
			im.line = hh - 1;

		for (auto j = 0; j <= cols; j++) {
			im.sample = j * grid;
			if (im.sample >= (ww - 1))
				im.sample = (ww - 1);

			for (auto k = 0; k < lyrs; k++) {
				double zz = _min_z + k * (_max_z - _min_z) / (lyrs - 1);

				ground_point_struct gp;
				gp.z = zz;

				_sm->vImgToGrnd(im, gp);
				vecgps.push_back(gp);
				vecips.push_back(im);
			}
		}
	}
}

double  RpcEstimator::polyTrans(double nv, double nu, double nw, double coe[20])
{
	double L[20];
	double Fi = 0;

	L[0] = 1;
	L[1] = nv;
	L[2] = nu;
	L[3] = nw;
	L[4] = nv * nu;
	L[5] = nv * nw;
	L[6] = nu * nw;
	L[7] = nv * nv;
	L[8] = nu * nu;
	L[9] = nw * nw;
	L[10] = nu * nv * nw;
	L[11] = nv * nv * nv;
	L[12] = nv * nu * nu;
	L[13] = nv * nw * nw;
	L[14] = nv * nv * nu;
	L[15] = nu * nu * nu;
	L[16] = nu * nw * nw;
	L[17] = nv * nv * nw;
	L[18] = nu * nu * nw;
	L[19] = nw * nw * nw;

	for (int i = 0; i < 20; i++)
		Fi += L[i] * coe[i];

	return Fi;
}

void    RpcEstimator::evalErrors()
{
	double _maxError = 0.0;
	double _rmsError = 0.0;
	for (auto i = 0; i < _latNorm.size(); i++)
	{
		double F1 = polyTrans(_lonNorm[i], _latNorm[i], _hgtNorm[i], _NumSamp.data());
		double F2 = polyTrans(_lonNorm[i], _latNorm[i], _hgtNorm[i], _DenSamp.data());
		double sampout = F1 / F2 * _samp_scale + _samp_offset;

		double F3 = polyTrans(_lonNorm[i], _latNorm[i], _hgtNorm[i], _NumLine.data());
		double F4 = polyTrans(_lonNorm[i], _latNorm[i], _hgtNorm[i], _DenLine.data());
		double lineout = F3 / F4 * _line_scale + _line_offset;

		double samp = _sampNorm[i] * _samp_scale + _samp_offset;
		double line = _lineNorm[i] * _line_scale + _line_offset;

		double lensqr = (samp - sampout) * (samp - sampout) + (line - lineout) * (line - lineout);

		_rmsError += lensqr;
		if (sqrt(lensqr) > _maxError)
			_maxError = sqrt(lensqr);
	}

	_rmsError = std::sqrt(_rmsError / _latNorm.size());
	_max_error = _maxError;
	_rms_error = _rmsError;
}

bool  RpcEstimator::toRPB(std::string rpcfilename)
{
	std::ofstream os;
	os.open(rpcfilename);
	if (os.fail())
		return false;

	os << "satId = \"NOT_ASSIGNED\";\n";
	os << "bandId = \"NOT_ASSIGNED\";\n";
	os << "SpecId = \"RPC00B\";\n";

	os << "BEGIN_GROUP = IMAGE\n";

	char buffer[255];
	sprintf(buffer, "\t%16s %+22.15f;\n", "errBias = ", _max_error);
	os << buffer;
	sprintf(buffer, "\t%16s %+22.15f;\n", "errRand = ", _rms_error);
	os << buffer;
	sprintf(buffer, "\t%16s %+22.15f;\n", "lineOffset = ", _line_offset);
	os << buffer;
	sprintf(buffer, "\t%16s %+22.15f;\n", "sampOffset = ", _samp_offset);
	os << buffer;
	sprintf(buffer, "\t%16s %+22.15f;\n", "latOffset = ", _lat_offset);
	os << buffer;
	sprintf(buffer, "\t%16s %+22.15f;\n", "longOffset = ", _lon_offset);
	os << buffer;
	sprintf(buffer, "\t%16s %+22.15f;\n", "heightOffset = ", _hgt_offset);
	os << buffer;

	sprintf(buffer, "\t%16s %+22.15f;\n", "lineScale = ", _line_scale);
	os << buffer;
	sprintf(buffer, "\t%16s %+22.15f;\n", "sampScale = ", _samp_scale);
	os << buffer;
	sprintf(buffer, "\t%16s %+22.15f;\n", "latScale = ", _lat_scale);
	os << buffer;
	sprintf(buffer, "\t%16s %+22.15f;\n", "longScale = ", _lon_scale);
	os << buffer;
	sprintf(buffer, "\t%16s %+22.15f;\n", "heightScale = ", _hgt_scale);
	os << buffer;

	os << "\tlineNumCoef = (\n\t\t\t";
	for (int i = 0; i < 20; ++i) {
		if (i > 0) os << "\t\t\t";
		sprintf(buffer, "%+22.15E", _NumLine[i]);
		os << buffer;
		if (i < 19) os << ";\n";
	}
	os << ");\n";

	os << "\tlineDenCoef = (\n\t\t\t";
	for (int i = 0; i < 20; ++i) {
		if (i > 0) os << "\t\t\t";
		sprintf(buffer, "%+22.15E", _DenLine[i]);
		os << buffer;
		if (i < 19) os << ";\n";
	}
	os << ");\n";

	os << "\tsampNumCoef = (\n\t\t\t";
	for (int i = 0; i < 20; ++i) {
		if (i > 0) os << "\t\t\t";
		sprintf(buffer, "%+22.15E", _NumSamp[i]);
		os << buffer;
		if (i < 19) os << ";\n";
	}
	os << ");\n";

	os << "\tsampDenCoef = (\n\t\t\t";
	for (int i = 0; i < 20; ++i) {
		if (i > 0) os << "\t\t\t";
		sprintf(buffer, "%+22.15E", _DenSamp[i]);
		os << buffer;
		if (i < 19) os << ";\n";
	}
	os << ");\n";

	os << "END_GROUP = IMAGE\n";
	os << "END;";

	os.close();
	return true;
}

void transformGridPointsToTSEA(std::vector<ground_point_struct>& vecGPs)
{
	for (auto i = 0; i < vecGPs.size(); i++) {
		GeoCoord gpt(vecGPs[i].y * DEG_TO_RAD, vecGPs[i].x * DEG_TO_RAD);

		const int level_ = 18;
		RtSphHexDGGS imaDGGS = RtSphHexDGGS(6371, level_);
		DGridCoord2D ixy = imaDGGS.RTForwISlice(gpt);

		ground_point_struct gp;
		gp.x = ixy.pt.x * 1000;
		gp.y = ixy.pt.y * 1000;
		gp.z = vecGPs[i].z;

		vecGPs[i] = gp;
	}
}

bool rpcmake(std::string inputImage, std::string rpbFile, std::string rpbexFile) {
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");

	std::shared_ptr<TifRPCDataSet> tifRpcDS = std::make_shared<TifRPCDataSet>();
	RPCSensorModel* sm = tifRpcDS->sensorModel();

	if (!tifRpcDS->loadImage(inputImage) && !tifRpcDS->readRPB(rpbFile)) {
		std::cerr << "Load source image failed: " << inputImage << std::endl;
		return EXIT_FAILURE;
	}
	else {
		tifRpcDS->readRPB(rpbFile);
		std::cout << "Source image loaded: " << inputImage << std::endl;
	}

	OGRSpatialReference eoSRS;
	if (GDAL_VERSION_MAJOR >= 3)
		eoSRS.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);

	double zmin = 0.0, zmax = 1000.0;

	double iZMin = zmin, iZMax = zmax;

	std::shared_ptr<RpcEstimator> rpcSolver = std::make_shared<RpcEstimator>();
	rpcSolver->setSensorModel(tifRpcDS);
	rpcSolver->setZRange(iZMin, iZMax);

	std::vector<ground_point_struct> vecgps;
	std::vector<image_point_struct> vecips;
	rpcSolver->generateGrids(500, 5, vecgps, vecips);

	GeoCoord gpt;
	gpt.lon = 0;
	gpt.lat = 0;

	const ground_point_struct* footprints = tifRpcDS->footPrints();
	for (int i = 0; i < 4; i++) {
		gpt.lon += footprints[i].x * DEG_TO_RAD;
		gpt.lat += footprints[i].y * DEG_TO_RAD;
	}

	gpt.lon /= 4;
	gpt.lat /= 4;

	transformGridPointsToTSEA(vecgps);

	std::cout << "Solve new rpc parameters......";
	if (rpcSolver->solveLSM(vecgps, vecips)) {
		std::cout << "finished." << std::endl;

		rpcSolver->toRPB(rpbexFile);
	}

	return 0;
}