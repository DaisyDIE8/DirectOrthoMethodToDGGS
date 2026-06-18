#pragma once
#pragma warning(disable:4996)

#include <string>

class RPCSensorModel
{
public:
	RPCSensorModel();
	RPCSensorModel(double* rfcs);
	virtual ~RPCSensorModel();

	bool	readRPC(const std::string rpcFileName);
	bool	readRPB(const std::string rpcFileName);

protected:
	double	_NumLine[20], _DenLine[20], _NumSamp[20], _DenSamp[20];

	double	_line_offset;
	double	_samp_offset;
	double  _long_offset;
	double  _lat_offset;
	double  _height_offset;

	double	_line_scale;
	double	_samp_scale;
	double  _long_scale;
	double  _lat_scale;
	double  _height_scale;

	std::string _rpcFileName;

protected:
	double  polyTrans(const double& nv, const double& nu, const double& nw, double coe[20]);
	void    RFCsParse(double* rfcs);

	double  polyLinearize2Lon(double coe[20], const double& nv, const double& nu, const double& nw);
	double  polyLinearize2Lat(double coe[20], const double& nv, const double& nu, const double& nw);

public:
	bool    Conv_Map_To_SL(const double& lon, const double& lat, const double& hz, double* samp, double* line);
	bool    Conv_SL_To_Map(const double& samp, const double& line, const double& hz, double* lon, double* lat);

public:
	double  lineOffset() { return _line_offset; }
	double  sampOffset() { return _samp_offset; }
	double  longOffset() { return _long_offset; }
	double  latOffset() { return _lat_offset; }
	double  heightOffset() { return _height_offset; }

	double  lineScale() { return _line_scale; }
	double  sampScale() { return _samp_scale; }
	double  longScale() { return _long_scale; }
	double  latScale() { return _lat_scale; }
	double  heightScale() { return _height_scale; }

	double* lineNumCoefficients() { return _NumLine; }
	double* lineDenCoefficients() { return _DenLine; }
	double* sampNumCoefficients() { return _NumSamp; }
	double* sampDenCoefficients() { return _DenSamp; }
};
