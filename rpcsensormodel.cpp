#include "rpcsensormodel.h"

#include <fstream>
#include <sstream>
#include <cstring>
#include <math.h>

RPCSensorModel::RPCSensorModel()
{
	_line_offset = 0;
	_samp_offset = 0;
	_lat_offset = 0;
	_long_offset = 0;
	_height_offset = 0;

	_line_scale = 1;
	_samp_scale = 1;
	_lat_scale = 1;
	_long_scale = 1;
	_height_scale = 1;

	memset(_NumLine, 0, sizeof(double) * 20);
	memset(_DenLine, 0, sizeof(double) * 20);
	memset(_NumSamp, 0, sizeof(double) * 20);
	memset(_DenSamp, 0, sizeof(double) * 20);

	_rpcFileName = "";
}

RPCSensorModel::RPCSensorModel(double* rfcs)
{
	RFCsParse(rfcs);
	_rpcFileName = "";
}

RPCSensorModel::~RPCSensorModel()
{
}

void  RPCSensorModel::RFCsParse(double* rfcs)
{
	_line_offset = rfcs[0];
	_samp_offset = rfcs[1];
	_lat_offset = rfcs[2];
	_long_offset = rfcs[3];
	_height_offset = rfcs[4];

	_line_scale = rfcs[5];
	_samp_scale = rfcs[6];
	_lat_scale = rfcs[7];
	_long_scale = rfcs[8];
	_height_scale = rfcs[9];

	for (int i = 0; i < 20; i++)
	{
		_NumLine[i] = rfcs[i + 10];
		_DenLine[i] = rfcs[i + 30];
		_NumSamp[i] = rfcs[i + 50];
		_DenSamp[i] = rfcs[i + 70];
	}
}

bool  RPCSensorModel::readRPC(const std::string rpcFileName)
{
	char line[256], dummy[64];
	std::string aline;

	std::ifstream is;
	is.open(rpcFileName);
	if (is.fail())
		return false;

	double  rpc[90];
	for (int i = 0; i < 90; i++)
	{
		is.getline(line, 256);
		sscanf(line, "%s %lf", dummy, &rpc[i]);
	}
	is.close();

	RFCsParse(rpc);
	_rpcFileName = rpcFileName;

	return true;
}

bool  RPCSensorModel::readRPB(const std::string  rpcFileName)
{
	std::ifstream is;
	is.open(rpcFileName);
	if (is.fail())
		return false;

	char dummy[50], equmark[8], aline[256];
	int i;

	while (!is.eof())
	{
		is.getline(aline, 256);

		if (strstr(aline, "lineOffset"))
		{
			sscanf(aline, "%s %s %lf", dummy, equmark, &_line_offset);
		}
		else if (strstr(aline, "sampOffset"))
		{
			sscanf(aline, "%s %s %lf", dummy, equmark, &_samp_offset);
		}
		else if (strstr(aline, "latOffset"))
		{
			sscanf(aline, "%s %s %lf", dummy, equmark, &_lat_offset);
		}
		else if (strstr(aline, "longOffset"))
		{
			sscanf(aline, "%s %s %lf", dummy, equmark, &_long_offset);
		}
		else if (strstr(aline, "heightOffset"))
		{
			sscanf(aline, "%s %s %lf", dummy, equmark, &_height_offset);
		}
		else if (strstr(aline, "lineScale"))
		{
			sscanf(aline, "%s %s %lf", dummy, equmark, &_line_scale);
		}
		else if (strstr(aline, "sampScale"))
		{
			sscanf(aline, "%s %s %lf", dummy, equmark, &_samp_scale);
		}
		else if (strstr(aline, "latScale"))
		{
			sscanf(aline, "%s %s %lf", dummy, equmark, &_lat_scale);
		}
		else if (strstr(aline, "longScale"))
		{
			sscanf(aline, "%s %s %lf", dummy, equmark, &_long_scale);
		}
		else if (strstr(aline, "heightScale"))
		{
			sscanf(aline, "%s %s %lf", dummy, equmark, &_height_scale);
		}
		else if (strstr(aline, "lineNumCoef"))
		{
			for (i = 0; i < 20; i++)
			{
				is.getline(aline, 256);
				sscanf(aline, "%lf%c", &_NumLine[i], equmark);
			}
		}
		else if (strstr(aline, "lineDenCoef"))
		{
			for (i = 0; i < 20; i++)
			{
				is.getline(aline, 256);
				sscanf(aline, "%lf%c", &_DenLine[i], equmark);
			}
		}
		else if (strstr(aline, "sampNumCoef"))
		{
			for (i = 0; i < 20; i++)
			{
				is.getline(aline, 256);
				sscanf(aline, "%lf%c", &_NumSamp[i], equmark);
			}
		}
		else if (strstr(aline, "sampDenCoef"))
		{
			for (i = 0; i < 20; i++)
			{
				is.getline(aline, 256);
				sscanf(aline, "%lf%c", &_DenSamp[i], equmark);
			}
		}
	}

	is.close();
	_rpcFileName = rpcFileName;
	return true;
}

bool  RPCSensorModel::Conv_Map_To_SL(const double& lon, const double& lat, const double& hz, double* samp, double* line)
{
	double tempU, tempV, tempW;
	double L[20], tempx, tempy;
	double F1, F2, F3, F4;

	tempV = (lon - _long_offset) / _long_scale;
	tempU = (lat - _lat_offset) / _lat_scale;
	tempW = (hz - _height_offset) / _height_scale;

	L[0] = 1;
	L[1] = tempV;
	L[2] = tempU;
	L[3] = tempW;
	L[4] = tempV * tempU;
	L[5] = tempV * tempW;
	L[6] = tempU * tempW;
	L[7] = tempV * tempV;
	L[8] = tempU * tempU;
	L[9] = tempW * tempW;
	L[10] = tempU * tempV * tempW;
	L[11] = tempV * tempV * tempV;
	L[12] = tempV * tempU * tempU;
	L[13] = tempV * tempW * tempW;
	L[14] = tempV * tempV * tempU;
	L[15] = tempU * tempU * tempU;
	L[16] = tempU * tempW * tempW;
	L[17] = tempV * tempV * tempW;
	L[18] = tempU * tempU * tempW;
	L[19] = tempW * tempW * tempW;

	F1 = 0; F2 = 0; F3 = 0; F4 = 0;
	for (int i = 0; i < 20; i++)
	{
		F1 += L[i] * _NumLine[i];
		F2 += L[i] * _DenLine[i];
		F3 += L[i] * _NumSamp[i];
		F4 += L[i] * _DenSamp[i];
	}

	tempy = F1 / F2;
	tempx = F3 / F4;

	*samp = tempx * _samp_scale + _samp_offset;
	*line = tempy * _line_scale + _line_offset;
	return true;
}

double  RPCSensorModel::polyTrans(const double& nv, const double& nu, const double& nw, double coe[20])
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

double  RPCSensorModel::polyLinearize2Lon(double coe[20], const double& nv, const double& nu, const double& nw)
{
	double dF_V;
	dF_V = coe[1] + coe[4] * nu + coe[5] * nw + 2 * coe[7] * nv + coe[10] * nu * nw + 3 * coe[11] * nv * nv + coe[12] * nu * nu +
		coe[13] * nw * nw + 2 * coe[14] * nv * nu + 2 * coe[17] * nv * nw;

	return dF_V;
}

double  RPCSensorModel::polyLinearize2Lat(double coe[20], const double& nv, const double& nu, const double& nw)
{
	double dF_U;
	dF_U = coe[2] + coe[4] * nv + coe[6] * nw + 2 * coe[8] * nu + coe[10] * nv * nw + 2 * coe[12] * nv * nu + coe[14] * nv * nv +
		3 * coe[15] * nu * nu + coe[16] * nw * nw + 2 * coe[18] * nu * nw;

	return dF_U;
}

bool RPCSensorModel::Conv_SL_To_Map(const double& samp, const double& line, const double& hz, double* lon, double* lat)
{
	double tempLon = _long_offset;
	double tempLat = _lat_offset;
	double tempW = (hz - _height_offset) / _height_scale;
	double deltaMax = 1e8;

	int    iterCnt = 0;

	do
	{
		double tempV = (tempLon - _long_offset) / _long_scale;
		double tempU = (tempLat - _lat_offset) / _lat_scale;

		double F1, F2, F3, F4;
		F1 = polyTrans(tempV, tempU, tempW, _NumLine);
		F2 = polyTrans(tempV, tempU, tempW, _DenLine);
		F3 = polyTrans(tempV, tempU, tempW, _NumSamp);
		F4 = polyTrans(tempV, tempU, tempW, _DenSamp);

		double tempy = F1 / F2;
		double tempx = F3 / F4;

		double dF1_U, dF1_V, dF2_U, dF2_V, dF3_U, dF3_V, dF4_U, dF4_V;
		dF1_V = polyLinearize2Lon(_NumLine, tempV, tempU, tempW);
		dF1_U = polyLinearize2Lat(_NumLine, tempV, tempU, tempW);

		dF2_V = polyLinearize2Lon(_DenLine, tempV, tempU, tempW);
		dF2_U = polyLinearize2Lat(_DenLine, tempV, tempU, tempW);

		dF3_V = polyLinearize2Lon(_NumSamp, tempV, tempU, tempW);
		dF3_U = polyLinearize2Lat(_NumSamp, tempV, tempU, tempW);

		dF4_V = polyLinearize2Lon(_DenSamp, tempV, tempU, tempW);
		dF4_U = polyLinearize2Lat(_DenSamp, tempV, tempU, tempW);

		double a0, a1, b0, b1, l[2];
		a0 = (dF1_V * F2 - F1 * dF2_V) * _line_scale / ((F2 * F2) * _long_scale);
		a1 = (dF1_U * F2 - F1 * dF2_U) * _line_scale / ((F2 * F2) * _lat_scale);
		b0 = (dF3_V * F4 - F3 * dF4_V) * _samp_scale / ((F4 * F4) * _long_scale);
		b1 = (dF3_U * F4 - F3 * dF4_U) * _samp_scale / ((F4 * F4) * _lat_scale);

		l[0] = line - (tempy * _line_scale + _line_offset);
		l[1] = samp - (tempx * _samp_scale + _samp_offset);

		double n11, n12, n21, n22, u0, u1;
		n11 = a0 * a0 + b0 * b0;
		n12 = a0 * a1 + b0 * b1;
		n21 = n12;
		n22 = a1 * a1 + b1 * b1;

		u0 = a0 * l[0] + b0 * l[1];
		u1 = a1 * l[0] + b1 * l[1];

		double deltalon = (n22 * u0 - n12 * u1) / (n11 * n22 - n21 * n12);
		double deltalat = (-n21 * u0 + n11 * u1) / (n11 * n22 - n21 * n12);

		tempLon += deltalon;
		tempLat += deltalat;

		deltaMax = fabs(deltalon) > fabs(deltalat) ? fabs(deltalon) : fabs(deltalat);
		if (deltaMax < 1e-8)
			break;

		iterCnt++;
	} while (iterCnt < 10);

	*lat = tempLat;
	*lon = tempLon;
	return true;
}