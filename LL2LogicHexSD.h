#pragma once

#include "RtA4HexDGGS.h"
#include <gdal.h>
#include "ogrsf_frmts.h"
#include <algorithm>

#include <iostream>
#include <fstream>
#include <sstream>

#include<assert.h>
#include <string.h>

bool ll2logicHexSD(std::string srcImagePath, std::string outFilename, std::string prjFile, InterpolateType t)
{
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");

	GDALDataset* poSrcDS = (GDALDataset*)GDALOpen(srcImagePath.c_str(), GA_ReadOnly);
	if (poSrcDS == NULL) {
		std::cout << "Source image can not be opened:" << srcImagePath << std::endl;
		return false;
	}
	else {
		std::cout << "Source image loaded: " << srcImagePath << std::endl;
	}

	int nWidth = poSrcDS->GetRasterXSize();
	int nHeight = poSrcDS->GetRasterYSize();
	int nBands = poSrcDS->GetRasterCount();
	GDALDataType nDataType = poSrcDS->GetRasterBand(1)->GetRasterDataType();
	int nDepth = GDALGetDataTypeSize((GDALDataType)nDataType) / 8;

	double srcImageTrans[6];
	poSrcDS->GetGeoTransform(srcImageTrans);

	double dLon = srcImageTrans[1] * DEG_TO_RAD;
	double dLat = srcImageTrans[5] * DEG_TO_RAD;
	double lonBegin = srcImageTrans[0] * DEG_TO_RAD + dLon * 0.5;
	double latBegin = srcImageTrans[3] * DEG_TO_RAD + dLat * 0.5;

	if (nDataType != GDT_Byte && nDataType != GDT_UInt16 && nDataType != GDT_Int16 && nDataType != GDT_Float32) {
		std::cout << "Unsupported data bits for source image." << std::endl;
		GDALClose(poSrcDS);
		return false;
	}

	float* pSrcData;
	pSrcData = (float*)VSIMalloc2(nWidth, nHeight * nBands * sizeof(float));
	if (pSrcData == NULL) {
		GDALClose(poSrcDS);
		return false;
	}

	poSrcDS->RasterIO(GF_Read,
		0, 0,
		nWidth, nHeight,
		(void*)pSrcData,
		nWidth, nHeight,
		GDT_Float32, nBands, NULL, 0, 0, 0);

	std::cout << "Source image read as float successfully " << std::endl;
	GDALClose(poSrcDS);

	GDALDriver* poDriver;
	poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
	if (poDriver == NULL) {
		CPLFree(pSrcData);
		return false;
	}

	const int level_ = 18;
	RtSphHexDGGS imaDGGS = RtSphHexDGGS(6371, level_);

	std::ifstream is;
	is.open(prjFile);
	if (is.fail()) {
		std::cout << "prj failed" << endl;
		return EXIT_FAILURE;
	}

	std::stringstream buffer;
	buffer << is.rdbuf();
	std::string srsWKT(buffer.str());
	is.close();

	unsigned short  rtFaceId = -1;
	std::string tesa;
	std::stringstream ss(srsWKT);
	ss >> tesa >> rtFaceId;
	assert(rtFaceId >= 0 && rtFaceId < 60);

	double orthLogicHexX = imaDGGS.lengthK_ / 2 * 1000;
	double orthLogicHexY = imaDGGS.lengthH_ / 2 * 1000;

	DVec2D footPr[4];
	footPr[0] = imaDGGS.RTForwISlice(GeoCoord(
		latBegin,
		lonBegin
	)).pt;
	footPr[1] = imaDGGS.RTForwISlice(GeoCoord(
		latBegin,
		lonBegin + dLon * nWidth
	)).pt;
	footPr[2] = imaDGGS.RTForwISlice(GeoCoord(
		latBegin + dLat * nHeight,
		lonBegin
	)).pt;
	footPr[3] = imaDGGS.RTForwISlice(GeoCoord(
		latBegin + dLat * nHeight,
		lonBegin + dLon * nWidth
	)).pt;

	for (auto i = 0; i < 4; i++) {
		footPr[i].x *= 1000;
		footPr[i].y *= 1000;
	}

	double orthoXmin = numeric_limits<double>::max();
	double orthoYmin = numeric_limits<double>::max();
	double orthoXmax = -numeric_limits<double>::max();
	double orthoYmax = -numeric_limits<double>::max();

	for (auto i = 0; i < 4; i++) {
		orthoXmin = std::min(footPr[i].x, orthoXmin);
		orthoXmax = std::max(footPr[i].x, orthoXmax);

		orthoYmin = std::min(footPr[i].y, orthoYmin);
		orthoYmax = std::max(footPr[i].y, orthoYmax);
	}

	DGridCoord2Dij first = imaDGGS.GetIJfromTricoord(DGridCoord2D(rtFaceId, orthoXmin * 0.001, orthoYmax * 0.001));
	orthoXmin = imaDGGS.GetTricoordfromIJ(first).pt.x * 1000 - orthLogicHexX;
	orthoYmax = imaDGGS.GetTricoordfromIJ(first).pt.y * 1000;

	DGridCoord2Dij last = imaDGGS.GetIJfromTricoord(DGridCoord2D(rtFaceId, orthoXmax * 0.001, orthoYmin * 0.001));
	orthoXmax = imaDGGS.GetTricoordfromIJ(last).pt.x * 1000 + orthLogicHexX;
	orthoYmin = imaDGGS.GetTricoordfromIJ(last).pt.y * 1000;

	int nOrthoXSize = static_cast<int>((orthoXmax - orthoXmin) / orthLogicHexX + 0.5);
	if (nOrthoXSize % 4 != 0)
	{
		nOrthoXSize = 4 * (nOrthoXSize / 4 + 1);
		orthoXmax = orthoXmin + nOrthoXSize * orthLogicHexX;
	}
	int nOrthoYSize = static_cast<int>((orthoYmax - orthoYmin) / orthLogicHexY + 0.5);

	void* pOrthoData;
	pOrthoData = (void*)VSIMalloc2(nOrthoXSize, nOrthoYSize * nBands * nDepth);
	if (pOrthoData == NULL) {
		VSIFree(pSrcData);
		return false;
	}

	std::cout << "Start to reproject to SD plane..." << std::endl;

	float* g = new float[nBands];
	int onePerc = static_cast<int>(nOrthoYSize * 0.01), nPrec = 0;

#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < nOrthoYSize; i++)
	{
		double northing = orthoYmax - orthLogicHexY * i;
		for (int j = 0; j < nOrthoXSize; j += 2)
		{
			int m = j;
			if (i % 2) {
				m = j + 1;
				if (j == nOrthoXSize - 2) { continue; }
			}

			double easting = orthoXmin + (m + 1) * orthLogicHexX;
			DGridCoord2D ccxy(rtFaceId, easting * 0.001, northing * 0.001);
			GeoCoord tempGeo = imaDGGS.RTInvSlice(ccxy);

			double samp = (tempGeo.lon - lonBegin) / dLon;
			double line = -(latBegin - tempGeo.lat) / dLat;

			int ii = static_cast<int>(line);
			int jj = static_cast<int>(samp);

			double g1, g2, g3, g4;
			memset(g, 0, sizeof(float) * nBands);
			if (jj >= 1 && jj < nWidth - 1 && ii >= 1 && ii < nHeight - 1)
			{
				double dx = samp - jj;
				double dy = line - ii;

				for (int k = 0; k < nBands; ++k)
				{
					size_t pos = k * nWidth * nHeight;
					g1 = pSrcData[pos + ii * nWidth + jj];
					g2 = pSrcData[pos + ii * nWidth + jj + 1];
					g3 = pSrcData[pos + (ii + 1) * nWidth + jj + 1];
					g4 = pSrcData[pos + (ii + 1) * nWidth + jj];

					switch (t)
					{
					case BILINEAR:
						g[k] = (1 - dx) * (1 - dy) * g1 + dx * (1 - dy) * g2 + dx * dy * g3 + (1 - dx) * dy * g4;
						break;
					case NEAREST:
						if (dx <= 0.5)
						{
							if (dy <= 0.5) { g[k] = g1; }
							else { g[k] = g4; }
						}
						else
						{
							if (dy <= 0.5) { g[k] = g2; }
							else { g[k] = g3; }
						}
						break;
					default:
						break;
					}
				}
			}
			else
			{
				for (int k = 0; k < nBands; ++k)
					g[k] = 0;
			}

			if (nDataType == GDT_Byte) {
				unsigned char* pByteData = (unsigned char*)pOrthoData;
				for (int k = 0; k < nBands; ++k)
				{
					size_t pos = k * nOrthoYSize * nOrthoXSize;
					pByteData[pos + i * nOrthoXSize + m] = static_cast<unsigned char>(g[k]);
					pByteData[pos + i * nOrthoXSize + m + 1] = static_cast<unsigned char>(g[k]);
				}
			}
			else if (nDataType == GDT_UInt16) {
				unsigned short* pUShortData = (unsigned short*)pOrthoData;
				for (int k = 0; k < nBands; ++k)
				{
					size_t pos = k * nOrthoYSize * nOrthoXSize;
					pUShortData[pos + i * nOrthoXSize + m] = static_cast<unsigned short>(g[k]);
					pUShortData[pos + i * nOrthoXSize + m + 1] = static_cast<unsigned short>(g[k]);
				}
			}
			else if (nDataType == GDT_Int16) {
				short* pShortData = (short*)pOrthoData;
				for (int k = 0; k < nBands; ++k)
				{
					size_t pos = k * nOrthoYSize * nOrthoXSize;
					pShortData[pos + i * nOrthoXSize + m] = static_cast<short>(g[k]);
					pShortData[pos + i * nOrthoXSize + m + 1] = static_cast<short>(g[k]);
				}
			}
			else if (nDataType == GDT_Float32) {
				float* pFloatData = (float*)pOrthoData;
				for (int k = 0; k < nBands; ++k)
				{
					size_t pos = k * nOrthoYSize * nOrthoXSize;
					pFloatData[pos + i * nOrthoXSize + m] = (g[k]);
					pFloatData[pos + i * nOrthoXSize + m + 1] = (g[k]);
				}
			}
		}

		if (i % 2)
		{
			if (nDataType == GDT_Byte) {
				unsigned char* pByteData = (unsigned char*)pOrthoData;
				for (int k = 0; k < nBands; ++k)
				{
					size_t pos = k * nOrthoYSize * nOrthoXSize;
					pByteData[pos + i * nOrthoXSize] = static_cast<unsigned char>(0);
					pByteData[pos + i * nOrthoXSize + nOrthoXSize - 1] = static_cast<unsigned char>(0);
				}
			}
			else if (nDataType == GDT_UInt16) {
				unsigned short* pUShortData = (unsigned short*)pOrthoData;
				for (int k = 0; k < nBands; ++k)
				{
					size_t pos = k * nOrthoYSize * nOrthoXSize;
					pUShortData[pos + i * nOrthoXSize] = static_cast<unsigned short>(0);
					pUShortData[pos + i * nOrthoXSize + nOrthoXSize - 1] = static_cast<unsigned short>(0);
				}
			}
			else if (nDataType == GDT_Int16) {
				short* pShortData = (short*)pOrthoData;
				for (int k = 0; k < nBands; ++k)
				{
					size_t pos = k * nOrthoYSize * nOrthoXSize;
					pShortData[pos + i * nOrthoXSize] = static_cast<short>(0);
					pShortData[pos + i * nOrthoXSize + nOrthoXSize - 1] = static_cast<short>(0);
				}
			}
			else if (nDataType == GDT_Float32) {
				float* pFloatData = (float*)pOrthoData;
				for (int k = 0; k < nBands; ++k)
				{
					size_t pos = k * nOrthoYSize * nOrthoXSize;
					pFloatData[pos + i * nOrthoXSize] = (0);
					pFloatData[pos + i * nOrthoXSize + nOrthoXSize - 1] = (0);
				}
			}
		}

		if (i >= (nPrec * onePerc))
		{
			nPrec++;
			if (nPrec % 10 == 0) {
				char lineInfo[32];
				sprintf_s(lineInfo, "%d%%", nPrec);
				std::cout << lineInfo;
				fflush(stdout);
			}
			else if (nPrec % 2 == 0) {
				std::cout << ".";
				fflush(stdout);
			}
		}
	}

	delete[]g;

	bool dark = 1;
	int up = 0, down = 0, left = 0, right = 0;

	for (int i = 0; i < nOrthoYSize; i++) {
		for (int j = 0; j < nOrthoXSize; j++) {
			if (nDataType == GDT_Byte) {
				unsigned char* pByteData = (unsigned char*)pOrthoData;
				for (int k = 0; k < nBands; ++k) {
					size_t pos = k * nOrthoYSize * nOrthoXSize;
					if (pByteData[pos + i * nOrthoXSize + j])
					{
						dark = 0;
						break;
					}
				}
			}
			else if (nDataType == GDT_UInt16) {
				unsigned short* pUShortData = (unsigned short*)pOrthoData;
				for (int k = 0; k < nBands; ++k) {
					size_t pos = k * nOrthoYSize * nOrthoXSize;
					if (pUShortData[pos + i * nOrthoXSize + j])
					{
						dark = 0;
						break;
					}
				}
			}
			else if (nDataType == GDT_Int16) {
				short* pShortData = (short*)pOrthoData;
				for (int k = 0; k < nBands; ++k) {
					size_t pos = k * nOrthoYSize * nOrthoXSize;
					if (pShortData[pos + i * nOrthoXSize + j])
					{
						dark = 0;
						break;
					}
				}
			}
			else if (nDataType == GDT_Float32) {
				float* pFloatData = (float*)pOrthoData;
				for (int k = 0; k < nBands; ++k) {
					size_t pos = k * nOrthoYSize * nOrthoXSize;
					if (pFloatData[pos + i * nOrthoXSize + j])
					{
						dark = 0;
						break;
					}
				}
			}
			if (!dark) { break; }
		}
		if (!dark) { break; }
		up++;
	}
	std::cout << std::endl << "Dark line for Up side done: " << up << std::endl;

	dark = 1;
	for (int i = nOrthoYSize - 1; i > -1; i--) {
		for (int j = 0; j < nOrthoXSize; j++) {
			if (nDataType == GDT_Byte) {
				unsigned char* pByteData = (unsigned char*)pOrthoData;
				for (int k = 0; k < nBands; ++k) {
					size_t pos = k * nOrthoYSize * nOrthoXSize;
					if (pByteData[pos + i * nOrthoXSize + j])
					{
						dark = 0;
						break;
					}
				}
			}
			else if (nDataType == GDT_UInt16) {
				unsigned short* pUShortData = (unsigned short*)pOrthoData;
				for (int k = 0; k < nBands; ++k) {
					size_t pos = k * nOrthoYSize * nOrthoXSize;
					if (pUShortData[pos + i * nOrthoXSize + j])
					{
						dark = 0;
						break;
					}
				}
			}
			else if (nDataType == GDT_Int16) {
				short* pShortData = (short*)pOrthoData;
				for (int k = 0; k < nBands; ++k) {
					size_t pos = k * nOrthoYSize * nOrthoXSize;
					if (pShortData[pos + i * nOrthoXSize + j])
					{
						dark = 0;
						break;
					}
				}
			}
			else if (nDataType == GDT_Float32) {
				float* pFloatData = (float*)pOrthoData;
				for (int k = 0; k < nBands; ++k) {
					size_t pos = k * nOrthoYSize * nOrthoXSize;
					if (pFloatData[pos + i * nOrthoXSize + j])
					{
						dark = 0;
						break;
					}
				}
			}
			if (!dark) { break; }
		}
		if (!dark) { break; }
		down++;
	}
	std::cout << "Dark line for Down side done: " << down << std::endl;

	dark = 1;
	for (int j = 0; j < nOrthoXSize; j++) {
		for (int i = 0; i < nOrthoYSize; i++) {
			if (nDataType == GDT_Byte) {
				unsigned char* pByteData = (unsigned char*)pOrthoData;
				for (int k = 0; k < nBands; ++k) {
					size_t pos = k * nOrthoYSize * nOrthoXSize;
					if (pByteData[pos + i * nOrthoXSize + j])
					{
						dark = 0;
						break;
					}
				}
			}
			else if (nDataType == GDT_UInt16) {
				unsigned short* pUShortData = (unsigned short*)pOrthoData;
				for (int k = 0; k < nBands; ++k) {
					size_t pos = k * nOrthoYSize * nOrthoXSize;
					if (pUShortData[pos + i * nOrthoXSize + j])
					{
						dark = 0;
						break;
					}
				}
			}
			else if (nDataType == GDT_Int16) {
				short* pShortData = (short*)pOrthoData;
				for (int k = 0; k < nBands; ++k) {
					size_t pos = k * nOrthoYSize * nOrthoXSize;
					if (pShortData[pos + i * nOrthoXSize + j])
					{
						dark = 0;
						break;
					}
				}
			}
			else if (nDataType == GDT_Float32) {
				float* pFloatData = (float*)pOrthoData;
				for (int k = 0; k < nBands; ++k) {
					size_t pos = k * nOrthoYSize * nOrthoXSize;
					if (pFloatData[pos + i * nOrthoXSize + j])
					{
						dark = 0;
						break;
					}
				}
			}
			if (!dark) { break; }
		}
		if (!dark) { break; }
		left++;
	}
	std::cout << "Dark line for Left side done: " << left << std::endl;

	dark = 1;
	for (int j = nOrthoXSize - 1; j > -1; j--) {
		for (int i = 0; i < nOrthoYSize; i++) {
			if (nDataType == GDT_Byte) {
				unsigned char* pByteData = (unsigned char*)pOrthoData;
				for (int k = 0; k < nBands; ++k) {
					size_t pos = k * nOrthoYSize * nOrthoXSize;
					if (pByteData[pos + i * nOrthoXSize + j])
					{
						dark = 0;
						break;
					}
				}
			}
			else if (nDataType == GDT_UInt16) {
				unsigned short* pUShortData = (unsigned short*)pOrthoData;
				for (int k = 0; k < nBands; ++k) {
					size_t pos = k * nOrthoYSize * nOrthoXSize;
					if (pUShortData[pos + i * nOrthoXSize + j])
					{
						dark = 0;
						break;
					}
				}
			}
			else if (nDataType == GDT_Int16) {
				short* pShortData = (short*)pOrthoData;
				for (int k = 0; k < nBands; ++k) {
					size_t pos = k * nOrthoYSize * nOrthoXSize;
					if (pShortData[pos + i * nOrthoXSize + j])
					{
						dark = 0;
						break;
					}
				}
			}
			else if (nDataType == GDT_Float32) {
				float* pFloatData = (float*)pOrthoData;
				for (int k = 0; k < nBands; ++k) {
					size_t pos = k * nOrthoYSize * nOrthoXSize;
					if (pFloatData[pos + i * nOrthoXSize + j])
					{
						dark = 0;
						break;
					}
				}
			}
			if (!dark) { break; }
		}
		if (!dark) { break; }
		right++;
	}
	std::cout << "Dark line for Right side done: " << right << std::endl;

	int FinalXSize = nOrthoXSize - left - right;
	int FinalYSize = nOrthoYSize - up - down;

	GDALDataset* poDstDS;
	poDstDS = poDriver->Create(outFilename.c_str(), FinalXSize, FinalYSize, nBands, nDataType, NULL);
	if (poDstDS == NULL) {
		GDALClose(poSrcDS);
		CPLFree(pSrcData);
		return false;
	}
	else {
		std::cout << "Output image prepared: " << outFilename << std::endl;
	}

	double adfGeoTransform[6];
	adfGeoTransform[0] = orthoXmin + orthLogicHexX * (left - 1  /*+ 2*/);
	adfGeoTransform[1] = orthLogicHexX;
	adfGeoTransform[2] = 0;
	adfGeoTransform[3] = orthoYmax - orthLogicHexY * (up + 1  /*+ 20*/);
	adfGeoTransform[4] = 0;
	adfGeoTransform[5] = -orthLogicHexY;
	poDstDS->SetGeoTransform(adfGeoTransform);

	void* pOutData;
	pOutData = (void*)VSIMalloc2(FinalXSize, FinalYSize * nBands * nDepth);
	if (pOutData == NULL) {
		GDALClose(poDstDS);
		VSIFree(pSrcData);
		return false;
	}

	for (int i = 0; i < FinalYSize; i++) {
		for (int j = 0; j < FinalXSize; j++) {
			if (nDataType == GDT_Byte) {
				unsigned char* pByteData = (unsigned char*)pOrthoData;
				unsigned char* pFiByteData = (unsigned char*)pOutData;
				for (int k = 0; k < nBands; ++k) {
					size_t posF = k * FinalYSize * FinalXSize;
					size_t pos = k * nOrthoYSize * nOrthoXSize;
					pFiByteData[posF + i * FinalXSize + j] = pByteData[pos + (i + up) * nOrthoXSize + j + left];
				}
			}
			else if (nDataType == GDT_UInt16) {
				unsigned short* pFiUShorData = (unsigned short*)pOutData;
				unsigned short* pUShortData = (unsigned short*)pOrthoData;
				for (int k = 0; k < nBands; ++k) {
					size_t posF = k * FinalYSize * FinalXSize;
					size_t pos = k * nOrthoYSize * nOrthoXSize;
					pFiUShorData[posF + i * FinalXSize + j] = pUShortData[pos + (i + up) * nOrthoXSize + j + left];
				}
			}
			else if (nDataType == GDT_Int16) {
				short* pFiShortData = (short*)pOutData;
				short* pShortData = (short*)pOrthoData;
				for (int k = 0; k < nBands; ++k) {
					size_t posF = k * FinalYSize * FinalXSize;
					size_t pos = k * nOrthoYSize * nOrthoXSize;
					pFiShortData[posF + i * FinalXSize + j] = pShortData[pos + (i + up) * nOrthoXSize + j + left];
				}
			}
			else if (nDataType == GDT_Float32) {
				float* pFiFloatData = (float*)pOutData;
				float* pFloatData = (float*)pOrthoData;
				for (int k = 0; k < nBands; ++k) {
					size_t posF = k * FinalYSize * FinalXSize;
					size_t pos = k * nOrthoYSize * nOrthoXSize;
					pFiFloatData[posF + i * FinalXSize + j] = pFloatData[pos + (i + up) * nOrthoXSize + j + left];
				}
			}
		}
	}

	CPLErr err = poDstDS->RasterIO(GF_Write,
		0, 0,
		FinalXSize, FinalYSize,
		(void*)pOutData,
		FinalXSize, FinalYSize,
		nDataType, nBands, NULL, 0, 0, 0);

	CPLFree(pSrcData);
	VSIFree(pOrthoData);
	VSIFree(pOutData);
	GDALClose(poDstDS);

	return 1;
}