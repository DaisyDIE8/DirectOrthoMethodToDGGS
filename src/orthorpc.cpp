#include <iostream>
#include <fstream>
#include <sstream>

#include "gdal.h"
#include "gdal_priv.h"
#include "cpl_conv.h"

#include <boost\filesystem.hpp>

#include "georaster.h"
#include "tifrpcdataset.h"
#include "rpcsensormodel.h"
#include "RtA4HexDGGS.h"
#include "LL2LogicHexSD.h"
#include "rpcestimator.h"
#include "experiment.h"

#include <omp.h>

bool  InterpolateZ(const double& EX, const double& NY, double* HZ, const DEM_GRD& deminfo, float* dem)
{
	int II, JJ, i, j;
	double dx, dy, z1, z2, z3, z4;

	II = static_cast<int>((NY - deminfo.Ymax) / deminfo.gridY);
	JJ = static_cast<int>((EX - deminfo.Xmin) / deminfo.gridX);

	if (II >= deminfo.row - 1 || II < 0 || JJ >= deminfo.col - 1 || JJ < 0)
	{
		if (II >= deminfo.row - 1)
			i = deminfo.row - 1;
		else if (II < 0)
			i = 0;
		else
			i = II;

		if (JJ >= deminfo.col - 1)
			j = deminfo.col - 1;
		else if (JJ < 0)
			j = 0;
		else
			j = JJ;

		*HZ = dem[i * deminfo.col + j];
		return false;
	}

	dx = (EX - deminfo.Xmin - JJ * deminfo.gridX) / deminfo.gridX;
	dy = (NY - deminfo.Ymax - II * deminfo.gridY) / deminfo.gridY;

	z1 = dem[II * deminfo.col + JJ];
	z2 = dem[II * deminfo.col + JJ + 1];
	z3 = dem[(II + 1) * deminfo.col + JJ + 1];
	z4 = dem[(II + 1) * deminfo.col + JJ];

	*HZ = (1 - dx) * (1 - dy) * z1 + dx * (1 - dy) * z2 + dx * dy * z3 + (1 - dx) * dy * z4;
	return true;
}

char* SanitizeSRS(const char* pszUserInput)
{
	OGRSpatialReferenceH hSRS;
	char* pszResult = NULL;

	CPLErrorReset();

	hSRS = OSRNewSpatialReference(NULL);
	if (OSRSetFromUserInput(hSRS, pszUserInput) == OGRERR_NONE)
		OSRExportToWkt(hSRS, &pszResult);
	else
	{
		CPLError(CE_Failure, CPLE_AppDefined,
			"Translating source or target SRS failed:\n%s",
			pszUserInput);
		exit(EXIT_FAILURE);
	}

	OSRDestroySpatialReference(hSRS);

	return pszResult;
}

bool StringContains(const std::string& str, const std::string& sub_str) {
	return str.find(sub_str) != std::string::npos;
}

std::string StringReplace(const std::string& str, const std::string& old_str,
	const std::string& new_str) {
	if (old_str.empty()) {
		return str;
	}
	size_t position = 0;
	std::string mod_str = str;
	while ((position = mod_str.find(old_str, position)) != std::string::npos) {
		mod_str.replace(position, old_str.size(), new_str);
		position += new_str.size();
	}
	return mod_str;
}

std::string StringGetFileName(const std::string& path) {
	std::string modPath = StringReplace(path, "\\", "/");
	std::string::size_type startpos = modPath.rfind('/') == std::string::npos ? modPath.length() : modPath.rfind('/') + 1;
	std::string::size_type endpos = modPath.rfind('.') == std::string::npos ? modPath.length() : modPath.rfind('.');
	return modPath.substr(startpos, endpos - startpos);
}

bool  orthoRectificationTra(
	const std::shared_ptr<TifRPCDataSet> tifDS,
	const std::string srsWkt,
	CGeoRaster* demDS,
	const std::string orthofilename,
	double gsd)
{
	std::string srcImageFile = tifDS->imageFileName();
	GDALDataset* poSrcDS = (GDALDataset*)GDALOpen(srcImageFile.c_str(), GA_ReadOnly);
	if (poSrcDS == NULL) {
		std::cout << "Source image can not be opened:" << srcImageFile << std::endl;
		return false;
	}

	int nWidth = poSrcDS->GetRasterXSize();
	int nHeight = poSrcDS->GetRasterYSize();
	int nBands = poSrcDS->GetRasterCount();
	GDALDataType nDataType = poSrcDS->GetRasterBand(1)->GetRasterDataType();
	int nDepth = GDALGetDataTypeSize((GDALDataType)nDataType) / 8;

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

	GDALClose(poSrcDS);

	GDALDriver* poDriver;
	poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
	if (poDriver == NULL) {
		CPLFree(pSrcData);
		return false;
	}

	double orthoGSDX = tifDS->level0GSD();
	double orthoGSDY = orthoGSDX;
	if (gsd > std::numeric_limits<float>::epsilon())
		orthoGSDX = orthoGSDY = gsd;

	double minLon = std::numeric_limits<double>::max();
	double minLat = std::numeric_limits<double>::max();
	double maxLon = -std::numeric_limits<double>::max();
	double maxLat = -std::numeric_limits<double>::max();
	for (auto i = 0; i < 4; i++) {
		ground_point_struct gp = tifDS->footPrints()[i];
		minLon = std::min(gp.x, minLon);
		maxLon = std::max(gp.x, maxLon);

		minLat = std::min(gp.y, minLat);
		maxLat = std::max(gp.y, maxLat);
	}

	int nOrthoXSize = static_cast<int>((maxLon - minLon) / orthoGSDX + 0.5);
	if (nOrthoXSize % 4 != 0)
	{
		nOrthoXSize = 4 * (nOrthoXSize / 4 + 1);
		maxLon = minLon + nOrthoXSize * orthoGSDX;
	}

	int nOrthoYSize = static_cast<int>((maxLat - minLat) / orthoGSDY + 0.5);

	OGRSpatialReference demSRS = demDS->referenceCS();

	OGRCoordinateTransformation* poTransEo2Dem = nullptr;
	OGRCoordinateTransformation* poTransLatLon2Dem = nullptr;
	OGRSpatialReference eoSRS;
	if (GDAL_VERSION_MAJOR >= 3)
		eoSRS.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
	eoSRS.importFromWkt(srsWkt.c_str());

	poTransEo2Dem = OGRCreateCoordinateTransformation(&eoSRS, &demSRS);
	if (poTransEo2Dem == nullptr) {
		std::cerr << "Failed to create transformation from rpc/rpb to DEM." << std::endl;
		CPLFree(pSrcData);
		return false;
	}

	double uldemx = minLon;
	double uldemy = maxLat;

	poTransEo2Dem->Transform(1, &uldemx, &uldemy);

	double lrdemx = maxLon;
	double lrdemy = minLat;

	poTransEo2Dem->Transform(1, &lrdemx, &lrdemy);

	double demMinX = std::min(uldemx, lrdemx);
	double demMaxX = std::max(uldemx, lrdemx);

	double demMinY = std::min(uldemy, lrdemy);
	double demMaxY = std::max(uldemy, lrdemy);

	double row, col;
	demDS->ground2Image(demMinX, demMaxY, col, row);
	int ulDemAOIRow = static_cast<int>(row - 3) >= 0 ? static_cast<int>(row - 3) : 0;
	int ulDemAOICol = static_cast<int>(col - 3) >= 0 ? static_cast<int>(col - 3) : 0;

	demDS->ground2Image(demMaxX, demMinY, col, row);
	int lrDemAOIRow = static_cast<int>(row + 3) > (demDS->height() - 1) ? (demDS->height() - 1) : static_cast<int>(row + 3);
	int lrDemAOICol = static_cast<int>(col + 3) > (demDS->width() - 1) ? (demDS->width() - 1) : static_cast<int>(col + 3);

	DEM_GRD demAOI;
	demDS->gsd(demAOI.gridX, demAOI.gridY);
	demDS->image2Ground(ulDemAOICol, ulDemAOIRow, demAOI.Xmin, demAOI.Ymax);
	demDS->image2Ground(lrDemAOICol, lrDemAOIRow, demAOI.Xmax, demAOI.Ymin);
	demAOI.row = (lrDemAOIRow - ulDemAOIRow + 1);
	demAOI.col = (lrDemAOICol - ulDemAOICol + 1);

	float* demDataAOI = (float*)VSIMalloc2(demAOI.row, demAOI.col * sizeof(float));
	if (demDataAOI == nullptr) {
		GDALClose(poSrcDS);
		CPLFree(pSrcData);
		return false;
	}

	if (!demDS->readBlock(ulDemAOIRow, ulDemAOICol, demAOI.col, demAOI.row, (void*)demDataAOI, GDT_Float32)) {
		std::cerr << "Failed to load the elevation data in dem aoi: " << demDS->imageFileName() << std::endl;
		GDALClose(poSrcDS);
		CPLFree(pSrcData);
		CPLFree(demDataAOI);
		return false;
	}
	else
		std::cout << "DEM aoi loaded. " << std::endl;

	boost::filesystem::path p1(orthofilename);
	if (boost::filesystem::exists(orthofilename))
		boost::filesystem::remove(p1);

	std::string ovrFile = orthofilename + ".ovr";
	boost::filesystem::path p2(ovrFile);
	if (boost::filesystem::exists(ovrFile))
		boost::filesystem::remove(p2);

	GDALDataset* poDstDS;
	poDstDS = poDriver->Create(orthofilename.c_str(), nOrthoXSize, nOrthoYSize, nBands, nDataType, NULL);
	if (poDstDS == NULL) {
		GDALClose(poSrcDS);
		CPLFree(pSrcData);
		return false;
	}

	if (!StringContains(srsWkt, "TSEA")) {
		poDstDS->SetProjection(srsWkt.c_str());
	}

	double adfGeoTransform[6];
	adfGeoTransform[0] = minLon;
	adfGeoTransform[1] = orthoGSDX;
	adfGeoTransform[2] = 0;
	adfGeoTransform[3] = maxLat;
	adfGeoTransform[4] = 0;
	adfGeoTransform[5] = -orthoGSDY;
	poDstDS->SetGeoTransform(adfGeoTransform);

	void* pOrthoData;
	pOrthoData = (void*)VSIMalloc2(nOrthoXSize, nOrthoYSize * nBands * nDepth);
	if (pOrthoData == NULL) {
		GDALClose(poDstDS);
		VSIFree(pSrcData);
		return false;
	}

	double avgh = tifDS->meanTerrainHeight();
	RPCSensorModel* sm = dynamic_cast<RPCSensorModel*>(tifDS->sensorModel());

	std::cout << "Start to ortho rectification..." << std::endl;

	int onePerc = static_cast<int>(nOrthoYSize * 0.01), nPrec = 0;
	float* g = new float[nBands];

#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < nOrthoYSize; i++)
	{
		double northing = maxLat - orthoGSDY * i;

		std::vector<float>aLine;
		aLine.resize(nOrthoXSize * nBands);
		for (int j = 0; j < nOrthoXSize; j++)
		{
			double easting = minLon + j * orthoGSDX;

			double demx = easting;
			double demy = northing;
			poTransEo2Dem->Transform(1, &demx, &demy);

			double hz = avgh;
			if (!InterpolateZ(demx, demy, &hz, demAOI, demDataAOI))
				hz = avgh;

			double samp = 0, line = 0;
			sm->Conv_Map_To_SL(easting, northing, hz, &samp, &line);

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
					pByteData[pos + i * nOrthoXSize + j] = static_cast<unsigned char>(g[k]);
				}
			}
			else if (nDataType == GDT_UInt16) {
				unsigned short* pUShortData = (unsigned short*)pOrthoData;
				for (int k = 0; k < nBands; ++k)
				{
					size_t pos = k * nOrthoYSize * nOrthoXSize;
					pUShortData[pos + i * nOrthoXSize + j] = static_cast<unsigned short>(g[k]);
				}
			}
			else if (nDataType == GDT_Int16) {
				short* pShortData = (short*)pOrthoData;
				for (int k = 0; k < nBands; ++k)
				{
					size_t pos = k * nOrthoYSize * nOrthoXSize;
					pShortData[pos + i * nOrthoXSize + j] = static_cast<short>(g[k]);
				}
			}
			else if (nDataType == GDT_Float32) {
				float* pFloatData = (float*)pOrthoData;
				for (int k = 0; k < nBands; ++k)
				{
					size_t pos = k * nOrthoYSize * nOrthoXSize;
					pFloatData[pos + i * nOrthoXSize + j] = (g[k]);
				}
			}
		}

		if (i >= (nPrec * onePerc)) {
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

	CPLErr err = poDstDS->RasterIO(GF_Write,
		0, 0,
		nOrthoXSize, nOrthoYSize,
		(void*)pOrthoData,
		nOrthoXSize, nOrthoYSize,
		nDataType, nBands, NULL, 0, 0, 0);

	if (err == CE_Failure)
		std::cout << "Write ortho image file failed: " << orthofilename << std::endl;

	if (poTransLatLon2Dem)
		OGRCoordinateTransformation::DestroyCT(poTransLatLon2Dem);

	if (poTransEo2Dem)
		OGRCoordinateTransformation::DestroyCT(poTransEo2Dem);

	CPLFree(pSrcData);
	VSIFree(pOrthoData);
	CPLFree(demDataAOI);
	GDALClose(poDstDS);

	return true;
}

bool  orthoRectification(
	const std::shared_ptr<TifRPCDataSet> tifDS,
	const std::string srsWkt,
	CGeoRaster* demDS,
	const std::string orthofilename,
	double gsd,
	InterpolateType t)
{
	std::string srcImageFile = tifDS->imageFileName();
	GDALDataset* poSrcDS = (GDALDataset*)GDALOpen(srcImageFile.c_str(), GA_ReadOnly);
	if (poSrcDS == NULL) {
		std::cout << "Source image can not be opened:" << srcImageFile << std::endl;
		return false;
	}

	int nWidth = poSrcDS->GetRasterXSize();
	int nHeight = poSrcDS->GetRasterYSize();
	int nBands = poSrcDS->GetRasterCount();
	GDALDataType nDataType = poSrcDS->GetRasterBand(1)->GetRasterDataType();
	int nDepth = GDALGetDataTypeSize((GDALDataType)nDataType) / 8;

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

	GDALClose(poSrcDS);

	GDALDriver* poDriver;
	poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
	if (poDriver == NULL) {
		CPLFree(pSrcData);
		return false;
	}

	const int level_ = 18;
	RtSphHexDGGS imaDGGS = RtSphHexDGGS(6371, level_);

	unsigned short  rtFaceId = -1;
	std::string tesa;
	std::stringstream ss(srsWkt);
	ss >> tesa >> rtFaceId;
	assert(rtFaceId >= 0 && rtFaceId < 60);

	double orthLogicHexX = imaDGGS.lengthK_ / 2;
	double orthLogicHexY = imaDGGS.lengthH_ / 2;

	double orthoXmin = std::numeric_limits<double>::max();
	double orthoYmin = std::numeric_limits<double>::max();
	double orthoXmax = -std::numeric_limits<double>::max();
	double orthoYmax = -std::numeric_limits<double>::max();
	for (auto i = 0; i < 4; i++) {
		ground_point_struct gp = tifDS->footPrints()[i];
		orthoXmin = std::min(gp.x * 0.001, orthoXmin);
		orthoXmax = std::max(gp.x * 0.001, orthoXmax);

		orthoYmin = std::min(gp.y * 0.001, orthoYmin);
		orthoYmax = std::max(gp.y * 0.001, orthoYmax);
	}

	DGridCoord2Dij first = imaDGGS.GetIJfromTricoord(DGridCoord2D(rtFaceId, orthoXmin, orthoYmax));
	orthoXmin = imaDGGS.GetTricoordfromIJ(first).pt.x - orthLogicHexX;
	orthoYmax = imaDGGS.GetTricoordfromIJ(first).pt.y;

	DGridCoord2Dij last = imaDGGS.GetIJfromTricoord(DGridCoord2D(rtFaceId, orthoXmax, orthoYmin));
	orthoXmax = imaDGGS.GetTricoordfromIJ(last).pt.x + orthLogicHexX;
	orthoYmin = imaDGGS.GetTricoordfromIJ(last).pt.y;

	int nOrthoXSize = static_cast<int>((orthoXmax - orthoXmin) / orthLogicHexX + 0.5);
	int nOrthoYSize = static_cast<int>((orthoYmax - orthoYmin) / orthLogicHexY + 0.5);

	OGRCoordinateTransformation* poTransEo2Dem = nullptr;
	OGRCoordinateTransformation* poTransLatLon2Dem = nullptr;
	OGRSpatialReference demSRS = demDS->referenceCS();
	OGRSpatialReference oLatLon;
	if (GDAL_VERSION_MAJOR >= 3)
		oLatLon.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);

	oLatLon.importFromEPSG(4326);

	poTransLatLon2Dem = OGRCreateCoordinateTransformation(&oLatLon, &demSRS);
	if (poTransLatLon2Dem == nullptr) {
		std::cerr << "Failed to crerate transformation from WGS84 to DEM." << std::endl;
		CPLFree(pSrcData);
		return false;
	}

	orthoXmin *= 1000;
	orthoYmin *= 1000;
	orthoXmax *= 1000;
	orthoYmax *= 1000;
	orthLogicHexY *= 1000;
	orthLogicHexX *= 1000;

	DGridCoord2D ixy(rtFaceId, orthoXmin * 0.001, orthoYmax * 0.001);
	double uldemx = imaDGGS.RTInvSlice(ixy).lon * RAD_TO_DEG;
	double uldemy = imaDGGS.RTInvSlice(ixy).lat * RAD_TO_DEG;
	poTransLatLon2Dem->Transform(1, &uldemx, &uldemy);

	DGridCoord2D lxy(rtFaceId, orthoXmax * 0.001, orthoYmin * 0.001);
	double lrdemx = imaDGGS.RTInvSlice(lxy).lon * RAD_TO_DEG;
	double lrdemy = imaDGGS.RTInvSlice(lxy).lat * RAD_TO_DEG;
	poTransLatLon2Dem->Transform(1, &lrdemx, &lrdemy);

	double demMinX = std::min(uldemx, lrdemx);
	double demMaxX = std::max(uldemx, lrdemx);

	double demMinY = std::min(uldemy, lrdemy);
	double demMaxY = std::max(uldemy, lrdemy);

	double row, col;
	demDS->ground2Image(demMinX, demMaxY, col, row);
	int ulDemAOIRow = static_cast<int>(row - 3) >= 0 ? static_cast<int>(row - 3) : 0;
	int ulDemAOICol = static_cast<int>(col - 3) >= 0 ? static_cast<int>(col - 3) : 0;

	demDS->ground2Image(demMaxX, demMinY, col, row);
	int lrDemAOIRow = static_cast<int>(row + 3) > (demDS->height() - 1) ? (demDS->height() - 1) : static_cast<int>(row + 3);
	int lrDemAOICol = static_cast<int>(col + 3) > (demDS->width() - 1) ? (demDS->width() - 1) : static_cast<int>(col + 3);

	DEM_GRD demAOI;
	demDS->gsd(demAOI.gridX, demAOI.gridY);
	demDS->image2Ground(ulDemAOICol, ulDemAOIRow, demAOI.Xmin, demAOI.Ymax);
	demDS->image2Ground(lrDemAOICol, lrDemAOIRow, demAOI.Xmax, demAOI.Ymin);
	demAOI.row = (lrDemAOIRow - ulDemAOIRow + 1);
	demAOI.col = (lrDemAOICol - ulDemAOICol + 1);

	float* demDataAOI = (float*)VSIMalloc2(demAOI.row, demAOI.col * sizeof(float));
	if (demDataAOI == nullptr) {
		GDALClose(poSrcDS);
		CPLFree(pSrcData);
		return false;
	}

	if (!demDS->readBlock(ulDemAOIRow, ulDemAOICol, demAOI.col, demAOI.row, (void*)demDataAOI, GDT_Float32)) {
		std::cerr << "Failed to load the elevation data in dem aoi: " << demDS->imageFileName() << std::endl;
		GDALClose(poSrcDS);
		CPLFree(pSrcData);
		CPLFree(demDataAOI);
		return false;
	}
	else
		std::cout << "DEM aoi loaded. " << std::endl;

	boost::filesystem::path p1(orthofilename);
	if (boost::filesystem::exists(orthofilename))
		boost::filesystem::remove(p1);

	std::string ovrFile = orthofilename + ".ovr";
	boost::filesystem::path p2(ovrFile);
	if (boost::filesystem::exists(ovrFile))
		boost::filesystem::remove(p2);

	GDALDataset* poDstDS;
	poDstDS = poDriver->Create(orthofilename.c_str(), nOrthoXSize, nOrthoYSize, nBands, nDataType, NULL);
	if (poDstDS == NULL) {
		GDALClose(poSrcDS);
		CPLFree(pSrcData);
		return false;
	}

	if (!StringContains(srsWkt, "TSEA")) {
		poDstDS->SetProjection(srsWkt.c_str());
	}

	double adfGeoTransform[6];
	adfGeoTransform[0] = orthoXmin;
	adfGeoTransform[1] = orthLogicHexX;
	adfGeoTransform[2] = 0;
	adfGeoTransform[3] = orthoYmax;
	adfGeoTransform[4] = 0;
	adfGeoTransform[5] = -orthLogicHexY;
	poDstDS->SetGeoTransform(adfGeoTransform);

	void* pOrthoData;
	pOrthoData = (void*)VSIMalloc2(nOrthoXSize, nOrthoYSize * nBands * nDepth);
	if (pOrthoData == NULL) {
		GDALClose(poDstDS);
		VSIFree(pSrcData);
		return false;
	}

	double avgh = tifDS->meanTerrainHeight();
	RPCSensorModel* sm = dynamic_cast<RPCSensorModel*>(tifDS->sensorModel());

	std::cout << "Start to ortho rectification..." << std::endl;

	int onePerc = static_cast<int>(nOrthoYSize * 0.01), nPrec = 0;
	float* g = new float[nBands];

#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < nOrthoYSize; i++)
	{
		double northing = orthoYmax - orthLogicHexY * i;

		std::vector<float>aLine;
		aLine.resize(nOrthoXSize * nBands);

		for (int j = 0; j < nOrthoXSize; j += 2)
		{
			int m = j;
			if (i % 2) {
				m = j + 1;
				if (j == nOrthoXSize - 2) { continue; }
			}

			double easting = orthoXmin + (m + 1) * orthLogicHexX;

			DGridCoord2D ccxy(rtFaceId, easting * 0.001, northing * 0.001);
			GeoCoord gpt = imaDGGS.RTInvSlice(ccxy);
			double demx = gpt.lon * RAD_TO_DEG;
			double demy = gpt.lat * RAD_TO_DEG;
			poTransLatLon2Dem->Transform(1, &demx, &demy);

			double hz = avgh;
			if (!InterpolateZ(demx, demy, &hz, demAOI, demDataAOI))
				hz = avgh;

			double samp = 0, line = 0;
			sm->Conv_Map_To_SL(easting, northing, hz, &samp, &line);

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

		if (i >= (nPrec * onePerc)) {
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

	CPLErr err = poDstDS->RasterIO(GF_Write,
		0, 0,
		nOrthoXSize, nOrthoYSize,
		(void*)pOrthoData,
		nOrthoXSize, nOrthoYSize,
		nDataType, nBands, NULL, 0, 0, 0);

	if (err == CE_Failure)
		std::cout << "Write ortho image file failed: " << orthofilename << std::endl;

	if (poTransLatLon2Dem)
		OGRCoordinateTransformation::DestroyCT(poTransLatLon2Dem);

	if (poTransEo2Dem)
		OGRCoordinateTransformation::DestroyCT(poTransEo2Dem);

	CPLFree(pSrcData);
	VSIFree(pOrthoData);
	CPLFree(demDataAOI);
	GDALClose(poDstDS);

	return true;
}

bool orthoRPC(bool ifDGGS, std::string inputImage, std::string rpbFile, std::string demFile, std::string prjFile, std::string outFile, double gsd, InterpolateType t)
{
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");

	std::shared_ptr<TifRPCDataSet> tifRpcDS = std::make_shared<TifRPCDataSet>();
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

	if (!StringContains(srsWKT, "TSEA")) {
		char* pszTargetSRS = NULL;
		pszTargetSRS = SanitizeSRS(srsWKT.c_str());
		if (pszTargetSRS) {
			if (pszTargetSRS[0] != '\0' && eoSRS.importFromWkt(pszTargetSRS) != OGRERR_NONE) {
				CPLError(CE_Failure, CPLE_AppDefined, "Failed to import coordinate system `%s'.", pszTargetSRS);
				return EXIT_FAILURE;
			}
		}
	}

	std::shared_ptr<CGeoRaster> demDS = std::make_shared<CGeoRaster>();
	if (!demDS->readMetaData(demFile)) {
		std::cerr << "Read dem metadata failed: " << demFile << std::endl;
		return EXIT_FAILURE;
	}

	if (ifDGGS) {
		if (orthoRectification(tifRpcDS, srsWKT, demDS.get(), outFile, gsd, t)) {
			std::cout << std::endl << "Ortho image generated: " << outFile << std::endl;

			std::string orthoProj = outFile + ".prj";

			std::ofstream os;
			os.open(orthoProj);
			if (os.fail()) {
				std::cerr << "Can not write projection to file: " << orthoProj << std::endl << std::endl;
				return EXIT_FAILURE;
			}
			else {
				os << srsWKT << std::endl;
				os.close();
			}
		}
	}
	else {
		if (orthoRectificationTra(tifRpcDS, srsWKT, demDS.get(), outFile, gsd)) {
			std::cout << std::endl << "Ortho image generated: " << outFile << std::endl;

			std::string orthoProj = outFile + ".prj";

			std::ofstream os;
			os.open(orthoProj);
			if (os.fail()) {
				std::cerr << "Can not write projection to file: " << orthoProj << std::endl << std::endl;
				return EXIT_FAILURE;
			}
			else {
				os << srsWKT << std::endl;
				os.close();
			}
		}
	}

	return true;
}

/*
 * @brief	 Orthorectification of hexagonal images using the proposed method.
 *
 * @param inputImage: The path of the original image.
 * @param rpbexFile: The path of the extended rpbex file.
 * @param demFile: The path of the digital elevation model (DEM) file.
 * @param prjexFile: The path of the extended projection file.
 * @param proposeOut: The path of the output orthophotograph file.
 *
 * @return	 Runtime(s).
 */
double proposedMethod(std::string inputImage, std::string rpbexFile, std::string demFile, std::string prjexFile, std::string proposeOut, InterpolateType type)
{
	double gsd = 16;

	clock_t startT = clock();
	if (!orthoRPC(true, inputImage, rpbexFile, demFile, prjexFile, proposeOut, gsd, type)) {
		std::cout << "ortho Wrong!" << std::endl;
		std::exit(EXIT_SUCCESS);
	}
	clock_t endT = clock();
	double totalSeconds = static_cast<double>((endT - startT) * 1.0f / CLOCKS_PER_SEC);
	std::cout << std::endl << "proposedMethod finished, total time used:" << std::setprecision(5) << totalSeconds << " seconds." << std::endl << std::endl;

	return totalSeconds;
}

/*
 * @brief	 Orthorectification of hexagonal images using the traditional method.
 *
 * @param inputImage: The path of the original image.
 * @param rpbFile: The path of the original rpb file.
 * @param demFile: The path of the digital elevation model (DEM) file.
 * @param prjexFile: The path of the extended projection file.
 * @param traditionalOut: The path of the output orthophotograph file.
 *
 * @return	 Runtime(s).
 */
double traditionalMethod(std::string inputImage, std::string rpbFile, std::string demFile, std::string prjexFile, std::string traditionOut, InterpolateType type)
{
	double gsd = 1.5e-04;
	std::string LLimage = inputImage;
	LLimage.insert(61, "_LL");
	LLimage.pop_back();

	std::string prjFile =
		".\\data\\RPC\\normal.prj";

	clock_t startT = clock();
	if (!orthoRPC(false, inputImage, rpbFile, demFile, prjFile, LLimage, gsd, type)) {
		std::cout << "ortho Wrong!" << std::endl;
		std::exit(EXIT_SUCCESS);
	}

	if (!ll2logicHexSD(LLimage, traditionOut, prjexFile, type)) {
		std::cout << "reprojection Wrong!" << std::endl;
		std::exit(EXIT_SUCCESS);
	}
	clock_t endT = clock();
	double totalSeconds = static_cast<double>((endT - startT) * 1.0f / CLOCKS_PER_SEC);
	std::cout << std::endl << "traditionalMethod finished, total time used:" << std::setprecision(5) << totalSeconds << " seconds." << std::endl << std::endl;

	return totalSeconds;
}

int main(int argc, char* argv[])
{
	std::string inputImage, rpbFile, rpbexFile, demFile, prjexFile, traditionOut, proposeOut;
	inputImage =
		".\\data\\original\\GF1_WFV1_E117.5_N24.6_20241129_L1A13664362001.tiff";
	rpbFile =
		".\\data\\RPC\\GF1_WFV1_E117.5_N24.6_20241129_L1A13664362001.rpb";
	rpbexFile =
		".\\data\\RPC\\GF1_WFV1_E117.5_N24.6_20241129_L1A13664362001.rpbex";
	demFile =
		".\\data\\fujian_DEM_30m.tif";
	prjexFile =
		".\\data\\RPC\\tsea.rpbex.prj";
	proposeOut =
		".\\data\\output\\GF1_WFV1_E117.5_N24.6_20241129_L1A13664362001_Hex_direct.tif";
	traditionOut =
		".\\data\\output\\GF1_WFV1_E117.5_N24.6_20241129_L1A13664362001_Hex_indire.tif";
	// Extended RPC generation and image orthorectification.
	//bool ifmake = rpcmake(inputImage, rpbFile, rpbexFile);
	double proposedTime = proposedMethod(inputImage, rpbexFile, demFile, prjexFile, proposeOut, BILINEAR);
	double traditionalTime = traditionalMethod(inputImage, rpbFile, demFile, prjexFile, traditionOut, BILINEAR);

	// Quality assessment.
	//double* globalResult = globalQA(inputImage, proposeOut, traditionOut);
	//double* localResult = localQA(inputImage, proposeOut, traditionOut, rpbexFile, demFile, prjexFile);
	//int* imageLocat = GeoCalculate_Direct(14213, 3226, proposeOut, rpbexFile, prjexFile, demFile);
	//int* imageLocat2 = GeoCalculate_viaLL(14211, 3214, traditionOut, rpbFile, prjexFile, demFile);

	return EXIT_SUCCESS;
}