#include "experiment.h"
#include <sstream>
#include<assert.h>

double* localQA(
	const std::string originImage,
	const std::string HexDirec,
	const std::string HexIndir,
	const std::string rpbFile,
	const std::string demFile,
	const std::string prjFile)
{
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");

	std::shared_ptr<TifRPCDataSet> tifRpcDS = std::make_shared<TifRPCDataSet>();
	if (!tifRpcDS->loadImage(originImage) && !tifRpcDS->readRPB(rpbFile)) {
		std::cerr << "Load source image failed: " << originImage << std::endl;
		exit(0);
	}
	else {
		tifRpcDS->readRPB(rpbFile);
		std::cout << "Source image loaded: " << originImage << std::endl;
	}

	std::shared_ptr<CGeoRaster> demDS = std::make_shared<CGeoRaster>();
	if (!demDS->readMetaData(demFile)) {
		std::cerr << "Read dem metadata failed: " << demFile << std::endl;
		exit(0);
	}

	OGRSpatialReference eoSRS;
	if (GDAL_VERSION_MAJOR >= 3)
		eoSRS.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);

	std::ifstream is;
	is.open(prjFile);
	if (is.fail()) {
		std::cout << "prj failed" << endl;
		exit(0);
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

	std::string srcImageFile = tifRpcDS->imageFileName();
	GDALDataset* poSrcDS = (GDALDataset*)GDALOpen(srcImageFile.c_str(), GA_ReadOnly);
	if (poSrcDS == NULL) {
		std::cout << "Source image can not be opened:" << srcImageFile << std::endl;
		exit(0);
	}

	int nWidth = poSrcDS->GetRasterXSize();
	int nHeight = poSrcDS->GetRasterYSize();
	int nBands = poSrcDS->GetRasterCount();
	GDALDataType nDataType = poSrcDS->GetRasterBand(1)->GetRasterDataType();
	int nDepth = GDALGetDataTypeSize((GDALDataType)nDataType) / 8;

	if (nDataType != GDT_Byte && nDataType != GDT_UInt16 && nDataType != GDT_Int16 && nDataType != GDT_Float32) {
		std::cout << "Unsupported data bits for source image." << std::endl;
		GDALClose(poSrcDS);
		exit(0);
	}

	float* pSrcData;
	pSrcData = (float*)VSIMalloc2(nWidth, nHeight * nBands * sizeof(float));
	if (pSrcData == NULL) {
		GDALClose(poSrcDS);
		exit(0);
	}

	poSrcDS->RasterIO(GF_Read,
		0, 0,
		nWidth, nHeight,
		(void*)pSrcData,
		nWidth, nHeight,
		GDT_Float32, nBands, NULL, 0, 0, 0);

	GDALClose(poSrcDS);

	GDALDataset* direcHexDS = (GDALDataset*)GDALOpen(HexDirec.c_str(), GA_ReadOnly);
	if (direcHexDS == NULL) {
		std::cout << "HexDirect image can not be opened:" << HexDirec << std::endl;
		exit(0);
	}

	int dWidth = direcHexDS->GetRasterXSize();
	int dHeight = direcHexDS->GetRasterYSize();
	int dBands = direcHexDS->GetRasterCount();

	double hexImageTrans[6];
	direcHexDS->GetGeoTransform(hexImageTrans);
	double xMin = hexImageTrans[0];
	double yMax = hexImageTrans[3];
	double dx = hexImageTrans[1];
	double dy = hexImageTrans[5];
	double xMax = xMin + dx * dWidth;
	double yMin = yMax + dy * dHeight;

	float* hDirecData;
	hDirecData = (float*)VSIMalloc2(dWidth, dHeight * dBands * sizeof(float));
	if (hDirecData == NULL) {
		GDALClose(direcHexDS);
		exit(0);
	}

	direcHexDS->RasterIO(GF_Read,
		0, 0,
		dWidth, dHeight,
		(void*)hDirecData,
		dWidth, dHeight,
		GDT_Float32, dBands, NULL, 0, 0, 0);

	std::cout << "Hex Direct image readed as float: " << HexDirec << std::endl;

	GDALClose(direcHexDS);

	GDALDataset* indirHexDS = (GDALDataset*)GDALOpen(HexIndir.c_str(), GA_ReadOnly);
	if (indirHexDS == NULL) {
		std::cout << "HexIndirect image can not be opened:" << HexIndir << std::endl;
		exit(0);
	}

	int iWidth = indirHexDS->GetRasterXSize();
	int iHeight = indirHexDS->GetRasterYSize();
	int iBands = indirHexDS->GetRasterCount();

	double hexIndirecImageTrans[6];
	indirHexDS->GetGeoTransform(hexIndirecImageTrans);
	double xMinIn = hexIndirecImageTrans[0];
	double yMaxIn = hexIndirecImageTrans[3];

	float* hIndirData;
	hIndirData = (float*)VSIMalloc2(iWidth, iHeight * iBands * sizeof(float));
	if (hIndirData == NULL) {
		GDALClose(indirHexDS);
		exit(0);
	}

	indirHexDS->RasterIO(GF_Read,
		0, 0,
		iWidth, iHeight,
		(void*)hIndirData,
		iWidth, iHeight,
		GDT_Float32, iBands, NULL, 0, 0, 0);

	std::cout << "Hex Indirect image readed as float: " << HexIndir << std::endl;

	GDALClose(indirHexDS);

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
		exit(0);
	}

	const int level_ = 18;
	RtSphHexDGGS imaDGGS = RtSphHexDGGS(6371, level_);

	DGridCoord2D ixy(rtFaceId, xMin * 0.001, yMax * 0.001);
	double uldemx = imaDGGS.RTInvSlice(ixy).lon * RAD_TO_DEG;
	double uldemy = imaDGGS.RTInvSlice(ixy).lat * RAD_TO_DEG;
	poTransLatLon2Dem->Transform(1, &uldemx, &uldemy);

	DGridCoord2D lxy(rtFaceId, xMax * 0.001, yMin * 0.001);
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
		exit(0);
	}

	if (!demDS->readBlock(ulDemAOIRow, ulDemAOICol, demAOI.col, demAOI.row, (void*)demDataAOI, GDT_Float32)) {
		std::cerr << "Failed to load the elevation data in dem aoi: " << demDS->imageFileName() << std::endl;
		GDALClose(poSrcDS);
		CPLFree(pSrcData);
		CPLFree(demDataAOI);
		exit(0);
	}
	else
		std::cout << "DEM aoi loaded. " << std::endl;

	double east = 0, south = 0;

	RPCSensorModel* sm = dynamic_cast<RPCSensorModel*>(tifRpcDS->sensorModel());

	double sumD = 0, ssumD = 0;
	double sum = std::numeric_limits<double>::max();
	double ssum = std::numeric_limits<double>::max();
	int tempS = 0, tempL = 0;

	double originSamp = 0, originLine = 0;
	int oSamp = 0, oLine = 0;

	int originI = 0, originJ = 0;

	int inHexJ = 0, inHexI = 0;

	int Sr[11];
	int Di[11];
	int In[11];

	double Miu0 = 0, Miu1 = 0, Miu2 = 0;
	double theta0 = 0, theta1 = 0, theta2 = 0, conv1 = 0, conv2 = 0;
	double C1 = 1;
	double C2 = 1;
	double SSim1 = 0, SSim2 = 0;

	double tempDiff = 0, Differ1 = 0, Differ2 = 0;
	double tempSqua = 0, SquarE1 = 0, SquarE2 = 0;

	int onePerc = static_cast<int>((13800 - 2900) * 0.01), nPrec = 0;

	std::cout << "Begin the experiment. " << std::endl;
	for (int i = 2900; i < 13800; i += 4)
	{
		south = yMax + dy * i;

		for (int j = 4400; j < 22800; j += 8)
		{
			if (i % 2)
			{
				if (j % 2) { east = xMin + dx * (j + 1); }
				else { east = xMin + dx * j; }
			}
			else
			{
				if (j % 2) { east = xMin + dx * j; }
				else { east = xMin + dx * (j + 1); }
			}

			DGridCoord2D ccxy(rtFaceId, east * 0.001, south * 0.001);
			GeoCoord gpt = imaDGGS.RTInvSlice(ccxy);
			double demx = gpt.lon * RAD_TO_DEG;
			double demy = gpt.lat * RAD_TO_DEG;
			poTransLatLon2Dem->Transform(1, &demx, &demy);

			double avgh = tifRpcDS->meanTerrainHeight();
			double hz = avgh;
			if (!InterpolateZ(demx, demy, &hz, demAOI, demDataAOI))
				hz = avgh;

			sm->Conv_Map_To_SL(east, south, hz, &originSamp, &originLine);
			oSamp = (int)originSamp + 1;
			oLine = (int)originLine + 1;

			for (int teI = -1; teI < 2; teI++)
				for (int teJ = -1; teJ < 2; teJ++)
				{
					tempS = oSamp + teJ;
					tempL = oLine + teI;

					sumD += fabs(hDirecData[j + i * dWidth] - pSrcData[tempS + tempL * nWidth]);
					sumD += fabs(hDirecData[j + (i - 2) * dWidth] - pSrcData[tempS + 2 + tempL * nWidth]);
					sumD += fabs(hDirecData[j + (i + 2) * dWidth] - pSrcData[tempS - 2 + tempL * nWidth]);

					sumD += fabs(hDirecData[j - 2 + i * dWidth] - pSrcData[tempS + (tempL - 1) * nWidth]);
					sumD += fabs(hDirecData[j - 2 + (i - 2) * dWidth] - pSrcData[tempS + 2 + (tempL - 1) * nWidth]);
					sumD += fabs(hDirecData[j - 2 + (i + 2) * dWidth] - pSrcData[tempS - 2 + (tempL - 1) * nWidth]);

					sumD += fabs(hDirecData[j + 2 + i * dWidth] - pSrcData[tempS + (tempL + 1) * nWidth]);
					sumD += fabs(hDirecData[j + 2 + (i - 2) * dWidth] - pSrcData[tempS + 2 + (tempL + 1) * nWidth]);
					sumD += fabs(hDirecData[j + 2 + (i + 2) * dWidth] - pSrcData[tempS - 2 + (tempL + 1) * nWidth]);
					if (sumD < sum)
					{
						sum = sumD;
						originI = tempL;
						originJ = tempS;
					}
					sumD = 0;
				}
			sum = std::numeric_limits<double>::max();

			for (int teI = -20; teI < 21; teI++)
				for (int teJ = -10; teJ < 11; teJ += 2)
				{
					tempS = int((east - xMinIn) / dx) + teJ;
					tempL = int((south - yMaxIn) / dy) + teI;
					ssumD += fabs(hDirecData[j + i * dWidth] - hIndirData[tempS + tempL * iWidth]);
					ssumD += fabs(hDirecData[j - 2 + i * dWidth] - hIndirData[tempS - 2 + tempL * iWidth]);
					ssumD += fabs(hDirecData[j + 2 + i * dWidth] - hIndirData[tempS + 2 + tempL * iWidth]);

					ssumD += fabs(hDirecData[j + (i - 1) * dWidth] - hIndirData[tempS + (tempL - 1) * iWidth]);
					ssumD += fabs(hDirecData[j - 2 + (i - 1) * dWidth] - hIndirData[tempS - 2 + (tempL - 1) * iWidth]);
					ssumD += fabs(hDirecData[j + 2 + (i - 1) * dWidth] - hIndirData[tempS + 2 + (tempL - 1) * iWidth]);

					ssumD += fabs(hDirecData[j + (i + 1) * dWidth] - hIndirData[tempS + (tempL + 1) * iWidth]);
					ssumD += fabs(hDirecData[j - 2 + (i + 1) * dWidth] - hIndirData[tempS - 2 + (tempL + 1) * iWidth]);
					ssumD += fabs(hDirecData[j + 2 + (i + 1) * dWidth] - hIndirData[tempS + 2 + (tempL + 1) * iWidth]);
					if (ssumD < ssum)
					{
						ssum = ssumD;
						inHexI = tempL;
						inHexJ = tempS;
					}
					ssumD = 0;
				}
			ssum = std::numeric_limits<double>::max();

			Sr[0] = pSrcData[originJ + (originI)*nWidth];
			Sr[1] = pSrcData[originJ + (originI - 1) * nWidth];
			Sr[2] = pSrcData[originJ + (originI - 2) * nWidth];
			Sr[3] = pSrcData[originJ + (originI + 1) * nWidth];
			Sr[4] = pSrcData[originJ + (originI + 2) * nWidth];

			Sr[5] = pSrcData[originJ + 2 + (originI)*nWidth];
			Sr[6] = pSrcData[originJ + 2 + (originI + 1) * nWidth];
			Sr[7] = pSrcData[originJ + 2 + (originI - 1) * nWidth];

			Sr[8] = pSrcData[originJ - 2 + (originI)*nWidth];
			Sr[9] = pSrcData[originJ - 2 + (originI + 1) * nWidth];
			Sr[10] = pSrcData[originJ - 2 + (originI - 1) * nWidth];

			Di[0] = hDirecData[j + (i)*dWidth];
			Di[1] = hDirecData[j - 2 + i * dWidth];
			Di[2] = hDirecData[j - 4 + i * dWidth];
			Di[3] = hDirecData[j + 2 + i * dWidth];
			Di[4] = hDirecData[j + 4 + i * dWidth];

			Di[5] = hDirecData[j + (i - 2) * dWidth];
			Di[6] = hDirecData[j + 2 + (i - 2) * dWidth];
			Di[7] = hDirecData[j - 2 + (i - 2) * dWidth];

			Di[8] = hDirecData[j + (i + 2) * dWidth];
			Di[9] = hDirecData[j + 2 + (i + 2) * dWidth];
			Di[10] = hDirecData[j - 2 + (i + 2) * dWidth];

			In[0] = hIndirData[inHexJ + (inHexI)*iWidth];
			In[1] = hIndirData[inHexJ - 2 + inHexI * iWidth];
			In[2] = hIndirData[inHexJ - 4 + inHexI * iWidth];
			In[3] = hIndirData[inHexJ + 2 + inHexI * iWidth];
			In[4] = hIndirData[inHexJ + 4 + inHexI * iWidth];

			In[5] = hIndirData[inHexJ + (inHexI - 2) * iWidth];
			In[6] = hIndirData[inHexJ + 2 + (inHexI - 2) * iWidth];
			In[7] = hIndirData[inHexJ - 2 + (inHexI - 2) * iWidth];

			In[8] = hIndirData[inHexJ + (inHexI + 2) * iWidth];
			In[9] = hIndirData[inHexJ + 2 + (inHexI + 2) * iWidth];
			In[10] = hIndirData[inHexJ - 2 + (inHexI + 2) * iWidth];

			for (int k = 0; k < 11; k++)
			{
				tempDiff += fabs(Di[k] - Sr[k]);
			}
			Differ1 += tempDiff / 11;
			tempDiff = 0;

			for (int k = 0; k < 11; k++)
			{
				tempDiff += fabs(In[k] - Sr[k]);
			}
			Differ2 += tempDiff / 11;
			tempDiff = 0;

			for (int k = 0; k < 11; k++)
			{
				tempSqua += pow(Di[k] - Sr[k], 2);
			}
			SquarE1 += tempSqua / 11;
			tempSqua = 0;

			for (int k = 0; k < 11; k++)
			{
				tempSqua += pow(In[k] - Sr[k], 2);
			}
			SquarE2 += tempSqua / 11;
			tempSqua = 0;

			for (int k = 0; k < 11; k++)
			{
				Miu0 += Sr[k];
				Miu1 += Di[k];
				Miu2 += In[k];
			}
			Miu0 /= 11; Miu1 /= 11; Miu2 /= 11;
			for (int k = 0; k < 11; k++)
			{
				theta0 += pow(Sr[k] - Miu0, 2);
				theta1 += pow(Di[k] - Miu1, 2);
				theta2 += pow(In[k] - Miu2, 2);
				conv1 += (Di[k] - Miu1) * (Sr[k] - Miu0);
				conv2 += (In[k] - Miu2) * (Sr[k] - Miu0);
			}
			theta0 /= 11; theta1 /= 11; theta2 /= 11; conv1 /= 11; conv2 /= 11;

			SSim1 += (2 * Miu0 * Miu1 + C1) * (2 * conv1 + C2) / ((pow(Miu0, 2) + pow(Miu1, 2) + C1) * (theta0 + theta1 + C2));
			SSim2 += (2 * Miu0 * Miu2 + C1) * (2 * conv2 + C2) / ((pow(Miu0, 2) + pow(Miu2, 2) + C1) * (theta0 + theta2 + C2));

			Miu0 = 0; Miu1 = 0; Miu2 = 0;
			theta0 = 0; theta1 = 0; theta2 = 0;
			conv1 = 0; conv2 = 0;
		}

		if (i >= (2900 + nPrec * onePerc)) {
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

	Differ1 /= (2725 * 2300);
	Differ2 /= (2725 * 2300);
	std::cout << endl;
	std::cout << "directMAE:" << Differ1 << endl;
	std::cout << "indireMAE:" << Differ2 << endl;

	SquarE1 /= (2725 * 2300);
	SquarE2 /= (2725 * 2300);
	std::cout << "directMSE:" << SquarE1 << endl;
	std::cout << "indireMSE:" << SquarE2 << endl;

	SSim1 /= (2725 * 2300);
	SSim2 /= (2725 * 2300);
	std::cout << "directSSIM:" << SSim1 << endl;
	std::cout << "indireSSIM:" << SSim2 << endl;

	double* local;
	local = (double*)VSIMalloc2(1, 6 * sizeof(double));
	local[0] = Differ1;
	local[1] = Differ2;
	local[2] = SquarE1;
	local[3] = SquarE2;
	local[4] = SSim1;
	local[5] = SSim2;

	CPLFree(pSrcData);
	CPLFree(hDirecData);
	CPLFree(hIndirData);
	CPLFree(demDataAOI);

	return local;
}

bool meanAndStan(const double* src, const int& width, const int& height, const int& bands, double* mean, double* stand)
{
	std::cout << "begin to compute mean: \n";

	double* sumBand;
	sumBand = (double*)VSIMalloc2(1, (bands + 1) * sizeof(double));

	memset(sumBand, 0, sizeof(double) * (bands + 1));

	int onePerc = static_cast<int>(height * 0.01), nPrec = 0;

	bool value = 0;
	int nonZeroNum = 0;

	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			for (int k = 0; k < bands; ++k) {
				size_t pos = k * height * width;
				if (src[pos + i * width + j]) { value = 1; break; }
			}
			if (value == 0) { continue; }
			nonZeroNum += 1;
			value = 0;

			for (int k = 0; k < bands; ++k) {
				size_t pos = k * height * width;
				sumBand[k] += src[pos + i * width + j];
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

	for (int i = 0; i < bands; i++) {
		sumBand[i] /= (nonZeroNum);
	}

	for (int i = 0; i < bands; i++) {
		sumBand[bands] += sumBand[i];
	}
	sumBand[bands] /= bands;

	for (int i = 0; i < bands + 1; i++) {
		mean[i] = sumBand[i];
	}
	VSIFree(sumBand);

	std::cout << "\ncompute mean successfully !\n";

	onePerc = static_cast<int>(height * 0.01);
	nPrec = 0;

	std::cout << "begin to compute standard: \n";

	double* standDev;
	standDev = (double*)VSIMalloc2(1, (bands + 1) * sizeof(double));
	memset(standDev, 0, sizeof(double) * (bands + 1));

	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			for (int k = 0; k < bands; ++k) {
				size_t pos = k * height * width;
				if (src[pos + i * width + j]) { value = 1; break; }
			}
			if (value == 0) { continue; }
			value = 0;

			for (int k = 0; k < bands; ++k) {
				size_t pos = k * height * width;
				standDev[k] += pow(src[pos + i * width + j] - mean[k], 2);
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

	for (int i = 0; i < bands; i++) {
		standDev[i] = sqrt(standDev[i] / (nonZeroNum));
	}

	for (int i = 0; i < bands; i++) {
		standDev[bands] += standDev[i];
	}
	standDev[bands] /= bands;

	for (int i = 0; i < bands + 1; i++) {
		stand[i] = standDev[i];
	}
	VSIFree(standDev);
	std::cout << "\ncompute standard successfully !\n";

	return true;
}

bool entropyMeanGrad(const bool ifhex, const double* src, const int& width, const int& height, const int& bands, double* mg)
{
	std::cout << "begin to compute meanGrad: \n";

	double* sumBand;
	sumBand = (double*)VSIMalloc2(1, (bands + 1) * sizeof(double));

	memset(sumBand, 0, sizeof(double) * (bands + 1));

	int onePerc = static_cast<int>(height * 0.01), nPrec = 0;

	bool value = 0;
	int nonZeroNum = 0;

	if (ifhex) {
		for (int i = 0; i < height - 1; i += 2) {
			for (int j = 0; j < width - 1; j += 2) {
				for (int k = 0; k < bands; ++k) {
					size_t pos = k * height * width;
					if (src[pos + i * width + j] && src[pos + i * width + j + 1] && src[pos + (i + 1) * width + j])
					{
						value = 1; break;
					}
				}
				if (value == 0) { continue; }
				nonZeroNum += 1;
				value = 0;

				bool jiou = 0;
				if (src[i * width + j] == src[i * width + j + 1]
					&& src[height * width + i * width + j] == src[height * width + i * width + j + 1]
					&& src[2 * height * width + i * width + j] == src[2 * height * width + i * width + j + 1])
				{
					jiou = 1;
				}
				for (int k = 0; k < bands; ++k) {
					size_t pos = k * height * width;
					double aver = 0;
					if (jiou) {
						aver = (src[pos + (i + 1) * width + j] + src[pos + (i + 1) * width + j + 1]) / 2;
					}
					else {
						aver = (src[pos + (i + 1) * width + j] + src[pos + (i + 1) * width + j - 1]) / 2;
					}
					sumBand[k] +=
						sqrt(
							(
								pow(src[pos + i * width + j + 2] - src[pos + i * width + j], 2)
								+ pow(aver - src[pos + i * width + j], 2)
								) / 2
						);
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
	}

	else {
		for (int i = 0; i < height - 1; i++) {
			for (int j = 0; j < width - 1; j++) {
				for (int k = 0; k < bands; ++k) {
					size_t pos = k * height * width;
					if (src[pos + i * width + j] && src[pos + i * width + j + 1] && src[pos + (i + 1) * width + j])
					{
						value = 1; break;
					}
				}
				if (value == 0) { continue; }
				nonZeroNum += 1;
				value = 0;

				for (int k = 0; k < bands; ++k) {
					size_t pos = k * height * width;
					sumBand[k] +=
						sqrt(
							(
								pow(src[pos + i * width + j + 1] - src[pos + i * width + j], 2)
								+ pow(src[pos + (i + 1) * width + j] - src[pos + i * width + j], 2)
								) / 2
						);
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
	}

	for (int i = 0; i < bands; i++) {
		sumBand[i] /= (nonZeroNum);
	}

	for (int i = 0; i < bands; i++) {
		sumBand[bands] += sumBand[i];
	}
	sumBand[bands] /= bands;

	for (int i = 0; i < bands + 1; i++) {
		mg[i] = sumBand[i];
	}
	VSIFree(sumBand);

	std::cout << "\ncompute mg successfully !\n";

	return true;
}

double* globalQA(
	const std::string srcImagePath,
	const std::string HexDirec,
	const std::string HexIndir)
{
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");

	const int level_ = 18;
	RtSphHexDGGS imaDGGS = RtSphHexDGGS(6371, level_);
	GDALDataset* poSrcDS = (GDALDataset*)GDALOpen(srcImagePath.c_str(), GA_ReadOnly);
	if (poSrcDS == NULL) {
		std::cout << "Source image can not be opened:" << srcImagePath << std::endl;
		exit(0);
	}
	else {
		std::cout << "Source image loaded: " << srcImagePath << std::endl;
	}

	int sWidth = poSrcDS->GetRasterXSize();
	int sHeight = poSrcDS->GetRasterYSize();
	int sBands = poSrcDS->GetRasterCount();

	double* pSrcData;
	pSrcData = (double*)VSIMalloc2(sWidth, sHeight * sBands * sizeof(double));
	if (pSrcData == NULL) {
		GDALClose(poSrcDS);
		exit(0);
	}

	poSrcDS->RasterIO(GF_Read,
		0, 0,
		sWidth, sHeight,
		(void*)pSrcData,
		sWidth, sHeight,
		GDT_Float64, sBands, NULL, 0, 0, 0);

	GDALDataset* poDirectDS = (GDALDataset*)GDALOpen(HexDirec.c_str(), GA_ReadOnly);
	if (poDirectDS == NULL) {
		std::cout << "Direct image can not be opened:" << HexDirec << std::endl;
		exit(0);
	}
	else {
		std::cout << "Direct image loaded: " << HexDirec << std::endl;
	}

	int fWidth = poDirectDS->GetRasterXSize();
	int fHeight = poDirectDS->GetRasterYSize();
	int fBands = poDirectDS->GetRasterCount();

	double* pDirectData;
	pDirectData = (double*)VSIMalloc2(fWidth, fHeight * fBands * sizeof(double));
	if (pDirectData == NULL) {
		GDALClose(poDirectDS);
		exit(0);
	}

	poDirectDS->RasterIO(GF_Read,
		0, 0,
		fWidth, fHeight,
		(void*)pDirectData,
		fWidth, fHeight,
		GDT_Float64, fBands, NULL, 0, 0, 0);

	GDALDataset* poIndireDS = (GDALDataset*)GDALOpen(HexIndir.c_str(), GA_ReadOnly);
	if (poIndireDS == NULL) {
		std::cout << "Indirect image can not be opened:" << HexIndir << std::endl;
		exit(0);
	}
	else {
		std::cout << "Indirect image loaded: " << HexIndir << std::endl;
	}

	int tWidth = poIndireDS->GetRasterXSize();
	int tHeight = poIndireDS->GetRasterYSize();
	int tBands = poIndireDS->GetRasterCount();

	double* pIndireData;
	pIndireData = (double*)VSIMalloc2(tWidth, tHeight * tBands * sizeof(double));
	if (pIndireData == NULL) {
		GDALClose(poIndireDS);
		exit(0);
	}

	poIndireDS->RasterIO(GF_Read,
		0, 0,
		tWidth, tHeight,
		(void*)pIndireData,
		tWidth, tHeight,
		GDT_Float64, tBands, NULL, 0, 0, 0);

	double* meanSrc;
	meanSrc = (double*)VSIMalloc2(1, (sBands + 1) * sizeof(double));
	double* standSrc;
	standSrc = (double*)VSIMalloc2(1, (sBands + 1) * sizeof(double));

	double* gradSrc;
	gradSrc = (double*)VSIMalloc2(1, (sBands + 1) * sizeof(double));

	if (meanAndStan(pSrcData, sWidth, sHeight, sBands, meanSrc, standSrc)) {
		std::cout << "Mean and standard deviation compute successfully!" << std::endl;

		std::cout << "Aver mean: " << meanSrc[sBands] << std::endl;

		std::cout << "Aver standard: " << standSrc[sBands] << std::endl;
	}

	if (entropyMeanGrad(0, pSrcData, sWidth, sHeight, sBands, gradSrc)) {
		std::cout << "Mean Grad compute successfully!" << std::endl;
		std::cout << "Mean Grad: " << gradSrc[sBands] << std::endl;
	}

	double* meanDirect;
	meanDirect = (double*)VSIMalloc2(1, (fBands + 1) * sizeof(double));
	double* standDirect;
	standDirect = (double*)VSIMalloc2(1, (fBands + 1) * sizeof(double));

	double* gradDirect;
	gradDirect = (double*)VSIMalloc2(1, (fBands + 1) * sizeof(double));

	if (meanAndStan(pDirectData, fWidth, fHeight, fBands, meanDirect, standDirect)) {
		std::cout << "Mean and standard deviation compute successfully!" << std::endl;

		std::cout << "Aver mean: " << meanDirect[fBands] << std::endl;

		std::cout << "Aver standard: " << standDirect[fBands] << std::endl;
	}

	if (entropyMeanGrad(1, pDirectData, fWidth, fHeight, fBands, gradDirect)) {
		std::cout << "Mean Grad compute successfully!" << std::endl;
		std::cout << "Mean Grad: " << gradDirect[fBands] << std::endl;
	}

	double direcMeanDif = abs(meanDirect[fBands] - meanSrc[sBands]);
	double direcStanDif = abs(standDirect[fBands] - standSrc[sBands]);
	double direcGradDif = abs(gradDirect[fBands] - gradSrc[sBands]);

	double* meanIndire;
	meanIndire = (double*)VSIMalloc2(1, (tBands + 1) * sizeof(double));
	double* standIndire;
	standIndire = (double*)VSIMalloc2(1, (tBands + 1) * sizeof(double));

	double* gradIndire;
	gradIndire = (double*)VSIMalloc2(1, (tBands + 1) * sizeof(double));

	if (meanAndStan(pIndireData, tWidth, tHeight, tBands, meanIndire, standIndire)) {
		std::cout << "Mean and standard deviation compute successfully!" << std::endl;

		std::cout << "Aver mean: " << meanIndire[tBands] << std::endl;

		std::cout << "Aver standard: " << standIndire[tBands] << std::endl;
	}

	if (entropyMeanGrad(1, pIndireData, tWidth, tHeight, tBands, gradIndire)) {
		std::cout << "Mean Grad compute successfully!" << std::endl;
		std::cout << "Mean Grad: " << gradIndire[tBands] << std::endl;
	}

	double indirMeanDif = abs(meanIndire[fBands] - meanSrc[sBands]);
	double indirStanDif = abs(standIndire[fBands] - standSrc[sBands]);
	double indirGradDif = abs(gradIndire[fBands] - gradSrc[sBands]);

	double* global;
	global = (double*)VSIMalloc2(1, 6 * sizeof(double));
	global[0] = direcMeanDif;
	global[1] = indirMeanDif;
	global[2] = direcStanDif;
	global[3] = indirStanDif;
	global[4] = direcGradDif;
	global[5] = indirGradDif;

	std::cout << "direcMeanDif: " << global[0] << std::endl;
	std::cout << "indirMeanDif: " << global[1] << std::endl;
	std::cout << "direcStanDif: " << global[2] << std::endl;
	std::cout << "indirStanDif: " << global[3] << std::endl;
	std::cout << "direcGradDif: " << global[4] << std::endl;
	std::cout << "indirGradDif: " << global[5] << std::endl;

	return global;
}