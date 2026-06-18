#include "imaBack.h"

int* GeoCalculate_Direct(
	int samp, int line,
	std::string inputImage,
	std::string rpbexFile,
	std::string prjexFile,
	std::string demFile)
{
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");

	std::shared_ptr<TifRPCDataSet> tifRpcDS = std::make_shared<TifRPCDataSet>();
	std::shared_ptr<CGeoRaster> demDS = std::make_shared<CGeoRaster>();
	if (!demDS->readMetaData(demFile)) {
		std::cerr << "Read dem metadata failed: " << demFile << std::endl;
		exit(0);
	}

	if (!tifRpcDS->loadImage(inputImage) && !tifRpcDS->readRPB(rpbexFile)) {
		std::cerr << "Load source image failed: " << inputImage << std::endl;
		exit(0);
	}
	else {
		tifRpcDS->readRPB(rpbexFile);
		std::cout << "Source image loaded: " << inputImage << std::endl;
	}

	OGRSpatialReference eoSRS;
	if (GDAL_VERSION_MAJOR >= 3)
		eoSRS.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);

	std::ifstream is;
	is.open(prjexFile);
	if (is.fail()) {
		std::cout << "prj failed" << endl;
		exit(0);
	}

	std::stringstream buffer;
	buffer << is.rdbuf();
	std::string srsWKT(buffer.str());
	is.close();

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

	double srcImageTrans[6];
	poSrcDS->GetGeoTransform(srcImageTrans);
	double xMin = srcImageTrans[0];
	double yMax = srcImageTrans[3];
	double dx = srcImageTrans[1];
	double dy = srcImageTrans[5];

	double xMax = xMin + dx * nWidth;
	double yMin = yMax + dy * nHeight;

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

	GDALDriver* poDriver;
	poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
	if (poDriver == NULL) {
		CPLFree(pSrcData);
		exit(0);
	}

	const int level_ = 18;
	RtSphHexDGGS imaDGGS = RtSphHexDGGS(6371, level_);

	unsigned short  rtFaceId = -1;
	std::string tesa;
	std::stringstream ss(srsWKT);
	ss >> tesa >> rtFaceId;
	assert(rtFaceId >= 0 && rtFaceId < 60);

	double orthLogicHexX = imaDGGS.lengthK_ / 2 * 1000;
	double orthLogicHexY = imaDGGS.lengthH_ / 2 * 1000;

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

	std::cout << "input pixel: samp=" << samp << ",line=" << line << std::endl;

	double east = 0;
	if (line % 2)
	{
		if (samp % 2) { east = xMin + dx * samp; }
		else { east = xMin + dx * (samp - 1); }
	}
	else
	{
		if (samp % 2) { east = xMin + dx * (samp - 1); }
		else { east = xMin + dx * samp; }
	}
	double south = yMax + dy * (line - 1);

	DGridCoord2D ccxy(rtFaceId, east * 0.001, south * 0.001);
	GeoCoord gpt = imaDGGS.RTInvSlice(ccxy);
	double demx = gpt.lon * RAD_TO_DEG;
	double demy = gpt.lat * RAD_TO_DEG;
	poTransLatLon2Dem->Transform(1, &demx, &demy);

	double avgh = tifRpcDS->meanTerrainHeight();
	double hz = avgh;
	if (!InterpolateZ(demx, demy, &hz, demAOI, demDataAOI))
		hz = avgh;

	RPCSensorModel* sm = dynamic_cast<RPCSensorModel*>(tifRpcDS->sensorModel());

	double outSamp = 0, outLine = 0;
	sm->Conv_Map_To_SL(east, south, hz, &outSamp, &outLine);

	int oSamp = (int)outSamp + 1;
	int oLine = (int)outLine + 1;

	std::cout << "source pixel: samp=" << oSamp << ",line=" << oLine << std::endl;

	int* location;
	location = (int*)VSIMalloc2(1, 2 * sizeof(int));
	location[0] = oSamp;
	location[1] = oLine;

	return location;
}

int* GeoCalculate_viaLL(int samp, int line, std::string inputImage, std::string rpbFile, std::string prjexFile, std::string demFile)
{
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");

	std::shared_ptr<TifRPCDataSet> tifRpcDS = std::make_shared<TifRPCDataSet>();
	std::shared_ptr<CGeoRaster> demDS = std::make_shared<CGeoRaster>();
	if (!demDS->readMetaData(demFile)) {
		std::cerr << "Read dem metadata failed: " << demFile << std::endl;
		exit(0);
	}

	if (!tifRpcDS->loadImage(inputImage) && !tifRpcDS->readRPB(rpbFile)) {
		std::cerr << "Load source image failed: " << inputImage << std::endl;
		exit(0);
	}
	else {
		tifRpcDS->readRPB(rpbFile);
		std::cout << "Source image loaded: " << inputImage << std::endl;
	}

	OGRSpatialReference eoSRS;
	if (GDAL_VERSION_MAJOR >= 3)
		eoSRS.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);

	std::ifstream is;
	is.open(prjexFile);
	if (is.fail()) {
		std::cout << "prj failed" << endl;
		exit(0);
	}

	std::stringstream buffer;
	buffer << is.rdbuf();
	std::string srsWKT(buffer.str());
	is.close();

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

	double srcImageTrans[6];
	poSrcDS->GetGeoTransform(srcImageTrans);
	double xMin = srcImageTrans[0];
	double yMax = srcImageTrans[3];
	double dx = srcImageTrans[1];
	double dy = srcImageTrans[5];

	double xMax = xMin + dx * nWidth;
	double yMin = yMax + dy * nHeight;

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

	GDALDriver* poDriver;
	poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
	if (poDriver == NULL) {
		CPLFree(pSrcData);
		exit(0);
	}

	const int level_ = 18;
	RtSphHexDGGS imaDGGS = RtSphHexDGGS(6371, level_);

	unsigned short  rtFaceId = -1;
	std::string tesa;
	std::stringstream ss(srsWKT);
	ss >> tesa >> rtFaceId;
	assert(rtFaceId >= 0 && rtFaceId < 60);

	double orthLogicHexX = imaDGGS.lengthK_ / 2 * 1000;
	double orthLogicHexY = imaDGGS.lengthH_ / 2 * 1000;

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

	std::cout << "input pixel: samp=" << samp << ",line=" << line << std::endl;

	double east = 0;
	if (line % 2)
	{
		if (samp % 2) { east = xMin + dx * samp; }
		else { east = xMin + dx * (samp - 1); }
	}
	else
	{
		if (samp % 2) { east = xMin + dx * (samp - 1); }
		else { east = xMin + dx * samp; }
	}
	double south = yMax + dy * (line - 1);

	DGridCoord2D ccxy(rtFaceId, east * 0.001, south * 0.001);
	GeoCoord gpt = imaDGGS.RTInvSlice(ccxy);
	double demx = gpt.lon * RAD_TO_DEG;
	double demy = gpt.lat * RAD_TO_DEG;
	poTransLatLon2Dem->Transform(1, &demx, &demy);

	double avgh = tifRpcDS->meanTerrainHeight();
	double hz = avgh;
	if (!InterpolateZ(demx, demy, &hz, demAOI, demDataAOI))
		hz = avgh;

	RPCSensorModel* sm = dynamic_cast<RPCSensorModel*>(tifRpcDS->sensorModel());

	double outSamp = 0, outLine = 0;
	sm->Conv_Map_To_SL(gpt.lon * RAD_TO_DEG, gpt.lat * RAD_TO_DEG, hz, &outSamp, &outLine);

	int oSamp = (int)outSamp + 1;
	int oLine = (int)outLine + 1;

	std::cout << "source pixel: samp=" << oSamp << ",line=" << oLine << std::endl;

	int* location;
	location = (int*)VSIMalloc2(1, 2 * sizeof(int));
	location[0] = oSamp;
	location[1] = oLine;

	return location;
}