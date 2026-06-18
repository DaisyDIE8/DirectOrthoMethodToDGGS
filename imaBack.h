#pragma once
#include "georaster.h"
#include "tifrpcdataset.h"
#include "rpcsensormodel.h"
#include "RtA4HexDGGS.h"

#include "gdal.h"
#include "gdal_priv.h"
#include "cpl_conv.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <ctime>
#include<assert.h>

#include<omp.h>

/*
 * @brief	 Coordinate inverse calculation function for evaluating geometric errors.
 *
 * @param samp: The column number of a certain pixel on the output hexagonal image.
 * @param line: The row number of a certain pixel on the output hexagonal image.
 * @param inputImage: The path of the hexagonal orthophotograph file.
 * @param rpbFile: The path of the extended rpbex file.
 * @param prjexFile: The path of the extended projection file.
 * @param demFile: The path of the digital elevation model (DEM) file.
 *
 * @return	 An int-type array containing 2 values,
			 which describes the row and column numbers of the pixel on the original image
			 obtained through the inverse calculation of the input pixel.
 */
int* GeoCalculate_Direct(
	int samp, int line,
	std::string inputImage,
	std::string rpbexFile,
	std::string prjexFile,
	std::string demFile);

/*
 * @brief	 Coordinate inverse calculation function for evaluating geometric errors.
 *			 (conventional method via longitude & latitude)
 *
 * @param samp: The column number of a certain pixel on the output hexagonal image.
 * @param line: The row number of a certain pixel on the output hexagonal image.
 * @param inputImage: The path of the hexagonal orthophotograph file.
 * @param rpbFile: The path of the rpbex file.
 * @param prjexFile: The path of the projection file.
 * @param demFile: The path of the digital elevation model (DEM) file.
 *
 * @return	 An int-type array containing 2 values,
			 which describes the row and column numbers of the pixel on the original image
			 obtained through the inverse calculation of the input pixel.
 */
int* GeoCalculate_viaLL(
	int samp, int line,
	std::string inputImage,
	std::string rpbexFile,
	std::string prjexFile,
	std::string demFile);