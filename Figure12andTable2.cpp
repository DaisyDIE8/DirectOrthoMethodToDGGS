#include "Figure12andTable2.h"

void Figure12andTable2(
	int sampDirec, int lineDirec,
	int sampIndir, int lineIndir,
	std::string proposedImage,
	std::string traditionImage,
	std::string rpbexFile,
	std::string rpbFile,
	std::string prjexFile,
	std::string demFile)
{
	int* out1 = GeoCalculate_Direct(sampDirec, lineDirec, proposedImage, rpbexFile, prjexFile, demFile);
	int* out2 = GeoCalculate_viaLL(sampIndir, lineIndir, traditionImage, rpbFile, prjexFile, demFile);
	ofstream fout;
	fout.open("Figure12andTable2.txt", ios::out);
	fout << "proposed method  : source pixel: samp=" << out1[0] << ", line=" << out1[1] << std::endl;
	fout << "traditional method : source pixel: samp=" << out2[0] << ", line=" << out2[1] << std::endl;
	fout.close();
}