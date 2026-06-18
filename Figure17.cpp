#include "Figure17.h"

void Figure17(
	const std::string originImage,
	const std::string HexDirec,
	const std::string HexIndir,
	const std::string rpbFile,
	const std::string demFile,
	const std::string prjFile)
{
	double* out = localQA(originImage, HexDirec, HexIndir, rpbFile, demFile, prjFile);
	ofstream fout;
	fout.open("Figure17.txt", ios::out);
	fout << "directMSE: " << out[2] << ", indireMSE: " << out[3] << std::endl;
	fout.close();
}