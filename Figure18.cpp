#include "Figure18.h"

void Figure18(
	const std::string originImage,
	const std::string HexDirec,
	const std::string HexIndir,
	const std::string rpbFile,
	const std::string demFile,
	const std::string prjFile)
{
	double* out = localQA(originImage, HexDirec, HexIndir, rpbFile, demFile, prjFile);
	ofstream fout;
	fout.open("Figure18.txt", ios::out);
	fout << "directSSIM: " << out[4] << ", indireSSIM: " << out[5] << std::endl;
	fout.close();
}