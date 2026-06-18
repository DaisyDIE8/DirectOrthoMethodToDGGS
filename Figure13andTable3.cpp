#include "Figure13andTable3.h"

void Figure13andTable3(
	const std::string srcImagePath,
	const std::string HexDirec,
	const std::string HexIndir)
{
	double* out = globalQA(srcImagePath, HexDirec, HexIndir);
	ofstream fout;
	fout.open("Figure13andTable3.txt", ios::out);
	fout << "direcMeanDif: " << out[0] << ", indirMeanDif: " << out[1] << std::endl;
	fout.close();
}