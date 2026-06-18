#include "Figure15andTable5.h"

void Figure15andTable5(
	const std::string srcImagePath,
	const std::string HexDirec,
	const std::string HexIndir)
{
	double* out = globalQA(srcImagePath, HexDirec, HexIndir);
	ofstream fout;
	fout.open("Figure15andTable5.txt", ios::out);
	fout << "direcGradDif: " << out[4] << ", indirGradDif: " << out[5] << std::endl;
	fout.close();
}