#include "Figure14andTable4.h"

void Figure14andTable4(
	const std::string srcImagePath,
	const std::string HexDirec,
	const std::string HexIndir)
{
	double* out = globalQA(srcImagePath, HexDirec, HexIndir);
	ofstream fout;
	fout.open("Figure14andTable4.txt", ios::out);
	fout << "direcVarDif: " << out[2] << ", indirVarDif: " << out[3] << std::endl;
	fout.close();
}