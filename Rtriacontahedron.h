#pragma once
#include "RtCommon.h"
#include "RtTriface.h"
#include "RtRhface.h"

class Rtriacontahedron
{
public:

	Rtriacontahedron(void);

	Rtriacontahedron(double l);

	void Setface();

protected:

	double L_;

	double a_;

	double b_;

	RtRhface rhomb_[30];

	RtTriface tri_[60];
};
