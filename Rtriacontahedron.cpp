#include "stdafx.h"
#include "Rtriacontahedron.h"

Rtriacontahedron::Rtriacontahedron(void)
{
	L_ = sqrt(M_PI / (3 * sqrt(5)));
	a_ = 2 * L_ / sqrt(1 + 2.0 / (3 + sqrt(5)));
	b_ = 2 * L_ / sqrt((5 + sqrt(5)) / 2.0);

	Setface();
}

Rtriacontahedron::Rtriacontahedron(double l)
{
	L_ = l;
	a_ = 2 * L_ / sqrt(1 + 2.0 / (3 + sqrt(5)));
	b_ = 2 * L_ / sqrt((5 + sqrt(5)) / 2.0);

	Setface();
}

void Rtriacontahedron::Setface()
{
	for (int i = 0; i < 30; i++)
		rhomb_[i].SetID(i);

	rhomb_[0].SetVertsID(0, 1, 2, 7);
	rhomb_[1].SetVertsID(0, 2, 3, 8);
	rhomb_[2].SetVertsID(0, 3, 4, 9);
	rhomb_[3].SetVertsID(0, 4, 5, 10);
	rhomb_[4].SetVertsID(0, 5, 1, 6);

	rhomb_[5].SetVertsID(6, 11, 1, 7);
	rhomb_[6].SetVertsID(7, 12, 2, 8);
	rhomb_[7].SetVertsID(8, 13, 3, 9);
	rhomb_[8].SetVertsID(9, 14, 4, 10);
	rhomb_[9].SetVertsID(10, 15, 5, 6);

	rhomb_[10].SetVertsID(7, 11, 16, 21);
	rhomb_[11].SetVertsID(7, 16, 12, 22);
	rhomb_[12].SetVertsID(8, 12, 17, 22);
	rhomb_[13].SetVertsID(8, 17, 13, 23);
	rhomb_[14].SetVertsID(9, 13, 18, 23);

	rhomb_[15].SetVertsID(9, 18, 14, 24);
	rhomb_[16].SetVertsID(10, 14, 19, 24);
	rhomb_[17].SetVertsID(10, 19, 15, 25);
	rhomb_[18].SetVertsID(6, 15, 20, 25);
	rhomb_[19].SetVertsID(6, 20, 11, 21);

	rhomb_[20].SetVertsID(21, 26, 16, 22);
	rhomb_[21].SetVertsID(22, 27, 17, 23);
	rhomb_[22].SetVertsID(23, 28, 18, 24);
	rhomb_[23].SetVertsID(24, 29, 19, 25);
	rhomb_[24].SetVertsID(25, 30, 20, 21);

	rhomb_[25].SetVertsID(22, 26, 27, 31);
	rhomb_[26].SetVertsID(23, 27, 28, 31);
	rhomb_[27].SetVertsID(24, 28, 29, 31);
	rhomb_[28].SetVertsID(25, 29, 30, 31);
	rhomb_[29].SetVertsID(21, 30, 26, 31);

	rhomb_[0].SetNeighFaceID(5, 4, 1, 6);
	rhomb_[1].SetNeighFaceID(6, 0, 2, 7);
	rhomb_[2].SetNeighFaceID(7, 1, 3, 8);
	rhomb_[3].SetNeighFaceID(8, 2, 4, 9);
	rhomb_[4].SetNeighFaceID(9, 3, 0, 5);

	rhomb_[5].SetNeighFaceID(10, 19, 4, 0);
	rhomb_[6].SetNeighFaceID(12, 11, 0, 1);
	rhomb_[7].SetNeighFaceID(14, 13, 1, 2);
	rhomb_[8].SetNeighFaceID(16, 15, 2, 3);
	rhomb_[9].SetNeighFaceID(18, 17, 3, 4);

	rhomb_[10].SetNeighFaceID(19, 5, 11, 20);
	rhomb_[12].SetNeighFaceID(11, 6, 13, 21);
	rhomb_[14].SetNeighFaceID(13, 7, 15, 22);
	rhomb_[16].SetNeighFaceID(15, 8, 17, 23);
	rhomb_[18].SetNeighFaceID(17, 9, 19, 24);

	rhomb_[11].SetNeighFaceID(20, 10, 6, 12);
	rhomb_[13].SetNeighFaceID(21, 12, 7, 14);
	rhomb_[15].SetNeighFaceID(22, 14, 8, 16);
	rhomb_[17].SetNeighFaceID(23, 16, 9, 18);
	rhomb_[19].SetNeighFaceID(24, 18, 5, 10);

	rhomb_[20].SetNeighFaceID(25, 29, 10, 11);
	rhomb_[21].SetNeighFaceID(26, 25, 12, 13);
	rhomb_[22].SetNeighFaceID(27, 26, 14, 15);
	rhomb_[23].SetNeighFaceID(28, 27, 16, 17);
	rhomb_[24].SetNeighFaceID(29, 28, 18, 19);

	rhomb_[25].SetNeighFaceID(29, 20, 21, 26);
	rhomb_[26].SetNeighFaceID(25, 21, 22, 27);
	rhomb_[27].SetNeighFaceID(26, 22, 23, 28);
	rhomb_[28].SetNeighFaceID(27, 23, 24, 29);
	rhomb_[29].SetNeighFaceID(28, 24, 20, 25);

	for (int i = 0; i < 60; i++)
	{
		tri_[i].SetId(i);

		int m;
		if (i % 2 == 0)
			m = i / 2;
		else
			m = (i - 1) / 2;

		int* id_set = new int[4];
		id_set = rhomb_[m].GetVertsId();

		if (i % 2 == 0)
			tri_[i].SetVertsId(id_set[0], id_set[1], id_set[2]);
		else
			tri_[i].SetVertsId(id_set[3], id_set[2], id_set[1]);
	}
}