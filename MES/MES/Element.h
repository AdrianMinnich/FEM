#pragma once
#include "Node.h"
class Element
{
public:
	int id;
	Node *nodes;
	int ID[4];
	double dXdKSI[4];
	double dXdETA[4];
	double dYdKSI[4];
	double dYdETA[4];
	double dNdX[4][4];
	double dNdY[4][4];
	double detJ[4];
	double detJ_1[4];
	double H[4][4];
	double C[4][4];
	double P[4];
	int areaIsBC[4];
	double H_BC[4][4];
	Element();
	~Element();
	void printElement();
};

