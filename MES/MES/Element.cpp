#include "Element.h"
#include <iostream>
using namespace std;

Element::Element()
{
	nodes = new Node[4];
	for (int i = 0; i < 4; i++) {
		ID[i] = 0;
		dXdKSI[i] = 0;
		dXdETA[i] = 0;
		dYdKSI[i] = 0;
		dYdETA[i] = 0;
		detJ[i] = 0;
		detJ_1[i] = 0;
		for (int j = 0; j < 4; j++) {
			dNdX[i][j] = 0;
			dNdY[i][j] = 0;
			H_BC[i][j] = 0;
		}
	}
}

Element::~Element()
{
}

void Element::printElement()
{
	cout << "Element " << id << endl;
	for (int i = 0; i < 4; i++) {
		nodes[i].showNode();
	}
	cout << endl;
}
