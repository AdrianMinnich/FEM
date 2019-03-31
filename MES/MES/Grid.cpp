#include "Grid.h"
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

void Grid::initData()
{
	fstream file;
	string data;
	file.open("data.txt", std::ios::in);
	if (file.is_open()) {
		getline(file, data);
		H = stod(data);
		getline(file, data);
		L = stod(data);
		getline(file, data);
		nH = stoi(data);
		getline(file, data);
		nL = stoi(data);
		getline(file, data);
		K = stod(data);
		getline(file, data);
		c = stod(data);
		getline(file, data);
		ro = stod(data);
		getline(file, data);
		alfa = stoi(data);
		getline(file, data);
		ambientT = stoi(data);
		getline(file, data);
		initialTemperature = stoi(data);
		getline(file, data);
		stepTime = stoi(data);
		getline(file, data);
		simTime = stoi(data);

		nN = nH * nL;
		nE = (nH - 1) * (nL - 1);
		nodes = new Node[nN];
		elements = new Element[nE];
	} 
	else {
		cout << "File not found";
	}
}

void Grid::generateGrid()
{
	initData();
	if (nN > 0 && nE > 0) {
		int n = 0;
		for (int i = 0; i < nL; i++) {
			for (int j = 0; j < nH; j++) {
				nodes[n].x = (double)(i * (L / (double)(nL - 1)));
				nodes[n].y = (double)(j * (H / (double)(nH - 1)));
				nodes[n].id = n;
				if (i == 0 || j == 0 || i == nL - 1 || j == nH - 1) // czy node lezy na krawedzi
					nodes[n].bc = 1;
				else
					nodes[n].bc = 0;
				n++;
			}
		}

		int e = 0;
		int first = 0;
		for (int i = 0; i < nL - 1; i++) {
			for (int j = 0; j < nH - 1; j++) {
				elements[e].ID[0] = first;
				elements[e].nodes[0] = nodes[first];
				elements[e].ID[1] = first + nH;
				elements[e].nodes[1] = nodes[first + nH];
				elements[e].ID[2] = first + nH + 1;
				elements[e].nodes[2] = nodes[first + nH + 1];
				elements[e].ID[3] = first + 1;
				elements[e].nodes[3] = nodes[first + 1];
				elements[e].id = e;
				e++;
				first++;
			}
			first++;
		}

		for (int i = 0; i < nE; i++) {
			if (elements[i].nodes[0].bc == 1 && elements[i].nodes[1].bc == 1) // czy 2 wezly kolo siebie w danym elemencie leza na krawedzi
				elements[i].areaIsBC[0] = 1;
			else
				elements[i].areaIsBC[0] = 0;

			if (elements[i].nodes[1].bc == 1 && elements[i].nodes[2].bc == 1)
				elements[i].areaIsBC[1] = 1;
			else
				elements[i].areaIsBC[1] = 0;

			if (elements[i].nodes[2].bc == 1 && elements[i].nodes[3].bc == 1)
				elements[i].areaIsBC[2] = 1;
			else
				elements[i].areaIsBC[2] = 0;

			if (elements[i].nodes[3].bc == 1 && elements[i].nodes[0].bc == 1)
				elements[i].areaIsBC[3] = 1;
			else
				elements[i].areaIsBC[3] = 0;
		}
	}
}

void Grid::showGrid()
{
	for (int i = 0; i < nE; i++) {
		elements[i].printElement();
	}
}

void Grid::calculate()
{
	calculateH();
	calculateC();
	calculateP();
	calculateHBC();
	calculateGlobal();
	calculateTemperature();
}

void Grid::calculateH()
{
	for (int i = 0; i < 4; i++) { // obliczanie pochodnych funkcji ksztaltu wzgledem ksi i eta dla 4pc
		dNdKSI[i][0] = -0.25 * (1.0 - eta[i]);
		dNdKSI[i][1] = 0.25 * (1.0 - eta[i]);
		dNdKSI[i][2] = 0.25 * (1.0 + eta[i]);
		dNdKSI[i][3] = -0.25 * (1.0 + eta[i]);

		dNdETA[i][0] = -0.25 * (1.0 - ksi[i]);
		dNdETA[i][1] = -0.25 * (1.0 + ksi[i]);
		dNdETA[i][2] = 0.25 * (1.0 + ksi[i]);
		dNdETA[i][3] = 0.25 * (1.0 - ksi[i]);
	}

	for (int i = 0; i < nE; i++) { // obliczanie pochodnych x/ksi, x/eta, y/ksi, y/eta w celu obliczenia jakobianu przeksztalcenia dla 4pc
		for (int j = 0; j < 4; j++) {
			elements[i].dXdKSI[j] += dNdKSI[j][0] * nodes[elements[i].ID[0]].x;
			elements[i].dXdKSI[j] += dNdKSI[j][1] * nodes[elements[i].ID[1]].x;
			elements[i].dXdKSI[j] += dNdKSI[j][2] * nodes[elements[i].ID[2]].x;
			elements[i].dXdKSI[j] += dNdKSI[j][3] * nodes[elements[i].ID[3]].x;

			elements[i].dXdETA[j] += dNdETA[j][0] * nodes[elements[i].ID[0]].x;
			elements[i].dXdETA[j] += dNdETA[j][1] * nodes[elements[i].ID[1]].x;
			elements[i].dXdETA[j] += dNdETA[j][2] * nodes[elements[i].ID[2]].x;
			elements[i].dXdETA[j] += dNdETA[j][3] * nodes[elements[i].ID[3]].x;

			elements[i].dYdKSI[j] += dNdKSI[j][0] * nodes[elements[i].ID[0]].y;
			elements[i].dYdKSI[j] += dNdKSI[j][1] * nodes[elements[i].ID[1]].y;
			elements[i].dYdKSI[j] += dNdKSI[j][2] * nodes[elements[i].ID[2]].y;
			elements[i].dYdKSI[j] += dNdKSI[j][3] * nodes[elements[i].ID[3]].y;

			elements[i].dYdETA[j] += dNdETA[j][0] * nodes[elements[i].ID[0]].y;
			elements[i].dYdETA[j] += dNdETA[j][1] * nodes[elements[i].ID[1]].y;
			elements[i].dYdETA[j] += dNdETA[j][2] * nodes[elements[i].ID[2]].y;
			elements[i].dYdETA[j] += dNdETA[j][3] * nodes[elements[i].ID[3]].y;

		}
	}


	calculateJacobian(); // przeksztalcanie elementow do ukladu lokalnego (ksi i eta), obliczane dla kazdego punktu calkowania osobno

	// tablice dla kazdego punktu calkowania
	double dNdXT1[4][4];
	double dNdXT2[4][4];
	double dNdXT3[4][4];
	double dNdXT4[4][4];
	double dNdYT1[4][4];
	double dNdYT2[4][4];
	double dNdYT3[4][4];
	double dNdYT4[4][4];
	double dH1[4][4];
	double dH2[4][4];
	double dH3[4][4];
	double dH4[4][4];

	for (int i = 0; i < nE; i++) { // obliczanie {dN/dx}*{dN/dx}T oraz {dN/dy}*{dN/dy}T
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 4; k++) {
				dNdXT1[j][k] = elements[i].dNdX[0][j] * elements[i].dNdX[0][k];
				dNdXT2[j][k] = elements[i].dNdX[1][j] * elements[i].dNdX[1][k];
				dNdXT3[j][k] = elements[i].dNdX[2][j] * elements[i].dNdX[2][k];
				dNdXT4[j][k] = elements[i].dNdX[3][j] * elements[i].dNdX[3][k];

				dNdYT1[j][k] = elements[i].dNdY[0][j] * elements[i].dNdY[0][k];
				dNdYT2[j][k] = elements[i].dNdY[1][j] * elements[i].dNdY[1][k];
				dNdYT3[j][k] = elements[i].dNdY[2][j] * elements[i].dNdY[2][k];
				dNdYT4[j][k] = elements[i].dNdY[3][j] * elements[i].dNdY[3][k];
			}
		}

		for (int j = 0; j < 4; j++) { // sumowanie ww w kazdym punkcie calkowania
			for (int k = 0; k < 4; k++) {
				dH1[j][k] = (dNdXT1[j][k] + dNdYT1[j][k]);
				dH2[j][k] = (dNdXT2[j][k] + dNdYT2[j][k]);
				dH3[j][k] = (dNdXT3[j][k] + dNdYT3[j][k]);
				dH4[j][k] = (dNdXT4[j][k] + dNdYT4[j][k]);
			}
		}

		for (int j = 0; j < 4; j++) { // mnozenie przez przewodnosc oraz wyznacznik
			for (int k = 0; k < 4; k++) {
				elements[i].H[j][k] = K * elements[i].detJ[k] * (dH1[j][k] + dH2[j][k] + dH3[j][k] + dH4[j][k]);
			}
		}

		cout << "Macierz H dla elementu " << i << endl;
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 4; k++) {
				cout << elements[i].H[j][k] << "\t";
			}
			cout << endl;
		}
		cout << endl;
	}
}

void Grid::calculateJacobian() // obliczanie detJ, detJ^-1 w celu obliczenia pochodnych funkcji ksztaltu po x i y
{
	for (int i = 0; i < nE; i++) { 
		for (int j = 0; j < 4; j++) {
			elements[i].detJ[j] = (elements[i].dXdKSI[j] * elements[i].dYdETA[j] - elements[i].dYdKSI[j] * elements[i].dXdETA[j]); // dla 4pc
			elements[i].detJ_1[j] = (1.0 / elements[i].detJ[j]);
		}

		for (int j = 0; j < 4; j++) { // dla 4pc 
			for (int k = 0; k < 4; k++) { 
				elements[i].dNdX[j][k] = elements[i].detJ_1[k] * (elements[i].dYdETA[k] * dNdKSI[j][k] - elements[i].dYdKSI[k] * dNdETA[j][k]);
				elements[i].dNdY[j][k] = elements[i].detJ_1[k] * (-(elements[i].dXdETA[k]) * dNdKSI[j][k] + elements[i].dXdKSI[k] * dNdETA[j][k]);
			}
		}
	}
}

void Grid::calculateC()
{
	double N[4][4]; // 4 funkcje ksztaltu dla 4pc
	double NN1[4][4];
	double NN2[4][4];
	double NN3[4][4];
	double NN4[4][4];
	double tmp;

	for (int i = 0; i < 4; i++) { // obliczanie funkcji ksztaltu dla 4 punktow calkowania
		N[i][0] = 0.25 * (1.0  -ksi[i]) * (1.0 - eta[i]);
		N[i][1] = 0.25 * (1.0 + ksi[i]) * (1.0 - eta[i]);
		N[i][2] = 0.25 * (1.0 + ksi[i]) * (1.0 + eta[i]);
		N[i][3] = 0.25 * (1.0  -ksi[i]) * (1.0 + eta[i]);
	}

	for (int i = 0; i < nE; i++) { // obliczanie {N}*{N}T oraz mnozenie c (cieplo wlasciwe) i ro (gestosc)
		for (int j = 0; j < 4; j++) { // dla 4 pc
			for (int k = 0; k < 4; k++) {
				tmp = elements[i].detJ[j] * c * ro; 
				NN1[j][k] = N[0][j] * N[0][k] * tmp; // macierze w 4 punktach calkowania
				NN2[j][k] = N[1][j] * N[1][k] * tmp;
				NN3[j][k] = N[2][j] * N[2][k] * tmp;
				NN4[j][k] = N[3][j] * N[3][k] * tmp;
			}
		}

		for (int j = 0; j < 4; j++) { // sumowanie macierzy z 4 punktow calkowania
			for (int k = 0; k < 4; k++) {
				elements[i].C[j][k] = NN1[j][k] + NN2[j][k] + NN3[j][k] + NN4[j][k];
			}
		}

		cout << "Macierz C dla elementu " << i << endl;
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 4; k++) {
				cout << elements[i].C[j][k] << "\t";
			}
			cout << endl;
		}
		cout << endl;

	}
}

void Grid::calculateP() // mnozone razy 2 po 2 punkty calkowania (suma wartosci funkcji ksztaltu)
{
	for (int i = 0; i < nE; i++) {
		double x0 = elements[0].nodes[0].x;
		double y0 = elements[0].nodes[0].y;
		double x1 = elements[0].nodes[1].x;
		double y1 = elements[0].nodes[1].y;
		double x2 = elements[0].nodes[2].x;
		double y2 = elements[0].nodes[2].y;
		double x3 = elements[0].nodes[3].x;
		double y3 = elements[0].nodes[3].y;

		double L[4];
		L[0] = sqrt(pow(x1 - x0, 2) + pow(y1 - y0, 2));
		L[1] = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
		L[2] = sqrt(pow(x3 - x2, 2) + pow(y3 - y2, 2));
		L[3] = sqrt(pow(x0 - x3, 2) + pow(y0 - y3, 2));

		for (int j = 0; j < 4; j++) { // obliczanie wektora {P} = alfa * {N} * ambientT
			if (elements[i].areaIsBC[j]) // czy bok lezy na krawedzi
				elements[i].P[j] = 2 * alfa * (L[j] / 2) * ambientT; // suma funkcji ksztaltu w 1pc = 1, wiec w 2pc = 2, (L[j] / 2) to detJ (stosunek dlugosci boku globalnego do lokalnego)
			else
				elements[i].P[j] = 0;
		}
	}
}

void Grid::calculateHBC() // zastosowanie metode z obliczaniem 2 odpowiednich funkcji ksztaltu zamiast 4 (jakobian 1D)
{
	// N*N*alfa dla 2 pc
	double areaPc1[2][2];
	double areaPc2[2][2];

	double area1[2][2];
	double area2[2][2];
	double area3[2][2];
	double area4[2][2];

	double x0 = elements[0].nodes[0].x;
	double y0 = elements[0].nodes[0].y;
	double x1 = elements[0].nodes[1].x;
	double y1 = elements[0].nodes[1].y;
	double x2 = elements[0].nodes[2].x;
	double y2 = elements[0].nodes[2].y;
	double x3 = elements[0].nodes[3].x;
	double y3 = elements[0].nodes[3].y;

	double L[4]; // dlugosci bokow
	L[0] = sqrt(pow(x1 - x0, 2) + pow(y1 - y0, 2));
	L[1] = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
	L[2] = sqrt(pow(x3 - x2, 2) + pow(y3 - y2, 2));
	L[3] = sqrt(pow(x0 - x3, 2) + pow(y0 - y3, 2));

	double tmpN[2][2]; // 2 funkcje ksztaltu dla 2 punktow calkowania po powierzchni
	for (int i = 0; i < 2; i++) {
		tmpN[i][0] = 0.5 * (1 - ksi[i]);
		tmpN[i][1] = 0.5 * (1 + ksi[i]);
	}

	for (int i = 0; i < nE; i++) { // dla wszystkich elementow
		for (int j = 0; j < 2; j++) { // N*N*alfa dla 2 punktow calkowania po powierzchni
			for (int k = 0; k < 2; k++) {
				areaPc1[j][k] = tmpN[0][j] * tmpN[0][k] * alfa;
				areaPc2[j][k] = tmpN[1][j] * tmpN[1][k] * alfa;
			}
		}

		for (int j = 0; j < 2; j++) { // sprawdzanie na ktore krawedzie nalozono warunek brzegowy (jezeli nalozono to *1 jezeli nie to *0)
			for (int k = 0; k < 2; k++) {
				area1[j][k] = (areaPc1[j][k] + areaPc2[j][k]) * elements[i].areaIsBC[0];
				area2[j][k] = (areaPc1[j][k] + areaPc2[j][k]) * elements[i].areaIsBC[1];
				area3[j][k] = (areaPc1[j][k] + areaPc2[j][k]) * elements[i].areaIsBC[2];
				area4[j][k] = (areaPc1[j][k] + areaPc2[j][k]) * elements[i].areaIsBC[3];
			}
		}

		for (int j = 0; j < 2; j++) { // mnozenie przez detJ (deltaX/2), stosunek dlugosci boku globalnego do lokalnego, wsadzanie wartosci z macierzy powierzchni 2x2 do macierzy elementu 4x4
			for (int k = 0; k < 2; k++) {
				elements[i].H_BC[j][k] += area1[j][k] * (L[0] / 2); // (0,0) (0,1) (1,0) (1,1)
				elements[i].H_BC[j + 1][k + 1] += area2[j][k] * (L[1] / 2); // (1,1) (1,2) (2,1) (2,2)
				elements[i].H_BC[j + 2][k + 2] += area3[j][k] * (L[2] / 2); // (2,2) (2,3) (3,2) (3,3)
				elements[i].H_BC[j * 3][k * 3] += area4[j][k] * (L[3] / 2); // (0,0) (0,3) (3,0) (3,3)
			}
		}

		for (int j = 0; j < 4; j++) { // dodawanie warunku brzegowego do macierzy H w kazdym elemencie
			for (int k = 0; k < 4; k++) {
				elements[i].H[j][k] += elements[i].H_BC[j][k];
			}
		}

		cout << "Macierz H z warunkiem brzegowym dla elementu " << i << endl;
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 4; k++) {
				cout << elements[i].H[j][k] << "\t";
			}
			cout << endl;
		}
		cout << endl;
	}
}

void Grid::calculateGlobal()
{
	globalH = new double*[nN];
	globalC = new double*[nN];
	globalP = new double[nN];

	for (int i = 0; i < nN; i++) {
		globalP[i] = 0;
		globalH[i] = new double[nN];
		globalC[i] = new double[nN];
		for (int j = 0; j < nN; j++) { // macierz jest wielkosci nN x nN
			globalH[i][j] = 0;
			globalC[i][j] = 0;
		}
	}
	for (int i = 0; i < nE; i++) { 
		int elementID[4];
		elementID[0] = elements[i].ID[0];
		elementID[1] = elements[i].ID[1];
		elementID[2] = elements[i].ID[2];
		elementID[3] = elements[i].ID[3];
		for (int j = 0; j < 4; j++) {
			globalP[elementID[j]] += elements[i].P[j];
			for (int k = 0; k < 4; k++) {
				globalH[elementID[j]][elementID[k]] += elements[i].H[j][k];
				globalC[elementID[j]][elementID[k]] += elements[i].C[j][k];
			}
		}
	}
	cout << "Globalna macierz H: " << endl;
	for (int i = 0; i < nN; i++) {
		for (int j = 0; j < nN; j++) {
			cout << globalH[i][j] << "   ";
		}
		cout << endl;
	}
	
	cout << endl;

	cout << "Globalna macierz C: " << endl;
	for (int i = 0; i < nN; i++) {
		for (int j = 0; j < nN  ; j++) {
			cout << globalC[i][j] << "     ";
		}
		cout << endl;
	}
}

void Grid::calculateTemperature()
{
	double **H_final;	
	double *p_final;
	double **HP;
	double *t1;
	double min, max;
	double MIN[10], MAX[10];
	
	H_final = new double*[nN];
	p_final = new double[nN];
	t1 = new double[nN];
	HP = new double*[nN];

	for (int i = 0; i < nN; i++)
		HP[i] = new double[nN + 1];

	for (int i = 0; i < nN; i++) {
		t1[i] = initialTemperature; // wektor temperatur poczatkowych
	}

	for (int i = 0; i < nN; i++) {
		for (int j = 0; j < nN; j++) {
			globalC[i][j] = globalC[i][j] / stepTime; // [C] = [C]/dT
		}
	}

	for (int i = 0; i < nN; i++) {
		H_final[i] = new double[nN];
		for (int j = 0; j < nN; j++) {
			H_final[i][j] = globalH[i][j] + globalC[i][j]; // [H] = [H] + [C]/dT
		}
	}

	cout << endl;

	for (int s = stepTime; s <= simTime; s+=stepTime) { // petla o krok czasowy
	
		p_final = calculatePFinal(globalP, globalC, t1); // {P} = {P} + ([C]/dT) * t0

		for (int i = 0; i < nN; i++) {
			for (int j = 0; j < nN; j++) {
				HP[i][j] = H_final[i][j];
			}
		}

		for (int i = 0; i < nN; i++) {
			HP[i][nN] = p_final[i]; // macierz HP to macierz H z dopisanym wektorem P w ostatniej kolumnie   A x X = B    => A|B - macierz rozszerzona, X - wektor nieznanych temperatur
		}

		t1 = gauss(HP, t1, nN); // obliczanie wektora temperatur za pomoca metody eliminacji gaussa
	
	
		for (int i = 0; i < nN; i++) { // wyszukiwanie min i max temp
			if (i == 0) {
				min = t1[i];
				max = t1[i];
			}
			else {
				if (t1[i] < min)	min = t1[i];
				if (t1[i] > max)	max = t1[i];
			}
		}
		MIN[s / 50 - 1] = min;
		MAX[s / 50 - 1] = max;

		cout << "ITERATION " << (s / 50) - 1 << endl;
		cout << "Matrix [H]+[C]/dT" << endl;
		for (int i = 0; i < nN; i++) {
			for (int j = 0; j < nN; j++) {
				cout << H_final[i][j] << "    ";
			}
			cout << endl;
		}

		cout << endl;

		cout << "VECTOR {P}+{[C]/dT}*{T0}" << endl;
		for (int i = 0; i < nN; i++) {
			cout << p_final[i] << "   ";
		}
		cout << endl;

		cout << "MIN: " << min << ", MAX: " << max << endl << endl << endl;

		min = 0;
		max = 0;
		
		for (int i = 0; i < nN; i++)
			p_final[i] = 0;
	}
	cout << "Time[s]\tMinTemp[s]\tMaxTemp[s]\n";
	for (int s = stepTime; s <= simTime; s += stepTime) {
		cout << s << "\t" << MIN[s / 50 - 1] << "\t\t" << MAX[s / 50 - 1] << endl;
	}
}

double * Grid::calculatePFinal(double * p, double ** C, double * t0)
{
	double *p_final;
	double *c_temp;

	p_final = new double[nN];
	c_temp = new double[nN];

	for (int i = 0; i < nN; i++) {
		for (int j = 0; j < nN; j++) {
			c_temp[i] = 0;
		}
	}

	for (int i = 0; i < nN; i++) {
		for (int j = 0; j < nN; j++) {
			c_temp[i] += C[i][j] * t0[j];
		}
	}

	for (int i = 0; i < nN; i++) {
		p_final[i] = p[i] + c_temp[i];
	}

	return p_final;
}

double * Grid::gauss(double ** HP, double * T, int n) // metoda eliminacji Gaussa
{
	const double eps = 1e-12; // dokladnosc porownania z zerem
	int i, j, k;
	double m, s; 
	for (i = 0; i < n - 1; i++) { // eliminacja wspó³czynników - przeksztalcanie macierzy w macierz trojkatna gorna
		for (j = i + 1; j < n; j++) {
			if (fabs(HP[i][i]) < eps) return NULL;
			m = -HP[j][i] / HP[i][i];
			for (k = i + 1; k <= n; k++)
				HP[j][k] += m * HP[i][k];
		}
	}
	for (i = n - 1; i >= 0; i--) { // wyliczanie niewiadomych
		s = HP[i][n];
		for (j = n - 1; j >= i + 1; j--)
			s -= HP[i][j] * T[j];
		if (fabs(HP[i][i]) < eps) return NULL;
		T[i] = s / HP[i][i];
	}
	return T;
}

Grid::Grid()
{
	ksi[0] = -(1.0 / sqrt(3));
	ksi[1] = -ksi[0];
	ksi[2] = ksi[1];
	ksi[3] = -ksi[2];
	eta[0] = ksi[0];
	eta[1] = eta[0];
	eta[2] = -eta[1];
	eta[3] = eta[2];
}

Grid::~Grid()
{
}
