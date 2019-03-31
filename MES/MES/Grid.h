#pragma once
#include "Element.h"
#include "Node.h"
class Grid
{
public:
	double H;
	double L;
	int nH;
	int nL;
	int nN;
	int nE;
	double K; // conductivity
	double c; // specific heat
	double ro; // density
	double alfa; // convection
	double ambientT;
	double initialTemperature; // initial temperature
	double stepTime; // simulation step Time
	double simTime; // simulation time
	Node *nodes;
	Element *elements;
	double ksi[4];
	double eta[4];
	double dNdKSI[4][4];
	double dNdETA[4][4];
	double **globalH;
	double **globalC;
	double *globalP;
	void initData();
	void generateGrid();
	void showGrid();
	void calculate();
	void calculateH();
	void calculateJacobian();
	void calculateC();
	void calculateP();
	void calculateHBC();
	void calculateGlobal();
	void calculateTemperature();
	double* calculatePFinal(double *p, double **C, double *t0);
	double* gauss(double **A, double *B, int n);
	Grid();
	~Grid();
};

