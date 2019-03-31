#pragma once
class Node
{
public:
	double x;
	double y;
	double t0;
	int id;
	int bc;
	Node();
	Node(double x, double y);
	Node(double x, double y, double t0);
	~Node();
	void showNode();
};

