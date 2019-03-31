#include "Node.h"
#include <iostream>
using namespace std;

Node::Node()
{
}

Node::Node(double x, double y)
{
	this->x = x;
	this->y = y;
}

Node::Node(double x, double y, double t0)
{
	this->x = x;
	this->y = y;
	this->t0 = t0;
}

Node::~Node()
{
}

void Node::showNode()
{
	cout << id << ".(" << x << ", " << y << ")" << endl;
}
