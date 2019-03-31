#include <iostream>
#include "Grid.h"

using namespace std;

int main() {
	Grid *grid = new Grid();
	grid->generateGrid();
	//grid->showGrid();
	grid->calculate();

	system("PAUSE");
}