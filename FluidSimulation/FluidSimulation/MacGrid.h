#pragma once
#include <glm/glm.hpp>
class MacGrid
{
public:
	MacGrid(int size);
	virtual ~MacGrid();
private:
	int dimensions;
	double *pressureX, *pressureY;
	double *xVelocity, *xVelocity0;
	double *yVelocity, *yVelocity0;
	
	
};

