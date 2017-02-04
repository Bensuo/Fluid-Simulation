#include "MacGrid.h"



MacGrid::MacGrid(int size) : dimensions(size)
{
	pressureX = new double[dimensions*dimensions];
	pressureY = new double[dimensions*dimensions];
	xVelocity = new double[dimensions*dimensions];
	xVelocity0 = new double[dimensions*dimensions];
	yVelocity = new double[dimensions*dimensions];
	yVelocity0 = new double[dimensions*dimensions];
}


MacGrid::~MacGrid()
{
}
