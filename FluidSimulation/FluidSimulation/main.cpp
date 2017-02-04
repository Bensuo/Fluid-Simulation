#include <GL\glew.h>

#include <GLFW\glfw3.h>
#include <glm/glm.hpp>
#include <iostream>

using namespace std;

#define IX(x,y) ((x)+(N+2)*(y))

struct FluidCube{
	int size; 
	float dt; //Delta Time (Step-Time)
	float diff; //Diffussion Amount
	float visc; // Viscosity, thickness of the fluid

	float *s; // Source, set dye?
	float *density;

	//Velocity current
	float *Vx;
	float *Vy;

	//Velocity previous
	float *Vx0;
	float *Vy0;
};
typedef struct FluidCube FluidCube; // Ask Marco

FluidCube *FluidCubeCreate(int size, int diffusion, int viscosity, float dt){
	FluidCube *fluidCube = new FluidCube;
	int N = size;

	fluidCube->size = size;
	fluidCube->dt = dt;
	fluidCube->diff = diffusion;
	fluidCube->visc = viscosity;

	fluidCube->s = new float[(N+2) * (N + 2)];
	fluidCube->density = new float[(N + 2) * (N + 2)];

	fluidCube->Vx = new float[(N + 2) * (N + 2)];
	fluidCube->Vy = new float[(N + 2) * (N + 2)];

	fluidCube->Vx0 = new float[(N + 2) * (N + 2)];
	fluidCube->Vy0 = new float[(N + 2) * (N + 2)];

	return fluidCube;
}

// Add Density (Dye)
void fluidCubeAddDensity(FluidCube *fluidCube, int x, int y, float amount) {
	int N = fluidCube->size;
	fluidCube->density[IX(x, y)] += amount;
}

//Add Velocity
void FluidCubeAddVelocity(FluidCube *fluidCube, int x, int y, float amountX, float amountY) {
	int N = fluidCube->size;
	int index = IX(x, y);

	fluidCube->Vx[index] += amountX;
	fluidCube->Vy[index] += amountY;

}

//b= It tells the function which array it's dealing with, and so whether each boundary cell 
//should be set equal or opposite its neighbor's value.

/*
i = x
j = y
k = z

*/

static void set_bnd(int b, float *x, int N) {
	//x axis
	for (int i = 1; i < N - 1; i++) {
		x[IX(i, 0)] = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
		x[IX(i, N-1)] = b == 2 ? -x[IX(i, N-2)] : x[IX(i, N-2)];
	}

	//y axis
	for (int j = 1; j < N - 1; j++) {
		x[IX(0, j)] = b == 1 ? -x[IX(1, j)] : x[IX(1, j)];
		x[IX(N - 1, j)] = b == 1 ? -x[IX(N - 2, j)] : x[IX(N - 2, j)];
	}

	//Determines all 4 corner boundary cell velocity
	x[IX(0, 0)] = 0.5f * (x[IX(1, 0)] + x[IX(0, 1)]);
	x[IX(0, N-1)] = 0.5f * (x[IX(1, N-1)] + x[IX(0, N-2)]);

	x[IX(N - 1, N - 1)] = 0.5f * (x[IX(N - 2, N - 1)] + x[IX(N - 1, N - 2)]);
	x[IX(N - 1, 0)] = 0.5f * (x[IX(N - 2, 0)] + x[IX(N - 1, 1)]);
}

//Implementation of Gauss-Seidel relaxation formula.
//Ask Marco for clarification
void diffuse(int N, int b, float *x, float *x0, float diff, float dt) {
	int i, j, k;
	float a = dt * diff * N * N;

	for (k = 0; k < 20; k++) {
		for (i = 1; i <= N; i++) {
			for (j = 1; j <= N; j++) {
				x[IX(i, j)] = (x0[IX(i, j)] + a*(x[IX(i - 1, j)] + x[IX(i + 1, j)] + +x[IX(i, j - 1)] + +x[IX(i, j + 1)])) / (1 + 4 * a);
			}
		}
		set_bnd(b, x, N);
	}
}

// Linearly interpolate from the grid of previous density values 
// and assign this value to the current grid cell.
// Ask Marco

void advect(int N, int b, float *d, float *d0, float *u, float*v, float dt) {
	int i, j, i0, j0, i1, j1;
	float x, y, s0, t0, s1, t1, dt0;

	dt0 = dt * N;

	for (i = 1; i <= N; i++) {
		for (j = 1; j <= N; j++) {
			x = i - dt0*u[IX(i, j)]; y = j - dt0*v[IX(i, j)];
			if (x < 0.5) {
				x = 0.5;
			}

			if (x > N + 0.5) {
				x = N + 0.5;
			}

			i0 = (int)x;
			i1 = i0 + 1;

			if (y > 0.5) {
				y = 0.5;
			}

			if (y > N + 0.5) {
				y = N + 0.5;
			}

			j0 = (int)y;
			j1 = j0 + 1;

			s1 = x - i0;
			s0 = 1 - s1;
			t1 = y - j0;
			t0 = 1 - t1;

			d[IX(i, j)] = s0 * (t0 * d0[IX(i0, j0)] + t1 * d0[IX(i0, j1)]) + s1 * (t0 * d0[IX(i1, j0)] + t1 * d0[IX(i1, j1)]);
		}
	}
	set_bnd(b, d, N);
}


// It is assume that the physical length of each side of the grid is 
// one so that the grid spacing is given by h = 1 / N.

void project(int N, float *u, float *v, float *p, float *div) { //what is p?
	int i, j, k;
	float h;

	h = 1.0 / N;
	for (i = 1; i <= N; i++) {
		for (j = 1; j <= N; j++)
		{
			div[IX(i, j)] = -0.5 * h * (u[IX(i + 1, j)] - u[IX(i - 1, j)] + v[IX(i, j + 1)] - v[IX(i, j - 1)]);
			p[IX(i, j)] = 0;
		}
	}

	set_bnd(0, div, N);
	set_bnd(0, p, N);

	for (int k = 0; k < 20; k++)
	{
		for (int i = 1; i <= N; i++)
		{
			for (int j = 1; j <= N; j++)
			{
				p[IX(i, j)] = (div[IX(i, j)] + p[IX(i - 1, j)] + p[IX(i + 1, j)] + p[IX(i, j - 1)] + p[IX(i, j + 1)]) / 4;
			}
		}
		set_bnd(0, p, N);
	}

	for (int i = 1; i <= N; i++)
	{
		for (int j = 1; j <= N; j++)
		{
			u[IX(i, j)] -= 0.5 * (p[IX(i + 1, j)] - p[IX(i - 1, j)]) / h;
			v[IX(i, j)] -= 0.5 * (p[IX(i, j + 1)] - p[IX(i, j - 1)]) / h;
		}
	}
	set_bnd(1, u, N);
	set_bnd(2, v, N);
}

void FluidCubeTimeStep(FluidCube *fluidCube) {
	int N = fluidCube->size;
	float visc = fluidCube->visc;
	float diff = fluidCube->diff;
	float dt = fluidCube->dt;
	float *Vx = fluidCube->Vx;
	float *Vy = fluidCube->Vy;
	float *Vx0 = fluidCube->Vx0;
	float *Vy0 = fluidCube->Vy0;
	float *s = fluidCube->s;
	float *density = fluidCube->density;

	diffuse(N, 1, Vx, Vx0, diff, dt);
	diffuse(N, 2, Vy, Vy0, diff, dt);

	project(N, Vx0, Vy0, Vx, Vy);

	advect(N, 1, Vx, Vx0, Vx0, Vy0, dt);
	advect(N, 2, Vy, Vy0, Vx0, Vy0, dt);

	project(N, Vx, Vy, Vx0, Vy0);

	diffuse(N, 0, s, density, diff, dt);
	advect(N, 0, density, s, Vx, Vy, dt);
}

void add_source(int N, float * x, float * s, float dt) {

}

int main()
{
	return 0; 
}