#include <glm/glm.hpp>
#include <SDL.h>
#include <math.h>
#include <GL\glew.h>
#undef main
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
using namespace std;

#define IX(x,y) ((x)+(N+2)*(y)) 

#define CELL_SIZE 10 //Distance between each vertex point building one square cell.
#define GRID_SIZE 100 //Number of cells in the grid. 
#define DISPLAY_SIZE_X 1280 // Size of display on x cords.
#define DISPLAY_SIZE_Y 720 // Size of display on y cords. 

//Global Variables

int x, y;

//Mouse position
int mousePosiX;
int mousePosiY;
int mousePosiZ;

struct Particle
{
	float x = 0;
	float y = 0;
	bool isActive = false;
};

enum CellState
{
	EMPTY,
	FLUID,
	WALL
};

struct FluidCube {
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

	Particle *particles;
	CellState *cellStates;
};
typedef struct FluidCube FluidCube;

//Set aspect ratio (world cords)
const static float ASPECT_RATIO = ((float)DISPLAY_SIZE_X / (float)DISPLAY_SIZE_Y);
const static float CELL_SIZE_X = (float)((1.0f / ((float)GRID_SIZE / (float)DISPLAY_SIZE_X))); //size of cells in grid
const static float CELL_SIZE_Y = (float)((1.0f / ((float)GRID_SIZE / (float)DISPLAY_SIZE_Y))); //size of cells in grid

FluidCube *FluidCubeCreate(int size, float diffusion, float viscosity, float dt) {
	FluidCube *fluidCube = new FluidCube;
	int N = size;

	fluidCube->size = size;
	fluidCube->dt = dt;
	fluidCube->diff = diffusion;
	fluidCube->visc = viscosity;

	fluidCube->s = new float[(N + 2) * (N + 2)]();
	fluidCube->density = new float[(N + 2) * (N + 2)]();

	fluidCube->Vx = new float[(N + 2) * (N + 2)]();
	fluidCube->Vy = new float[(N + 2) * (N + 2)]();

	fluidCube->Vx0 = new float[(N + 2) * (N + 2)]();
	fluidCube->Vy0 = new float[(N + 2) * (N + 2)]();

	//Bens stuff
	fluidCube->cellStates = new CellState[(N + 2) * (N + 2)]();
	fluidCube->particles = new Particle[N*N*2]();

	return fluidCube;
}


// Add Density (Dye)
void fluidCubeAddDensity(FluidCube *fluidCube, int x, int y, float amount) {
	int N = fluidCube->size;
	fluidCube->s[IX(x, y)] += amount;
	//cout << "Density amount at point: " << (float)fluidCube->s[IX(x, y)] << endl;
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
NOTE TO SELF
i = x
j = y
k = z

b = ?
x = float based on velocity reference so it can be changed
N = size of grid
*/

static void set_bnd(int b, float *x, int N) {


	int i;
	for (i = 1; i <= N; i++) {
		// set points for the 4 edges of the boundry based on the size of grid for next step based on current step
		// bounces velocity by inversing if boundry hits the boundry
		x[IX(0, i)] = b == 1 ? -x[IX(1, i)] : x[IX(1, i)];
		x[IX(N + 1, i)] = b == 1 ? -x[IX(N, i)] : x[IX(N, i)];
		x[IX(i, 0)] = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
		x[IX(i, N + 1)] = b == 2 ? -x[IX(i, N)] : x[IX(i, N)];
	}
	/*
	set_bnd(0, div, N);
	set_bnd(0, p, N);

	set points for the 4 corners of the boundry based on size of grid for the current step
	*/

	x[IX(0, 0)] = 0.5*(x[IX(1, 0)] + x[IX(0, 1)]);
	x[IX(0, N + 1)] = 0.5*(x[IX(1, N + 1)] + x[IX(0, N)]);
	x[IX(N + 1, 0)] = 0.5*(x[IX(N, 0)] + x[IX(N + 1, 1)]);
	x[IX(N + 1, N + 1)] = 0.5*(x[IX(N, N + 1)] + x[IX(N + 1, N)]);
}


/*
Function finds previous velocity

Implementation of Gauss-Seidel relaxation formula.
diffuse(N, 1, Vx0, Vx, visc, dt);

N = Size
b = checks if previous velocity x or y as 1 or 2
x = previous velocity
x0 = current velocity
diff = viscocity
dt = deltaTime / Timestep

diffuse(N, 0, s, density, diff, dt); //computes diffusion for next step in grid
*/
void diffuse(int N, int b, float *x, float *x0, float diff, float dt) {
	int i, j, k;
	float a = dt * diff * N * N; // Delta Time * Viscocity * Size * Size

	for (k = 0; k < 20; k++) {
		for (i = 1; i <= N; i++) {
			for (j = 1; j <= N; j++) {
				x[IX(i, j)] = (x0[IX(i, j)] + a*(x[IX(i - 1, j)] + x[IX(i + 1, j)] + +x[IX(i, j - 1)] + +x[IX(i, j + 1)])) / (1 + 4 * a);
			}
		}
		set_bnd(b, x, N);
	}
}

/*	
Linearly interpolate from the grid of previous density values
and assign this value to the current grid cell.

advect(N, 1, Vx, Vx0, Vx0, Vy0, dt);
advect(N, 2, Vy, Vy0, Vx0, Vy0, dt);

N = size
b = if x or y
d* = current velocity of b
d0* = previous velocity of b
u* = previous velocity of x
v* = previous velocity of y
dt = delta time

advection definition: the transfer of heat or matter by the flow of a fluid,
especially horizontally in the atmosphere or the sea.

advect(N, 0, density, s, Vx, Vy, dt); //computes advection of current density for next step.

*/


void advect(int N, int b, float *d, float *d0, float *u, float*v, float dt) {
	int i, j, i0, j0, i1, j1;
	float x, y, s0, t0, s1, t1, dt0;

	dt0 = dt * N; // delta time * size

	/*
	
	*/
	for (i = 1; i <= N; i++) {
		for (j = 1; j <= N; j++) {
			x = i - dt0*u[IX(i, j)]; y = j - dt0*v[IX(i, j)];
			if (x < 0.5) { 
				x = 0.5;
			}
			/*
			-----above-----
			clamp the value of x between 0.5 and N + 0.5
			-----below-----
			*/
			if (x > N + 0.5) {
				x = N + 0.5;
			}

			i0 = (int)x; //if x = 1.5, x = 1, wrap to int.
			i1 = i0 + 1;

			if (y < 0.5) {
				y = 0.5;
			}
			/*
			---- - above---- -
				clamp the value between 0.5 and N + 0.5
			---- - below---- -
			*/
			if (y > N + 0.5) {
				y = N + 0.5;
			}

			j0 = (int)y; //if y = 1.5, y = 1, wrap to int.
			j1 = j0 + 1;

			// math.h round function <---- similar----below for x and y position.
			s1 = x - i0;  
			s0 = 1 - s1;
			t1 = y - j0;
			t0 = 1 - t1;

			//computes current velocity of x or y
			d[IX(i, j)] = s0 * (t0 * d0[IX(i0, j0)] + t1 * d0[IX(i0, j1)]) + s1 * (t0 * d0[IX(i1, j0)] + t1 * d0[IX(i1, j1)]);
		}
	}
	set_bnd(b, d, N); // ensures advection is not going outside of the boundry, if not fix it.
}


/*
project(N, Vx0, Vy0, Vx, Vy);

It is assume that the physical length of each side of the grid is
one so that the grid spacing is given by h = 1 / N.

N = size
U = array previous velocity X
v = array previous velocity Y
p = array of current velocity X / possibly pressure as well
div = array  of current velocity Y
h = square root of cell (half a cell)

*/


void project(int N, float *u, float *v, float *p, float *div) { 
	int i, j, k;
	float h;

	h = 1.0 / N; 
	for (i = 1; i <= N; i++) {
		for (j = 1; j <= N; j++)
		{// takes all 4 values to...    
			/*
			-0.5 * h * (left cell X - right cell X + top cell Y - bottom cell Y)
			compute and sets array of current velocity for Y 
			
			possibly sets magnitude and angle - > computes vector?
			*/
			div[IX(i, j)] = -0.5 * h * (u[IX(i + 1, j)] - u[IX(i - 1, j)] + v[IX(i, j + 1)] - v[IX(i, j - 1)]); 
			p[IX(i, j)] = 0; //sets current velocity of X to 0
		}
	}

	set_bnd(0, div, N);
	set_bnd(0, p, N);

	// why 20?
	for (int k = 0; k < 20; k++)
	{
		for (int i = 1; i <= N; i++)
		{
			for (int j = 1; j <= N; j++)
			{
				// why are we doing same calc 20x? We reset I and J to 1, 20x.
				// compute and sets array of current velocity for X
				//P value changes as we progress through computation.
				p[IX(i, j)] = (div[IX(i, j)] + p[IX(i - 1, j)] + p[IX(i + 1, j)] + p[IX(i, j - 1)] + p[IX(i, j + 1)]) / 4;
			}
		}
		set_bnd(0, p, N);
	}

	for (int i = 1; i <= N; i++)
	{
		for (int j = 1; j <= N; j++)
		{	// compute and sets array of previous velocity for X and Y based on current velocity X
			u[IX(i, j)] -= 0.5 * (p[IX(i + 1, j)] - p[IX(i - 1, j)]) / h;
			v[IX(i, j)] -= 0.5 * (p[IX(i, j + 1)] - p[IX(i, j - 1)]) / h;
		}
	}
	set_bnd(1, u, N);
	set_bnd(2, v, N);
}

//add_source(N, Vx, Vx0, dt);
/*
N = size
Vx = velocity
Vx0 = previous velocity
dt = time step

add_source(N, density, s, dt); //compute density into grid

*/
void add_source(int N, float * Vx, float * Vx0, float dt) {
	int i, size = (N + 2) * (N + 2);
	for (i = 0; i < size; i++)
	{
		Vx[i] += dt*Vx0[i];
	}

}

void getForces(int N, float * s, float * Vx0, float * Vy0) //reset Force to 0
{
	int i, size = (N + 2) * (N + 2);
	for (i = 0; i < size; i++)
	{
		s[i] = 0;
		Vx0[i] = 0;
		Vy0[i] = 0.0f;
	}
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
	float *s = fluidCube->s; // density of previous step
	float *density = fluidCube->density; // density of current step
	Particle* particles = fluidCube->particles;
	CellState *cellStates = fluidCube->cellStates;


	//Velocity step
	add_source(N, Vx, Vx0, dt); //compute velocity for X
	add_source(N, Vy, Vy0, dt); //compute velocity for Y
	diffuse(N, 1, Vx0, Vx, visc, dt); //computes current velocity for next step for X
	diffuse(N, 2, Vy0, Vy, visc, dt); //computes current velocity for next step for Y

	project(N, Vx0, Vy0, Vx, Vy);  //Set boundry as a square/cube and computes reaction of velocity if the boundry is hit. 

	advect(N, 1, Vx, Vx0, Vx0, Vy0, dt); //computes advection of current velocity for next step for X
	advect(N, 2, Vy, Vy0, Vx0, Vy0, dt); //computes advection of current velocity for next step for Y

	project(N, Vx, Vy, Vx0, Vy0); //Set boundry as a square/cube and computes reaction of velocity if the boundry is hit.

	//Density step
	add_source(N, density, s, dt); //compute density into grid
	diffuse(N, 0, s, density, diff, dt); //computes diffusion for next step in grid
	advect(N, 0, density, s, Vx, Vy, dt); //computes advection of current density for next step. 

	getForces(N, s, Vx0, Vy0); //Do all of above and then reset velocity and density to 0 for x and y
}

void drawFluidVelocity(FluidCube* fluidCube) {
	int i, j;
	int N = fluidCube->size;
	//cout << "N " << N << endl;
	float x, y;
	float h = 1.0f;
	float *Vx = fluidCube->Vx;
	float *Vy = fluidCube->Vy;
	float cellSize = CELL_SIZE;

	glLineWidth(2.0f);

	int velocityMultiplier = 10;


	glBegin(GL_LINES);
	
	for (i = 1; i <= N; i++) {
		x = (i) * h;
		for (j = 1; j <= N; j++) {
			y = (j) * h;

			int pointVeloX = Vx[IX(i, j)];
			int pointVeloY = Vy[IX(i, j)];

			int actualX = x - Vx[IX(i, j)]; 
			int actualY = y - Vy[IX(i, j)];

			
			int addPointVelo = pointVeloX + pointVeloY;

			glColor3f(1.0f, 1.0f, 1.0f);
			glVertex2f(x*CELL_SIZE_X, y*CELL_SIZE_Y);
			glColor3f(1.0, 0.0, (addPointVelo / 10.0));
			glVertex2f((x * CELL_SIZE_X + (Vx[IX(i, j)] * velocityMultiplier)) , (y * CELL_SIZE_Y + (Vy[IX(i, j)] * velocityMultiplier)) );
		}
	}
	glEnd();

}

void drawFluidDensity(FluidCube* fluidCube) {

	int N = fluidCube->size;
	int size = N*N;
	float h = 1.0f;
	float sourceAlpha = 0.5;
	float sourceAlpha00 = 0.05;
	float sourceAlpha01 = 0.05;
	float sourceAlpha10 = 0.05;
	float sourceAlpha11 = 0.05;
	float d00, d01, d10, d11;
	float d00R, d01R, d10R, d11R;
	float d00G, d01G, d10G, d11G;
	float d00B, d01B, d10B, d11B;

	for (int i = 0; i < N; i++)
	{
		float x = (i )*h;//- 0.5f
		for (int j = 0; j < N; j++) {
			float y = (j )*h;//- 0.5f
			//float density = fluidCube->density[i+j];

			//calculate position of the fluid in the grid
			d00 = fluidCube->density[IX(i, j)]; 
			d01 = fluidCube->density[IX(i, j + 1 )];
			d11 = fluidCube->density[IX(i + 1, j + 1)];
			d10 = fluidCube->density[IX(i + 1, j)];



			//Add Density color and hue if greater than 0.1
			if (fluidCube->density[IX(i, j)] > 0.1) {
				d00R = d00 + 0.0f;
				d00G = d00 + 1.0f;
				d00B = d00 + 0.0f;
				sourceAlpha00 = 1.0;
			}
			else {
				d00R = d00 + 1.0f;
				d00G = d00 + 1.0f;
				d00B = d00 + 1.0f;
				sourceAlpha00 = 1.0;
			}

			if (fluidCube->density[IX(i, j + 1)] > 0.1) {
				d01R = d01 + 0.0f;
				d01G = d01 + 1.0f;
				d01B = d01 + 0.0f;
				sourceAlpha01 = 1.0;
			}
			else {
				d01R = d01 + 1.0f;
				d01G = d01 + 1.0f;
				d01B = d01 + 1.0f;
				sourceAlpha01 = 1.0;
			}

			if (fluidCube->density[IX(i + 1, j + 1)] > 0.1) {
				d11R = d11 + 0.0f;
				d11G = d11 + 1.0f;
				d11B = d11 + 0.0f;
				sourceAlpha11 = 1.0;
			}
			else {
				d11R = d11 + 1.0f;
				d11G = d11 + 1.0f;
				d11B = d11 + 1.0f;
				sourceAlpha11 = 1.0;
			}

			if (fluidCube->density[IX(i + 1, j)] > 0.1) {
				d10R = d10 + 0.0f;
				d10G = d10 + 1.0f;
				d10B = d10 + 0.0f;
				sourceAlpha10 = 1.0;
			}
			else {
				d10R = d10 + 1.0f;
				d10G = d10 + 1.0f;
				d10B = d10 + 1.0f;
				sourceAlpha10 = 1.0;
			}

			//Add Density color and hue if greater than 0.2
			if (fluidCube->density[IX(i, j)] >= 0.5) {
				d00R = d00 + 0.0f;
				d00G = d00 + 0.0f;
				d00B = d00 + 1.0f;
				sourceAlpha00 = 1.0;
			}
			else {
				d00R = d00 + 0.0f;
				d00G = d00 + 0.0f;
				d00B = d00 + 0.0f;
				sourceAlpha00 = 1.0;
			}


			if (fluidCube->density[IX(i, j + 1)] >= 1.0) {
				d01R = d01 + 0.0f;
				d01G = d01 + 0.0f;
				d01B = d01 + 1.0f;
				sourceAlpha01 = 1.0;
			}
			else {
				d01R = d01 + 0.0f;
				d01G = d01 + 0.0f;
				d01B = d01 + 0.0f;
				sourceAlpha01 = 1.0;
			}

			if (fluidCube->density[IX(i + 1, j + 1)] >= 1.0) {
				d11R = d11 + 0.0f;
				d11G = d11 + 0.0f;
				d11B = d11 + 1.0f;
				sourceAlpha11 = 1.0;
			}
			else {
				d11R = d11 + 0.0f;
				d11G = d11 + 0.0f;
				d11B = d11 + 0.0f;
				sourceAlpha11 = 1.0;
			}

			if (fluidCube->density[IX(i + 1, j)] >= 1.0) {
				d10R = d10 + 0.0f;
				d10G = d10 + 0.0f;
				d10B = d10 + 1.0f;
				sourceAlpha10 = 1.0;
			}
			else {
				d10R = d10 + 0.0f;
				d10G = d10 + 0.0f;
				d10B = d10 + 0.0f;
				sourceAlpha10 = 1.0;
			}

			// draw density as a cube of quads

			glBegin(GL_QUADS);
			glColor4f(d00R, d00G, d00B, sourceAlpha00);
			glVertex3f(x * CELL_SIZE_X, y * CELL_SIZE_Y, 0);

			glColor4f(d01R, d01G, d01B, sourceAlpha01);
			glVertex3f(x * CELL_SIZE_X, (y + h) * CELL_SIZE_Y, 0);

			glColor4f(d11R, d11G, d11B, sourceAlpha11);
			glVertex3f((x + h) * CELL_SIZE_X, (y + h) * CELL_SIZE_Y, 0);

			glColor4f(d10R, d10G, d10B, sourceAlpha10);
			glVertex3f((x + h) * CELL_SIZE_X, y * CELL_SIZE_Y, 0);

			glEnd();

		}
	}


}

//Window settings
const int SCREEN_WIDTH = DISPLAY_SIZE_X;
const int SCREEN_HEIGHT = DISPLAY_SIZE_Y;
SDL_Window* gWindow = NULL;
SDL_Renderer* gRenderer = NULL;
//TTF_Font* gFont = NULL;

int main()
{
	SDL_Init(SDL_INIT_VIDEO);
	//TTF_Init();

	// Request an OpenGL 3.0 context.
	SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
	SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 0);
	SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);

	SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);  // double buffering on
	SDL_GL_SetAttribute(SDL_GL_MULTISAMPLEBUFFERS, 1);
	SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES, 4); // Turn on x4 multisampling anti-aliasing (MSAA)

	gWindow = SDL_CreateWindow("Fluid Simulator - IPM", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, SCREEN_WIDTH, SCREEN_HEIGHT, SDL_WINDOW_OPENGL | SDL_WINDOW_SHOWN);
	SDL_GLContext context = SDL_GL_CreateContext(gWindow); //Create OpenGL Context and binds the gWindow to SDL.
	gRenderer = SDL_CreateRenderer(gWindow, -1, SDL_RENDERER_ACCELERATED);
	glEnable(GL_BLEND);
	FluidCube * fluidCube = FluidCubeCreate(GRID_SIZE, 0.00001f, 0.1f, 0.1f);
	float totalFluid = 5000;
	float currentFluid = totalFluid;
	bool running = true; // set running to true

	SDL_Event sdlEvent;  // variable to detect SDL events
	while (running) {	// the event loop
		while (SDL_PollEvent(&sdlEvent)) {
			if (sdlEvent.type == SDL_QUIT)
				running = false;
		}
		const Uint8 *keys = SDL_GetKeyboardState(NULL);
		
		//Release Dam
		if (running)
		{
			if (currentFluid >= totalFluid) {
				for (int j = 1; j <= GRID_SIZE; j++) {
					FluidCubeAddVelocity(fluidCube, 1, j, 100, 0);
					fluidCubeAddDensity(fluidCube, 1, 25, 1000);
					fluidCubeAddDensity(fluidCube, 1, 50, 1000);
					fluidCubeAddDensity(fluidCube, 1, 75, 1000);
				}
				currentFluid -= 25;
			}
			else {
					FluidCubeAddVelocity(fluidCube, 1, 25, 100, 0);
					FluidCubeAddVelocity(fluidCube, 1, 50, 100, 0);
					FluidCubeAddVelocity(fluidCube, 1, 75, 100, 0);
					fluidCubeAddDensity(fluidCube, 0, 0, 0);
					currentFluid += 0.1;
				
			}
			//Add density from bottom to top of the very left hand side of the screen w/ a constant velocity.
			/*for (int i = 1; i <= 10; i++) {
				for (int j = 1; j <= GRID_SIZE; j++) {
					
				}
			}*/
			
		}

		//Exit Simulation
		if (sdlEvent.type == SDL_KEYDOWN) {
			if (sdlEvent.key.keysym.sym == SDLK_ESCAPE) {
				running = false;
			}
		}

		//Get density & velocity at mouse posi and cout + get mouse posi in world cords
		SDL_GetMouseState(&mousePosiX, &mousePosiY); 
		int gridX = (float)(((GRID_SIZE)/ (float)DISPLAY_SIZE_X) * mousePosiX);
		int gridY = (float)(((GRID_SIZE) / (float)DISPLAY_SIZE_Y) * mousePosiY);

		//Get mouse cell position on grid (pixel cords)
		int adjustedMousePosiY = (DISPLAY_SIZE_Y - mousePosiY);
		int mouseGridPosiX = round((float)mousePosiX / CELL_SIZE_X);
		int mouseGridPosiY = round((float)adjustedMousePosiY / CELL_SIZE_Y);

		//If mouse button is pressed
		if (sdlEvent.type == SDL_MOUSEBUTTONDOWN || SDL_KEYDOWN) {

			//-----------------KEYBOARD CONTROLS

			//If Left Ctrl then add density to area
			if (sdlEvent.button.button == SDL_SCANCODE_LCTRL) {
				fluidCubeAddDensity(fluidCube, mouseGridPosiX, mouseGridPosiY, 1000);
			}

			//If Left Shift then deduct density from area (seen as black fluid on display
			if (sdlEvent.button.button == SDL_SCANCODE_LSHIFT) {
				fluidCubeAddDensity(fluidCube, mouseGridPosiX, mouseGridPosiY, -1000);
			}

			//-----------------MOUSE CONTROLS

			//If the left mouse button is pressed then velocity direction goes left
			if (sdlEvent.button.button == SDL_BUTTON_LEFT) {
				FluidCubeAddVelocity(fluidCube, mouseGridPosiX, mouseGridPosiY, -1000.0f, 0);
			}
			
			//If the right mouse button is pressed then velocity direction goes right
			if(sdlEvent.button.button == SDL_BUTTON_RIGHT){
				FluidCubeAddVelocity(fluidCube, mouseGridPosiX, mouseGridPosiY, 1000.0f, 0);

			}

		}

		//TODO: Some timestep stuff
		FluidCubeTimeStep(fluidCube);

		
		//TODO: Proper drawing code
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glClearColor(0.0f, 0.0f, 1.0f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT);
		glMatrixMode(GL_MODELVIEW); //Model Matrix
		glLoadIdentity(); //Identity Matrix
		glOrtho(0, DISPLAY_SIZE_X, 0, DISPLAY_SIZE_Y, 0.0f, 1.0f); //Projection Matrix

		glPointSize(10.0f);
		glLineWidth(5.0f);

		drawFluidDensity(fluidCube);
		drawFluidVelocity(fluidCube);

		SDL_GL_SwapWindow(gWindow); // swap buffers
	}
	return 0;
}