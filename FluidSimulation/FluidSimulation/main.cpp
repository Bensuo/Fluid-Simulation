#include <glm/glm.hpp>
#include <SDL.h>
#include <SDL_ttf.h>
#include <GL\glew.h>
#undef main
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
using namespace std;

#define IX(x,y) ((x)+(N+2)*(y))

#define CELL_SIZE 80
#define GRID_SIZE 10


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
typedef struct FluidCube FluidCube; // Ask Marco

int cellOffSetX = 8;
int cellOffSetY = 2;

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

	fluidCube->cellStates = new CellState[(N + 2) * (N + 2)]();
	fluidCube->particles = new Particle[N*N*2]();

	return fluidCube;
}


// Add Density (Dye)
void fluidCubeAddDensity(FluidCube *fluidCube, int x, int y, float amount) {
	int N = fluidCube->size;
	fluidCube->s[IX(x, y)] += amount;
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


	int i;
	for (i = 1; i <= N; i++) {
		x[IX(0, i)] = b == 1 ? -x[IX(1, i)] : x[IX(1, i)];
		x[IX(N + 1, i)] = b == 1 ? -x[IX(N, i)] : x[IX(N, i)];
		x[IX(i, 0)] = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
		x[IX(i, N + 1)] = b == 2 ? -x[IX(i, N)] : x[IX(i, N)];
	}
	x[IX(0, 0)] = 0.5*(x[IX(1, 0)] + x[IX(0, 1)]);
	x[IX(0, N + 1)] = 0.5*(x[IX(1, N + 1)] + x[IX(0, N)]);
	x[IX(N + 1, 0)] = 0.5*(x[IX(N, 0)] + x[IX(N + 1, 1)]);
	x[IX(N + 1, N + 1)] = 0.5*(x[IX(N, N + 1)] + x[IX(N + 1, N)]);
}

//Implementation of Gauss-Seidel relaxation formula.
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

			if (y < 0.5) {
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
void add_source(int N, float * x, float * s, float dt) {
	int i, size = (N + 2) * (N + 2);
	for (i = 0; i < size; i++)
	{
		x[i] += dt*s[i];
	}

}

void getForces(int N, float * s, float * Vx0, float * Vy0)
{
	int i, size = (N + 2) * (N + 2);
	for (i = 0; i < size; i++)
	{
		s[i] = 0;
		Vx0[i] = 0;
		Vy0[i] = 10.0f;
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
	float *s = fluidCube->s;
	float *density = fluidCube->density;
	Particle* particles = fluidCube->particles;
	CellState *cellStates = fluidCube->cellStates;


	//Velocity step
	add_source(N, Vx, Vx0, dt);
	add_source(N, Vy, Vy0, dt);
	diffuse(N, 1, Vx0, Vx, visc, dt);
	diffuse(N, 2, Vy0, Vy, visc, dt);

	project(N, Vx0, Vy0, Vx, Vy);

	advect(N, 1, Vx, Vx0, Vx0, Vy0, dt);
	advect(N, 2, Vy, Vy0, Vx0, Vy0, dt);

	project(N, Vx, Vy, Vx0, Vy0);

	//Density step
	add_source(N, density, s, dt);
	diffuse(N, 0, s, density, diff, dt);
	advect(N, 0, density, s, Vx, Vy, dt);


	getForces(N, s, Vx0, Vy0);
}



std::string ValuesToText(FluidCube* fluidCube)
{
	std::string value = "";
	stringstream ss;
	ss << fixed << setprecision(2);
	int N = fluidCube->size + 2;
	int size = N*N;

	for (int i = 0; i < size; i++)
	{
		ss << fluidCube->density[i] << " ";
		if ((i + 1) % (N) == 0 && i != 0)
		{
			ss << '\n';
		}
	}
	value = ss.str();
	return value;
}

const int SCREEN_WIDTH = 1900;
const int SCREEN_HEIGHT = 1060;
SDL_Window* gWindow = NULL;
SDL_Renderer* gRenderer = NULL;
TTF_Font* gFont = NULL;

void drawFluidVelocity(FluidCube* fluidCube) {
	int i, j;
	int N = fluidCube->size;
	float x, y;
	float h = 1.0f;
	float *Vx = fluidCube->Vx;
	float *Vy = fluidCube->Vy;
	int offSetX = cellOffSetX * CELL_SIZE;
	int offSetY = cellOffSetY * CELL_SIZE;
	float cellSize = CELL_SIZE;

	glColor3f(1.0f, 1.0f, 1.0f);
	glLineWidth(5.0f);

	glBegin(GL_LINES);

	for (i = 1; i <= N; i++) {
		x = (i - 1.0f) * h;
		for (j = 1; j <= N; j++) {
			y = (j - 1.0f) * h;

			glVertex2f(x*cellSize + offSetX, y*cellSize + offSetY);
			glVertex2f((x + Vx[IX(i, j)]) * cellSize + offSetX, (y + Vy[IX(i, j)]) * cellSize + offSetY);
		}
	}
	glEnd();
}

void drawFluidDensity(FluidCube* fluidCube) {

	int N = fluidCube->size;
	int size = N*N;
	float h = 1.0f;
	int offSetX = cellOffSetX * CELL_SIZE;
	int offSetY = cellOffSetY * CELL_SIZE;
	
	int cellSize = CELL_SIZE; //size of cells in grid
	for (int i = 0; i < N; i++)
	{
		float x = (i - 0.5f)*h;
		for (int j = 0; j < N; j++) {
			float y = (j - 0.5f)*h;
			//float density = fluidCube->density[i+j];
			float sourceAlpha = 0.1;
			float d00, d01, d10, d11;

			//calculate position of the fluid in the grid
			d00 = fluidCube->density[IX(i, j)]; 
			d01 = fluidCube->density[IX(i, j +1 )];
			d11 = fluidCube->density[IX(i + 1, j + 1)];
			d10 = fluidCube->density[IX(i + 1, j)];

			// draw density as a cube of quads

			glBegin(GL_QUADS);
			glColor4f(d00 + 0.4, d00, d00, sourceAlpha);
			glVertex3f(x*cellSize + offSetX, y*cellSize + offSetY, 0);

			glColor4f(d01 + 0.4, d01, d01, sourceAlpha);
			glVertex3f(x*cellSize + offSetX, (y + h)*cellSize + offSetY, 0);

			glColor4f(d11 + 0.4, d11, d11, sourceAlpha);
			glVertex3f((x + h)*cellSize + offSetX, (y + h)*cellSize + offSetY, 0);

			glColor4f(d10 + 0.4, d10, d10, sourceAlpha);
			glVertex3f((x + h)*cellSize + offSetX, y*cellSize + offSetY, 0);

			glEnd();

		}
	}


}

int main()
{
	SDL_Init(SDL_INIT_VIDEO);
	TTF_Init();

	// Request an OpenGL 3.0 context.
	SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
	SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 0);
	SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);

	SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);  // double buffering on
	SDL_GL_SetAttribute(SDL_GL_MULTISAMPLEBUFFERS, 1);
	SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES, 4); // Turn on x4 multisampling anti-aliasing (MSAA)

	gFont = TTF_OpenFont("Verdana.ttf", 12);
	SDL_Color White = { 0, 0, 0 };
	gWindow = SDL_CreateWindow("SDL Tutorial", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, SCREEN_WIDTH, SCREEN_HEIGHT, SDL_WINDOW_OPENGL | SDL_WINDOW_SHOWN);
	SDL_GLContext context = SDL_GL_CreateContext(gWindow); //Create OpenGL Context and binds the gWindow to SDL.
	gRenderer = SDL_CreateRenderer(gWindow, -1, SDL_RENDERER_ACCELERATED);
	SDL_SetRenderDrawColor(gRenderer, 0xFF, 0xFF, 0xFF, 0xFF);
	glEnable(GL_BLEND);
	FluidCube * fluidCube = FluidCubeCreate(GRID_SIZE, 0.001f, 0.1f, 0.1f);

	//FluidCubeAddVelocity(fluidCube, 1, 1, 100.0f, 0);
	//fluidCubeAddDensity(fluidCube, 1, 1, 10);
	//fluidCubeAddDensity(fluidCube, 4, 4, 10000);
	//fluidCubeAddDensity(fluidCube, 6, 6, 10000);
	//fluidCubeAddDensity(fluidCube, 10, 10, 10000);
	bool running = true; // set running to true

	SDL_Event sdlEvent;  // variable to detect SDL events
	while (running) {	// the event loop
		while (SDL_PollEvent(&sdlEvent)) {
			if (sdlEvent.type == SDL_QUIT)
				running = false;
		}
		const Uint8 *keys = SDL_GetKeyboardState(NULL);

		//Add Velocity
		if (keys[SDL_SCANCODE_Z])
		{
			FluidCubeAddVelocity(fluidCube, 1, 1, 1000.0f, 0);
		}

		if (keys[SDL_SCANCODE_X])
		{
			FluidCubeAddVelocity(fluidCube, 1, 1, 100.0f, 0);
		}

		if (keys[SDL_SCANCODE_C])
		{
			FluidCubeAddVelocity(fluidCube, 1, 1, 10.0f, 0);
		}


		if (keys[SDL_SCANCODE_V])
		{
			FluidCubeAddVelocity(fluidCube, 1, 1, 1.0f, 0);
		}

		//Add Density (Fluid)
		if (keys[SDL_SCANCODE_A])
		{
			fluidCubeAddDensity(fluidCube, 1, 1, 1000);
		}

		if (keys[SDL_SCANCODE_S])
		{
			fluidCubeAddDensity(fluidCube, 1, 1, 100);
		}

		if (keys[SDL_SCANCODE_D])
		{
			fluidCubeAddDensity(fluidCube, 1, 1, 10);
		}

		if (keys[SDL_SCANCODE_F])
		{
			fluidCubeAddDensity(fluidCube, 1, 1, 1);
		}

		//Add Density and Fluid

		if (keys[SDL_SCANCODE_Q])
		{
			FluidCubeAddVelocity(fluidCube, 1, 1, 1000.0f, 0);
			fluidCubeAddDensity(fluidCube, 1, 1, 1000);
		}

		if (keys[SDL_SCANCODE_W])
		{
			FluidCubeAddVelocity(fluidCube, 1, 1, 100.0f, 0);
			fluidCubeAddDensity(fluidCube, 1, 1, 100);
		}

		if (keys[SDL_SCANCODE_E])
		{
			FluidCubeAddVelocity(fluidCube, 1, 1, 10.0f, 0);
			fluidCubeAddDensity(fluidCube, 1, 1, 10);
		}

		if (keys[SDL_SCANCODE_R])
		{
			FluidCubeAddVelocity(fluidCube, 1, 1, 1.0f, 0);
			fluidCubeAddDensity(fluidCube, 1, 1, 1);
		}

		//Move Grid (Camera)

		if (keys[SDL_SCANCODE_I])
		{
			cellOffSetY -= 1;
		}

		if (keys[SDL_SCANCODE_K])
		{
			cellOffSetY += 1;
		}

		if (keys[SDL_SCANCODE_J])
		{
			cellOffSetX += 1;
		}

		if (keys[SDL_SCANCODE_L])
		{
			cellOffSetX -= 1;
		}

		//Exit Simulation
		if (sdlEvent.type == SDL_KEYDOWN) {
			if (sdlEvent.key.keysym.sym == SDLK_ESCAPE) {
				running = false;
			}
		}

		//Mouse position
		int x;
		int y; 

		SDL_GetMouseState(&x, &y);
		int gridX = x / 10;
		int gridY = y / 10;
		
		cout << fluidCube->density[gridX + gridY] << endl;

		//TODO: Some timestep stuff
		FluidCubeTimeStep(fluidCube);

		//TODO: Proper drawing code
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glOrtho(0, 1900, 0, 1060, 0.0f, 1.0f); // Reference system of the simulation
		glPointSize(10.0f);
		glLineWidth(5.0f);

		drawFluidDensity(fluidCube);
		drawFluidVelocity(fluidCube);

		SDL_GL_SwapWindow(gWindow); // swap buffers

		
		//Ben's Code to draw numbers
		/*SDL_RenderClear(gRenderer);
		string output = ValuesToText(fluidCube);
		SDL_Surface* textSurface = TTF_RenderText_Blended_Wrapped(gFont, output.c_str(), White, SCREEN_WIDTH);
		SDL_Texture* text = SDL_CreateTextureFromSurface(gRenderer, textSurface);
		int text_width = textSurface->w;
		int text_height = textSurface->h;
		SDL_FreeSurface(textSurface);
		SDL_Rect textQuad = { 20, 30, text_width, text_height };
		SDL_RenderCopy(gRenderer, text, NULL, &textQuad);
		SDL_RenderPresent(gRenderer);
		SDL_DestroyTexture(text);*/
	}
	return 0;
}