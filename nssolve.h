#include <iostream>
#include "bcneumann.h"
#include "bcdirichlet.h"
typedef struct {
bctype *u;
bctype *v;
bctype *p;
} nsbc;

class nssolve {
protected:
int mg;
int ng;
char * filename;

double eps; // epsilon

arc *u, *v, *p, *x, *y;
double error_global_u;
double error_global_v;
double error_global_p; 
nsbc *bc[4]; // boundary conditions

public:
nssolve(char * grid_x,char* grid_y,double eps,double dt); // read x, y coordinates
nssolve();
void c_resid(int ni, int nj, double *u, double *v, double *p,double *res, double *x, double *y, double eps, double dt);

void timesteps(double dt, int nsteps); // time stepping
void write_solution(char * filename, double t1); // obvious
virtual void initial_cond() = 0; // prepare u, v, p at t = 0
virtual void init_boundary() = 0; // prepare bc
~nssolve();
};
