#include "nssolve.h"
#include<cmath>

extern "C" int resid_(int*, int*, double*, double*, double*, double*,double*, double*, double*, double*);

nssolve::nssolve(char * grid_x,char* grid_y,double eps,double dt) // read x, y coordinates
{
    //     Open file and read dimensions, mg, ng
    int nr;
    FILE *fp;
    fp = fopen (grid_x , "r" );
    nr  = fread(& mg, sizeof(int ), 1, fp );
    nr  = fread(& ng, sizeof(int ), 1, fp );
    fclose(fp);
    // Define grid points
    x=new arc(grid_x,2);
    y=new arc(grid_y,2);
}
nssolve::nssolve() {} // Default constructor
void nssolve::timesteps(double dt, int nsteps) // time stepping
{
int * grid;
grid=x->getlocal_grid();
int ni=grid[0];
int nj=grid[1];
double res[3*ni*nj];
for(int st=1;st<=nsteps;st++)
{
arc *un=u;
arc *vn=v;
arc *pn=p;
c_resid(ni, nj,u->getv(),v->getv(),p->getv(), res,x->getv(),y->getv(),eps,dt);

arc *u1=u;
arc *v1=v;
arc *p1=p;

c_resid(ni, nj,u->getv(),v->getv(),p->getv(), res,x->getv(),y->getv(),eps,dt);


arc *u2=u;
arc *v2=v;
arc *p2=p;

u->fun(un,u2);
v->fun(vn,v2);
p->fun(pn,p2);

}
// error monitoring:
double U_err=res[0];
double V_err=res[1];
double P_err=res[2];
MPI_Allreduce(&U_err, &error_global_u, 1,MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
MPI_Allreduce(&V_err, &error_global_v, 1,MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
MPI_Allreduce(&P_err, &error_global_p, 1,MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

}

void nssolve::c_resid(int ni, int nj, double *u, double *v, double *p,
double *res, double *x, double *y, double eps, double dt)
{
int ni_tmp = ni;
int nj_tmp = nj;
double eps_tmp = eps;
double dt_tmp = dt;
resid_(&ni_tmp, &nj_tmp, u, v, p, res, x, y, &eps_tmp, &dt_tmp);
}

void nssolve::write_solution(char * filename,double t1)
{
int myproc,nproc;
MPI_Comm_rank(MPI_COMM_WORLD, &myproc);
double rtmin,rtmax;
double t2=MPI_Wtime()-t1;
MPI_Reduce(&t2,&rtmax,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
MPI_Reduce(&t2,&rtmin,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);


if(myproc ==0)  {
 cout << "Real time (s) : " << rtmin<< "  " << rtmax << endl;
 cout << "Global error of u : " <<error_global_u << endl;
 cout << "Global error of v : " <<error_global_v << endl;
 cout << "Global error of p : " <<error_global_p << endl;
 }

// writing to files
u->write_file("u.dat");
v->write_file("v.dat");
p->write_file("p.dat");
} // obvious

nssolve::~nssolve()  {}
