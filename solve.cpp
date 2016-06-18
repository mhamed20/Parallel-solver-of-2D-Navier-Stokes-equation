#include "solve.h"

solve::solve(char * grid_x,char* grid_y,double eps,double dt) :nssolve(grid_x,grid_y,eps,dt) {}
void solve::initial_cond() // prepare u, v, p at t = 0
{
u=new arc(mg,ng,2);
v=new arc(mg,ng,2);
p=new arc(mg,ng,2);

u->g(y);
v->zeros();
p->zeros();
MPI_Barrier(MPI_COMM_WORLD);
}

void solve::init_boundary() // prepare bc
{
for(int s=0;s<=3;s++) bc[s]= new nsbc;;
// Left boundary
bc[0]->u=new bcdirichlet(u,0,"u",x,y);
bc[0]->v=new bcdirichlet(v,0,"v",x,y);
bc[0]->p=new bcneumann(p,0,"p",x,y);

 // top boundary
bc[1]->u=new bcdirichlet(u,1,"u",x,y);
bc[1]->v=new bcdirichlet(v,1,"v",x,y);
bc[1]->p=new bcneumann(p,1,"p",x,y);

// right boundary
bc[2]->u=new bcneumann(u,2,"u",x,y);
bc[2]->v=new bcdirichlet(v,2,"v",x,y);
bc[2]->p=new bcdirichlet(p,2,"p",x,y);

// down boundary
bc[3]->u=new bcdirichlet(u,3,"u",x,y);
bc[3]->v=new bcdirichlet(v,3,"v",x,y);
bc[3]->p=new bcneumann(p,3,"p",x,y);
Impose_boundary();
}
void solve::Impose_boundary() // impose bc
{
bc[0]->u->impose(u);bc[0]->v->impose(v);bc[0]->p->impose(p);
bc[1]->u->impose(u);bc[1]->v->impose(v);bc[1]->p->impose(p);
bc[2]->u->impose(u);bc[2]->v->impose(v);bc[2]->p->impose(p);
bc[3]->u->impose(u);bc[3]->v->impose(v);bc[3]->p->impose(p);
}
solve::~solve()  {}

