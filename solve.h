#include <iostream>
#include "nssolve.h"

class solve :public   nssolve
{
public:
solve(char * grid_x,char* grid_y,double eps,double dt);
void initial_cond(); // prepare u, v, p at t = 0
void init_boundary(); // prepare bc
void Impose_boundary();
~solve();
};
