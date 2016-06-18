#include <iostream>
#include "bcneumann.h"

bcneumann::bcneumann(arc *u,int side,char * vectortype,arc *x,arc *y): bctype(u,side,vectortype,x,y) {}
void bcneumann::impose(arc *u)
{u->impose_neumann(side);}



