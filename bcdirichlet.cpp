#include <iostream>
#include "bcdirichlet.h"

bcdirichlet::bcdirichlet(arc *u,int side,char * vectortype,arc *x,arc *y): bctype(u,side,vectortype,x,y) {}
void bcdirichlet::impose(arc *u)
{

  if (side>0) u->zeros(side);
  else
{if (vectortype=="u")    {u->g(y,side);}}
}

bcdirichlet::~bcdirichlet()  {}