#include <iostream>

#ifndef BCTYPE_H
#define BCTYPE_H
#include "bctype.h"

class bcneumann :public bctype
{
public:
bcneumann(arc *u,int side,char * vectortype,arc *x,arc *y);
void impose(arc *u);
~bcneumann();
};
#endif