#include <iostream>
#include "bctype.h"
class bcdirichlet :public bctype
{
public:
bcdirichlet(arc *u,int side,char * vectortype,arc *x,arc *y);

void impose(arc *u);
~bcdirichlet();
};
