#include "arc.h"
#include "fstream" 
#include <math.h>
#include <iostream>

arc::arc (int i, int j, int olp )
{
 set_up( i, j, olp );
}

void arc::set_up( int i, int j, int olp )
{
mg=i;
ng=j;
olap=olp;
// Determine the number of processors and my identity .
MPI_Comm_size(MPI_COMM_WORLD,&nproc);
MPI_Comm_rank(MPI_COMM_WORLD, &myproc);
// // Initialize dimensions.
for (int z=0;z<2;z++) {dims[z]=0;periods[z]=0;}
// // Decompose the number of processors.
MPI_Dims_create(nproc,2,dims);
np1=dims[0];
np2=dims[1];
// Arange a cartesian topology.
MPI_Cart_create(MPI_COMM_WORLD,2,dims,periods,0,&comm_cart);
// //  Determine my coordinates in the two-dimensional processor array.
MPI_Cart_coords (comm_cart,myproc,2,mycoords);
p1=mycoords[0];
p2=mycoords[1];
// // Determine my identity in the two-dimensional processor array.
MPI_Cart_rank(comm_cart,mycoords,&My_ID);
// // Determine the number of grid points in my processor (m,n)
// // Get the local to global index transform,
// // iglobal = mmap + ilocal, jglobal = nmap + jlocal   (mmap,nmap)
get_index_transform();
// // Allocate the array of local size, as a one dimensional array
v = new double[m*n];
// // Determine and store the identity of the neighbouring
set_neighbours_IDs();
}

void  arc::set_neighbours_IDs()
{
MPI_Cart_shift(comm_cart, 1, 1, &nm, &np);
MPI_Cart_shift(comm_cart, 0, 1, &mm, &mp);
}

void arc::get_index_transform()
{
int mremain,nremain,s1,s2;
s1=0; s2=0;
s1=  (mg+(np1-1)*olap)/(np1);     m=s1;
s2=  (ng+(np2-1)*olap)/(np2);     n=s2;

mremain=(mg+(np1-1)*olap)  % np1;
nremain=(ng+(np2-1)*olap)  % np2;

mmap=(p1)*(s1-olap)+min(p1,mremain);
nmap=(p2)*(s2-olap)+min(p2,nremain);

if (p1+1<=mremain) {m= m+1;}
if (p2+1<=nremain) {n= n+1;}

}

arc arc::operator=(const arc  a)
{
delete [] v;
mg = a.mg;
ng = a.ng;
olap=a.olap;
set_up( mg, ng, olap );

for (int i = 0; i < m*n; i++)  v[i] = a.v[i];
return *this;
}
arc::arc( const arc& a )
{
mg = a.mg;
ng = a.ng;
olap=a.olap;
set_up( mg, ng, olap );
for (int i = 0; i < m*n; i++) 
  v[i] = a.v[i];
}

arc arc::operator+( arc  a)
{
arc res(a.mg,a.ng,a.olap);
for (int i = 0; i < m*n; i++)
  res.v[i] = a.v[i]+v[i];
return res;
}

arc arc::f( arc  a)
{
arc res(a.mg,a.ng,a.olap);
for (int i = 0; i < m*n; i++)
  res.v[i] = sin(a.v[i]/10)*sin(a.v[i]/10)*cos(a.v[i]/10);
return res;
}

int * arc::getlocal_grid()
{
  int * grid;
  grid=new int[1];
  grid[0]=m;
  grid[1]=n;

return grid;
}

double * arc::getv()
{
return v;
}
void arc::g(arc * a,int side)
{
int i,j;

  if (side==0)
   {
  //left
  if (mm==-1)  // left boundary
   {
   i=1;
   for(j=1; j<=n ;j++)
   {
   v[i-1+m*(j-1)]=.1*a->v[i-1+m*(j-1)]*(3-a->v[i-1+m*(j-1)]);
   }
   }
   }

    if (side==1)
   {
    if (nm==-1)  // top boundary
   {
  //top
    j=1;
    for(i=1; i<=i ;i++)
   {
   v[i-1+m*(j-1)]=.1*a->v[i-1+m*(j-1)]*(3-a->v[i-1+m*(j-1)]);
   }
   }
   }
     if (side==2)
   {
  //right
   if (mp==-1)  // right boundary
   {
   i=m;
   for(j=1; j<=n ;j++)
   {
   v[i-1+m*(j-1)]=.1*a->v[i-1+m*(j-1)]*(3-a->v[i-1+m*(j-1)]);
   }
   }
   }
     if (side==3)
   {
  //down
   if (np==-1)  // down boundary
   {
  //top
    j=n;
    for(i=1; i<=i ;i++)
   {
   v[i-1+m*(j-1)]=.1*a->v[i-1+m*(j-1)]*(3-a->v[i-1+m*(j-1)]);
   }
   }
   }

}

void arc::g(arc * a)
{
int i;
for(i=0; i<=m*n ;i++){v[i]=.1*a->v[i]*(3-a->v[i]);}
}

void arc::impose_neumann(int side)
{
  int i,j;

  if (side==0)
   {
  //left
  if (mm==-1)  // left boundary
   {
   i=1;
   for(j=1; j<=n ;j++)
   {
   v[i-1+m*(j-1)]=v[i+m*(j-1)];
   }
   }
   }

    if (side==1)
   {
    if (nm==-1)  // top boundary
   {
  //top
    j=1;
    for(i=1; i<=m ;i++)
   {
   v[i-1+m*(j-1)]=v[i-1+m*(j)];
   }
   }
   }
     if (side==2)
   {
  //right
   if (mp==-1)  // right boundary
   {
   i=m;
   for(j=1; j<=n ;j++)
   {
   v[i-1+m*(j-1)]=v[i-1+m*(j-1)];
   }
   }
   }
     if (side==3)
   {
  //down
   if (np==-1)  // down boundary
   {
  //top
    j=n;
    for(i=1; i<=m ;i++)
   {
   v[i-1+m*(j-1)]=v[i-1+m*(j-2)];
   }
   }
   }
}

void arc::zeros()
{
  int i,j;
  for (int i = 0; i < m*n; i++)   v[i] = 0;
}
void arc::fun(arc *a,arc *b)
{
  int i,j;
  for (int i = 0; i < m*n; i++)   v[i] = .5*(a->v[i]+b->v[i]);
}

void arc::zeros(int side)
{

  int i,j;

  if (side==0)
   {//left
  if (mm==-1)  // left boundary
   {i=1;for(j=1; j<=n ;j++){v[i-1+m*(j-1)]=0;}}
   }

if (side==1)
   {//top

    if (nm==-1)  // top boundary
   {j=1;for(i=1; i<=m ;i++)     {v[i-1+m*(j-1)]=0; }}
   }

if (side==2)
   {//right
   if (mp==-1)  // right boundary
   {i=m;for(j=1; j<=n ;j++)     {v[i-1+m*(j-1)]=0;}}
   }

if (side==3)
   {//down
   if (np==-1)  // down boundary
   {j=n;for(i=1; i<=m ;i++)   {v[i-1+m*(j-1)]=0;}}
   }

}

arc::arc( char* file, int olp )
{
MPI_Comm_size(MPI_COMM_WORLD, &nproc);
MPI_Comm_rank(MPI_COMM_WORLD, &myproc);
if( myproc == 0 )
  {
    //     Open file and read dimensions, mg, ng 
    int nr;
    FILE *fp;
    fp = fopen (file , "r" );
    nr  = fread(& mg, sizeof(int ), 1, fp );
    nr  = fread(& ng, sizeof(int ), 1, fp );

    fclose(fp);

    //     Send (mg,ng) to all processors
    for(int des=1;des<nproc;des++)
	{      
	  MPI_Send(&mg,1,MPI_INT,des,99,MPI_COMM_WORLD);
	  MPI_Send(&ng,1,MPI_INT,des,99,MPI_COMM_WORLD);
	}
  }

else

  {
    //     receive (mg,ng) from processor zero.
    MPI_Status status;
    MPI_Recv(&mg,1,MPI_INT,0,99,MPI_COMM_WORLD,&status);
    MPI_Recv(&ng,1,MPI_INT,0,99,MPI_COMM_WORLD,&status);
  } 

set_up( mg, ng, olp);

if( myproc == 0 )
  {
    //  read_my_part
    int position,nr,i,j,ind;
    FILE *fp ;
    fp=fopen(file,"r+");
    for( int j=1 ; j<=n; j++ )
	{
	  position=(mmap)+(j-1+nmap)*mg+1;
	  fseek( fp, position*sizeof(double), SEEK_SET );
	  nr  = fread( v + m*(j-1), sizeof(double), m, fp );
	}
    fclose(fp);
  }
else
  {
    //  receive a message from myproc-1
    MPI_Status status;
    MPI_Recv(&message,1,MPI_INT,myproc-1,99,MPI_COMM_WORLD,&status);
    //   read_my_part
    read__a_part(file);
  }
if( myproc < nproc-1 )
  {
    message=2;
    //  send message to myproc + 1
    MPI_Send(&message,1,MPI_INT,myproc+1,99,MPI_COMM_WORLD);
  }
}

void arc::read__a_part(char* file)
{
//Open file
int position,nr,i,j,ind;
FILE *fp ;
fp=fopen(file,"r+");
for( int j=1 ; j<=n; j++ )
  {
    position=(mmap)+(j-1+nmap)*mg+1;
    fseek( fp, position*sizeof(double), SEEK_SET );
    nr  = fread( v + m*(j-1), sizeof(double), m, fp );
  }
fclose(fp);
}

void arc::write_file( char* file )
{

if( myproc == 0 )
  {
    //     Open file and write dimensions, mg, ng 
    int nr;
    FILE *fp;
    fp = fopen (file , "w" );
    nr  = fwrite(&mg, sizeof(int ), 1, fp );
    nr  = fwrite(&ng, sizeof(int ), 1, fp );
    fclose(fp);
  }

 if( myproc == 0 )
  { 
    //  write_my_part
    int position,nr,i,j,ind;
    FILE *fp ;
     fp=fopen(file,"r+");
      for( int j=1 ; j<=n; j++ )
	{
	  position=mmap+(j-1+nmap)*mg+1;
	  fseek( fp, position*sizeof(double), SEEK_SET );
	  nr  = fwrite( v + m*(j-1), sizeof(double), m, fp );
	}
    fclose(fp);
  }
else
  {
    //  receive a message from myproc-1
    MPI_Status status;
    MPI_Recv(&message,1,MPI_INT,myproc-1,99,MPI_COMM_WORLD,&status);
    //   write_my_part
     write__a_part(file);
  }
if( myproc < nproc-1 )
  {
    message=1;
    //  send message to myproc + 1
    MPI_Send(&message,1,MPI_INT,myproc+1,99,MPI_COMM_WORLD);
  }
}

void arc::write__a_part(char* file)
{
// Open file
int position,nr,i,j,ind;
FILE *fp ;
 fp=fopen(file,"r+");
 for( int j=1 ; j<=n; j++ )
  {
    position=(mmap)+(j-1+nmap)*mg+1;
    fseek( fp, position*sizeof(double), SEEK_SET );
    nr  = fwrite( v + m*(j-1), sizeof(double), m, fp );
  } 
fclose(fp);
}

arc arc::D0x(arc &x, arc &y)
{
int i,j;
double x_xi,x_etha,y_xi,y_etha,u_xi,u_etha,J;
double *temp=new double[m*n+2];
double *xx=new double[m*n];
double *yy=new double[m*n];
// Allocate new arc object
arc vprime(mg,ng,olap);
// define grid points

for( j= 1 ; j<=n ; j++ ) {
 for( i= 1 ;i<= m ; i++ ){
  xx[i-1 + m*(j-1)]=x.v[i-1 + m*(j-1)];
  yy[i-1 + m*(j-1)]=y.v[i-1 + m*(j-1)];
}}

 //  calculate D0x in my part
 for(j=1 ; j<=n-2 ;j++){
 for(i=1 ;i<= m-2; i++ )
 {
     x_xi=(xx[i+1 + m*(j)]-xx[i-1 + m*(j)])/2;
     y_xi=(yy[i+1 + m*(j)]-yy[i-1 + m*(j)])/2;
     u_xi=(v[i+1+ m*(j)]-v[i-1 + m*(j)])/2;
     x_etha=(xx[i + m*(j+1)]-xx[i + m*(j-1)])/2;
     y_etha=(yy[i + m*(j+1)]-yy[i + m*(j-1)])/2;
     u_etha=(v[i + m*(j+1)]-v[i + m*(j-1)])/2;
     J=(x_xi*y_etha-x_etha*y_xi);
     temp[i+ m*(j)+2]=(1/J)*(u_xi*y_etha-u_etha*y_xi);
 }}

//  Boundary conditions treatment

   if (mm==-1)  // left boundary
   {
    i=0;
   for(j=1; j<=n-2 ;j++){
    x_xi=(-3*xx[i + m*(j)]+4*xx[i+1 + m*(j)]-xx[i+2 + m*(j)])/2;
    y_xi=(-3*yy[i + m*(j)]+4*yy[i+1 + m*(j)]-yy[i+2 + m*(j)])/2;
    u_xi=(-3*v[i + m*(j)]+4*v[i+1 + m*(j)]-v[i+2 + m*(j)])/2;
    x_etha=(xx[i + m*(j+1)]-xx[i + m*(j-1)])/2;
    y_etha=(yy[i + m*(j+1)]-yy[i + m*(j-1)])/2;
    u_etha=(v[i + m*(j+1)]-v[i + m*(j-1)])/2;
    J=(x_xi*y_etha-x_etha*y_xi);
    temp[i+ m*(j)+2]=(1/J)*(u_xi*y_etha-u_etha*y_xi);
   }
    if (nm==-1) 
    {
    j=0;
    x_etha=(-3*xx[i + m*(j)]+4*xx[i + m*(j+1)]-xx[i + m*(j+2)])/2;
    y_etha=(-3*yy[i + m*(j)]+4*yy[i + m*(j+1)]-yy[i + m*(j+2)])/2;
    u_etha=(-3*v[i + m*(j)]+4*v[i + m*(j+1)]-v[i + m*(j+2)])/2;
    J=(x_xi*y_etha-x_etha*y_xi);
    temp[i+ m*(j)+2]=(1/J)*(u_xi*y_etha-u_etha*y_xi);
    }
    if (np==-1)
    {
    j=n-1;
    J=(x_xi*y_etha-x_etha*y_xi);
    temp[i+ m*(j)+2]=(1/J)*(u_xi*y_etha-u_etha*y_xi);
    x_etha=(-3*xx[i + m*(j)]+4*xx[i + m*(j-1)]-xx[i + m*(j-2)])/2;
    y_etha=(-3*yy[i + m*(j)]+4*yy[i + m*(j-1)]-yy[i + m*(j-2)])/2;
    u_etha=(-3*v[i + m*(j)]+4*v[i + m*(j-1)]-v[i + m*(j-2)])/2;
    }
   }

 if (mp==-1) // Right Boundary
  {
     i=m-1;
    for(j=1; j<=n-2 ;j++){
     x_xi=(-3*xx[i + m*(j)]+4*xx[i-1 + m*(j)]-xx[i-2 + m*(j)])/2;
     y_xi=(-3*yy[i + m*(j)]+4*yy[i-1 + m*(j)]-yy[i-2 + m*(j)])/2;
     u_xi=(-3* v[i + m*(j)]+4* v[i-1 + m*(j)]- v[i-2 + m*(j)])/2;
     x_etha=(xx[i + m*(j+1)]-xx[i + m*(j-1)])/2;
     y_etha=(yy[i + m*(j+1)]-yy[i + m*(j-1)])/2;
     u_etha=(v[i + m*(j+1)]-v[i + m*(j-1)])/2;
     J=(x_xi*y_etha-x_etha*y_xi);
     temp[i+ m*(j)+2]=(1/J)*(u_xi*y_etha-u_etha*y_xi);
  }
    if (nm==-1)
    {
    j=0;
    x_etha=(-3*xx[i + m*(j)]+4*xx[i + m*(j+1)]-xx[i + m*(j+2)])/2;
    y_etha=(-3*yy[i + m*(j)]+4*yy[i + m*(j+1)]-yy[i + m*(j+2)])/2;
    u_etha=(-3*v[i + m*(j)]+4*v[i + m*(j+1)]-v[i + m*(j+2)])/2;
    J=(x_xi*y_etha-x_etha*y_xi);
    temp[i+ m*(j)+2]=(1/J)*(u_xi*y_etha-u_etha*y_xi);
    }
    if (np==-1)
    {
    j=n-1;
    J=(x_xi*y_etha-x_etha*y_xi);
    temp[i+ m*(j)+2]=(1/J)*(u_xi*y_etha-u_etha*y_xi);
    x_etha=(-3*xx[i + m*(j)]+4*xx[i + m*(j-1)]-xx[i + m*(j-2)])/2;
    y_etha=(-3*yy[i + m*(j)]+4*yy[i + m*(j-1)]-yy[i + m*(j-2)])/2;
    u_etha=(-3*v[i + m*(j)]+4*v[i + m*(j-1)]-v[i + m*(j-2)])/2;
    }
 }

  if (nm==-1)   // north boundary
  {
    j=0;
  for(i=1 ;i<= m-2; i++ )
  {
   x_xi=(xx[i+1 + m*(j)]-xx[i-1 + m*(j)])/2;
   y_xi=(yy[i+1 + m*(j)]-yy[i-1 + m*(j)])/2;
   u_xi=(v[i+1+ m*(j)]-v[i-1 + m*(j)])/2;
   x_etha=(-3*xx[i + m*(j)]+4*xx[i + m*(j+1)]-xx[i + m*(j+2)])/2;
   y_etha=(-3*yy[i + m*(j)]+4*yy[i + m*(j+1)]-yy[i + m*(j+2)])/2;
   u_etha=(-3*v[i + m*(j)]+4*v[i + m*(j+1)]-v[i + m*(j+2)])/2;
   J=(x_xi*y_etha-x_etha*y_xi);
   temp[i+ m*(j)+2]=(1/J)*(u_xi*y_etha-u_etha*y_xi);
  }
  }

 if (np==-1)  // south boundary
 {
    j=n-1;
   for(i=1 ;i<= m-2; i++ )
  {
    x_xi=(xx[i+1 + m*(j)]-xx[i-1 + m*(j)])/2;
    y_xi=(yy[i+1 + m*(j)]-yy[i-1 + m*(j)])/2;
    u_xi=(v[i+1+ m*(j)]-v[i-1 + m*(j)])/2;
    x_etha=(-3*xx[i + m*(j)]+4*xx[i + m*(j-1)]-xx[i + m*(j-2)])/2;
    y_etha=(-3*yy[i + m*(j)]+4*yy[i + m*(j-1)]-yy[i + m*(j-2)])/2;
    u_etha=(-3*v[i + m*(j)]+4*v[i + m*(j-1)]-v[i + m*(j-2)])/2;
    J=(x_xi*y_etha-x_etha*y_xi);
    temp[i+ m*(j)+2]=(1/J)*(u_xi*y_etha-u_etha*y_xi);
  }
 }

for( j= 0 ; j<=n-1 ; j++ ) {for( i= 0 ;i<= m-1 ; i++ ){
   vprime. v[i+ m*(j)]=temp[i+ m*(j)+2]; ;
}}
delete [] temp;
// exchange ghost points with neighbors
v=comm_boundary(v);
return vprime;
}

double * arc::comm_boundary(double * v)
{
// // Boundary communication for non-overlapping decomposition
MPI_Status status;
int i,j;
double *vv = new double[m*n];
double * tempv = new double[m*n+2];
// -------------------   copy internal data  --------------------------------
for( int j= 1 ; j<=n ; j++ ){for( int i= 1 ; i<=m ; i++ )
{tempv[i-1 + (m)*(j-1)+2]=v[i-1 + m*(j-1)];}}

// -------------------   exchange ghost points in     x direction  ----------
for( j= 1 ; j<=n ; j++ ) {
MPI_Sendrecv(v+1+m*(j-1), 1,MPI_DOUBLE,mm,100,
tempv+m-1+(m)*(j-1)+2,1,MPI_DOUBLE,mp,100,comm_cart,&status);}

for( j= 1 ; j<=n ; j++ ) {
MPI_Sendrecv(v+m-2+m*(j-1), 1,MPI_DOUBLE,mp,100
,tempv+0+(m)*(j-1)+2,1,MPI_DOUBLE,mm,100,comm_cart,&status);}

// -------------------   exchange ghost points in     y direction  ------------------
for( i= 1 ; i<=m ; i++ ){
MPI_Sendrecv(v+i-1+m*(1), 1,MPI_DOUBLE,nm,100,
tempv+i-1+(m)*(n-1)+2,1,MPI_DOUBLE,np,100,comm_cart,&status);}

for( i= 1 ; i<=m ; i++ )
{MPI_Sendrecv(v+i-1 + m* (n-2),1,MPI_DOUBLE,np,100,
tempv+i-1+(m)*(0)+2,1,MPI_DOUBLE,nm,100,comm_cart,&status );};
for( int j= 1 ; j<=n ; j++ ){for( int i= 1 ; i<=m ; i++ )
{vv[i-1 + (m)*(j-1)]=tempv[i-1 + m*(j-1)+2];
}}
delete [] tempv;
return vv;
}
void arc::out()
{
  int i,j;
 if(myproc ==0)
  {
    int i,j;
   cout << "process   " << myproc << endl;
   for( j=1; j<=n ;j++){for( i= 1;i<= m; i++ )
  {cout << v[i-1+(m)*(j-1)]<<" | ";} cout << "== " << nmap+j << endl;  }
  }}

arc::~arc()
{
delete[] v;
}
