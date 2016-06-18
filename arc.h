#include "mpi.h"
#include <iostream.h>
using namespace std;

class arc
{
private:
// Number of grid points in my processor
  int m,  n,mloc,nloc;

// Number of grid points globally
   int mg, ng;

// Number of  processors, and my processor ID.   
   int nproc, myproc;

// Mapping from global to local coordinates iglobal = mmap + ilocal   
   int mmap, nmap;

// IDs of neighbour processors, in the two dimensional processor enumeration.
   int mp, mm, np, nm;

// Overlap size (should be even, see Lecture Notes)
   int olap;

// The array
   double *v;
   double *vp;
   double *tempv;
// Setup function, called by the constructors
   void set_up( int , int , int );

// Hide standard constructor
arc() {};

// Topology
int  message;
MPI_Comm comm_cart;
int My_ID;
int np1,np2;
int p1,p2;
int dims[2];
int mycoords[2];
int periods[2];
//output 
char* file;

public:
// Constructor, giving global number of points and overlap.
arc( int i, int j, int olp );

// Copy constructor
arc( const arc& a ); 

// Construct by reading the array from a file
arc( char* filename, int ol );

// Destructor, ( gives back memory )
~arc();

// Write the array to a file.
void write_file( char* filename );
// Update points in overlap.
double * comm_boundary(double * v);
// Copy assignment
arc operator=(const arc a); 
arc operator+( arc  a);
inline double &operator()( int i, int j ){return v[i + m*j];}
arc  f( arc a );
void g(arc *a);
void g(arc * a,int side);
void  zeros(int side);
void  zeros();
void fun(arc *a,arc *b);
double * getv();
int * getlocal_grid();
arc  u(arc &a,arc &b );
arc D0x(arc &x, arc &y);
void out();
void set_neighbours_IDs();
void set_neighbours_coordinates();
void get_index_transform();
void read__a_part(char*  file);
void write__a_part(char* file);

void impose_neumann(int side);
int ind( int i , int j );
};









