// numpy1.cpp : Defines the exported functions for the DLL application.
//

/*

// The following ins a traversal of all element in a 2D array, acessing them in the order
// i which they are stored in memoer (i.e. acessing each element with a stride of 1 which
// is along the x-axis.
//
// a numpy array of size(y,x) has the xdim in shape[1], ydim in shape[0] (by default.. 
//		you can swap these which is called putting the indices in 'Fortran' order)
//
//		This looks strange at first, but it makes sense of it like these are in C order.. or the indices
//		are written in the common MAtrix convention that (i,j) is (row, column)

// data(y,x)
int stride1 = shape[1];
for(int y=0;y<shape[0];y++) {
	for(int x=0;x<shape[1];x++) {
		data[y*stride1 + x] += (y+1);
	}
}

*/
// TODO: this configuration needs to go out to Cmake file.
#define RESTRICT __restrict__

#if defined(_WIN32) || defined(_WIN64)
	// No heavy WinAPI stuff that needs precompiled headers..
	//#include "stdafx.h"
	#define RESTRICT __restrict
#endif

	//#define DO_DEBUG
	#define USE_C_MATRIX_SOLVERS
 //USE_C_MATRIX_SOLVERS_AND_PYTHON_AND_COMPARE does the matrix math using the C and python versions and then compares the results for testing.
// Only usfull if USE_C_MATRIX_SOLVERS is also defined
//	#define USE_C_MATRIX_SOLVERS_AND_PYTHON_AND_COMPARE

//#define USE_EIGEN 1


#include "numpy1.h"
#include <stdio.h>
#include <stdlib.h> 
#include <cstring>			// for memcpy
#include <vector>
#include <math.h>
#include <float.h>

#ifdef USE_EIGEN
	#include <Eigen/Dense>
	using namespace Eigen;
#endif


#ifdef DO_DEBUG
	#define DEBUG(x) printf(x);
#else
	#define DEBUG(x) 
#endif

using namespace std;
#define PI 3.14159265358979323846


// HACK CHEAT GLOBAL just to have a place to put some statcs that can be printed in python.
//	TODO: make this nicer.. have a better stats reporting.. cython.
double * solverStatsGlobal = 0;


// This is an example of an exported variable
NUMPY1_API int nnumpy1=0;

// This is an example of an exported function.
NUMPY1_API int fnnumpy1(void)
{
	return 42;
}


typedef int (*matrix_solver_t)(double*, int, int, double);
typedef int (*no_arg_function_t)(int);

extern "C" {
	extern NUMPY1_API int basic_array_write( double* data, int len) {
		int i;
		printf("data = %p\n", (void*) data);
		for (i = 0; i < len; i++) {
			printf("data[%d] = %f\n", i, data[i]);
		}
		printf("len = %d\n", len);
		return len + 1;
	}




	
	typedef void*(*allocator_t)(int, int*, int);
	enum numpy_type { D, I, F};
	allocator_t numpy_allocator = NULL;
	int dim_temp1[1];
	int dim_temp2[2];
	int dim_temp3[3];

	extern NUMPY1_API void numpy_oshun_setup(allocator_t allocator) {
		numpy_allocator = allocator;
		//printf("numpyallocator set to = 0x%p\n", numpy_allocator);
	}

	extern void* numpy_malloc( numpy_type data_type, int dim, int dim1, int dim2, int dim3) {
		int* shape_array = 0;
		
		if(dim == 1) {
			dim_temp1[0] = dim1;
			shape_array = dim_temp1;

		} else if(dim == 2) {

			dim_temp2[0] = dim2;
			dim_temp2[1] = dim1;
			shape_array =  dim_temp2;

		} else if(dim == 3) {

			dim_temp3[0] = dim3;
			dim_temp3[1] = dim2;
			dim_temp3[2] = dim1;
			shape_array = dim_temp3;

		} else {
			return NULL;
		}		
		//printf("numpyallocater = 0x%p\n", numpy_allocator);
		void * memory_pointer = numpy_allocator(dim, shape_array, data_type);
		//printf("Mem was allocated! \n");
		//printf("data = 0x%p\n", memory_pointer);
		return memory_pointer;
	}

}

/*
dim1 = x, dim2=y, dim3=z

*/
template<class T>
class numpy_ref {
	public:
		int dim;
		int dim1; int dim2; int dim3; int dims[3];
		T* data;
		numpy_type numpy_data_type;

		numpy_ref(numpy_type numpy_type_data_type, int dim1) :						dim(1),dim1(dim1),dim2(0),
																					dim3(0),numpy_data_type(numpy_type_data_type) {
			data = (T*) numpy_malloc(numpy_data_type, 1, dim1,0,0);
		};
		numpy_ref(numpy_type numpy_type_data_type, int dim1, int dim2) :			dim(2),dim1(dim1),dim2(dim2),
																					dim3(0),numpy_data_type(numpy_type_data_type) {
			data = (T*) numpy_malloc(numpy_data_type,2,dim1,dim2,0);
		};
		numpy_ref(numpy_type numpy_type_data_type, int dim1, int dim2, int dim3) :	dim(3),dim1(dim1),dim2(dim2),
																					dim3(dim3),numpy_data_type(numpy_type_data_type) {
			data = (T*) numpy_malloc(numpy_data_type,3,dim1,dim2,dim3);
		};

		~numpy_ref() {
			// cause havoc!!!!!
			printf("Destruction and meyhem!!! \n");
			printf("We are deleting a Numpy Reference from c++.. This is dangerous and things will go wrong. It shouldn't happen.. All malloc/demalloc new/delete should be done from Python. \n");
			delete data;
		}

};


class numpy_functions_test {

	public:
		numpy_ref<double>* numpy_array1;

		void init() {
			numpy_array1= new numpy_ref<double>(D, 1000,200);
		}

		void add_one() {

			double * data = numpy_array1->data;
			if( numpy_array1->dim == 1) {
				for(int i=0;i<numpy_array1->dim1;i++) {
					data[i] += 1.0;
				}
			}
			if( numpy_array1->dim == 2) {
				int stride1 = numpy_array1->dim1;
				for(int y=0;y<numpy_array1->dim2;y++) {
					for(int x=0;x<numpy_array1->dim1;x++) {
						//printf("%f ", data[y*stride1 + x]);
						data[y*stride1 + x] += 1;
					}
					//printf("\n");
				}
				//printf("\n");
			}

		}

		void add_one_flat() {
			double * data = numpy_array1->data;
			if( numpy_array1->dim == 1) {
				for(int i=0;i<numpy_array1->dim1;i++) {
					data[i] += 1.0;
				}
			}
			if( numpy_array1->dim == 2) {

			}			
		}

};

numpy_functions_test npy_functs;

extern "C" {
	extern  NUMPY1_API void test_alloc() {
		printf("c++: I am gonna allocate a numpy array!!!");
		npy_functs.init();
		printf("c++: done!\n\n");
	}
	
	extern  NUMPY1_API void add_one() {
		npy_functs.add_one();
	}
	extern NUMPY1_API void add_one_to(double *data, int* shape, int dim ) {
			
			if( dim == 1) {
				for(int i=0;i<shape[0];i++) {
					data[i] += 1.0;
				}
			}
			// a numpy array of size(y,x) has the xdim in shape[1], ydim in shape[0]
			if( dim == 2) {
				int stride1 = shape[1];
				for(int y=0;y<shape[0];y++) {
					for(int x=0;x<shape[1];x++) {
						//printf("%f ", data[y*stride1 + x]);
						data[y*stride1 + x] += (y+1);
					}
					//printf("\n");
				}
				//printf("\n");
			}
	}
	

}
struct numpy_array {
	double* data;
	int dim;
	int shape[3];
};



/*
=====================================================================
Lifted from Numerical Recepies...

no consideration given for cahce coherence etc.. mainly, this is here as
	a resonable alternative to calling pyython if there it is not desierable to
	link to ATLAS or BLAS etc...
=====================================================================
*/

#define A(i,j) a[i*stride+j]


int ludcmp(	double * RESTRICT a, int stride, int n, int* RESTRICT indx, 
			double * RESTRICT temp, double * RESTRICT d) {
	double big,sum,dum, crap;
	int imax;
	
	*d = 1.0;
	imax = 0;
	for(int i=0;i<n;i++) {
		big = 0.0;
		//int linear_index = stride*i;
		// find the largest element in the row to use as a scale factor (and remeber it)
		for(int j=0;j<n;j++) {
			crap = A(i,j);
			if(crap < 0.0) crap *= -1.0;
			if(crap > big)  big =  crap;
		}
		// if all elements were zero (i.e. we have a dengerenate filth-bag matrix) then return an error code telling the row it occured upon.
		if(big == 0.0) return -1*(i+1);
		temp[i] = 1.0/big;
	}
	
	for(int j=0;j<n;j++) {
		for(int i=0;i<j;i++) {
			sum = A(i,j); //a[stride*i + j];
			for(int k=0;k<i;k++) sum -= A(i,k)*A(k,j); //a[stride*i + k]*a[stride*k+j];
			A(i,j) = sum;	//a[stride*i + j] = sum;
		}
		big=0.0;
		for(int i=j; i< n;i++) {
			sum = A(i,j); //a[stride*i + j];
			for(int k=0;k<j;k++) sum -= A(i,k)*A(k,j); //a[stride*i + k]*a[stride*k+j];
			A(i,j)=sum; //a[stride*i + j] = sum;
			dum = temp[i]*fabs( sum );
			if(dum >= big ) {
				big=dum;
				imax = i;
				
			}
		}
		
		if(j != imax) {
			for(int k=0;k<n;k++) {
				dum = A(imax,k); //a[stride*imax + k];
				A(imax,k) = A(j,k); //a[stride*imax + k] = a[stride*j +k];
				A(j,k)=dum; //a[stride*j +k] = dum;
			}
			*d = -(*d);
			temp[imax] = temp[j];
		}
		indx[j] = imax;
		
		// put the following line in to help wit singularities.. but pur marices should never even come close to needed this.
		// if(A(j,j) == 0.0) A(j,j) = TINY;

		if(j != n) {
			dum = 1.0 / (A(j,j)); //(a[stride*j + j]);
			for(int i=j+1;i<n;i++) A(i,j) *= dum; //a[stride*i + j] *= dum;
		}
	}

	return 0;
}



int lubksb( double* RESTRICT a, int stride, int n, int* RESTRICT indx, double* RESTRICT b) {
	int ii=-1; int ip;
	double sum;
	
	for(int i=0;i<n;i++) {
		ip = indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if(ii>-1) {
			for(int j=ii;j<=(i-1);j++){
				sum -= A(i,j)*b[j]; //a[stride*i + j]*b[j];
			}
		}
		else if(sum) ii = i;
		b[i] = sum;
	}
	for(int i=n-1;i>=0;i--) {
		sum=b[i];
		for(int j=i+1;j<n;j++) {
			sum -= A(i,j)*b[j]; //a[stride*i - j]*b[j];
		}
		b[i] = sum/ A(i,i); //a[stride*i + i];
	}
	return 0;
}

int ldu_nonfancy(	double* RESTRICT a, int stride, int n, 
					double* RESTRICT right_hand_side, int* RESTRICT indx, double * RESTRICT temp) {
	//a:		matrix to solvify
	//stride: 	usually the same a 'n'.. but it is the number of 'doubles' in physical memory for each row of 'a'
	//n: 		dimension of a (a must be square)
	//			sometimes the length of a row will be larger then the actual number of elements for memory alignment reasons)
	//right_hand_side:	The right hand side of the matrix equation. Overwirrent with the final answer.
	//indx:		an integer array that has length at least as large as 'n' used for bookkeeping by the solver
	//temp:		an array that has length at least as large as 'n' use for temporary storage by the solver.
	
	// solution returned in 'right_hand_side' (so 'right_hand_side' is over written)
	double d = 0.0;
	int err = ludcmp(a, stride, n, indx, temp, &d);
	if(err < 0) return -1;
	err = lubksb(a, stride, n, indx, right_hand_side);
	if(err < 0) return -2;
	return 0;
}

int tridiag__c( double* RESTRICT a, double* RESTRICT b, double* RESTRICT c, double* RESTRICT r, 
				double* RESTRICT u, double* RESTRICT temp, int length) {
	// a is lower diagonal.. a[0] is unused
	// b is diagonal
	// c is upper diagonal.  c[length-1] is unused.
	// r is thr right hand side.
	// temp must be temp array of at least length 'length'
	// u is solution...
	
	if(b[0] == 0.0) { 
		printf("\n b[0] = 0.0!!! \n\n");
		return -1;
	}
	
	double bet =  b[0];
	double bet_old = bet;
	
	u[0] = r[0] / ( bet );
	for(int j=1;j<length;j++) {
		temp[j] = c[j-1]/bet;
		bet_old = bet;
		bet=b[j]-a[j]*temp[j];
		if(bet == 0.0)  {
			printf("\n A zero pivot.. wierd.. porbably not diagonallly dominate.!!! \n\n");
			printf("j=%d, bet was %f, c[j-1]=%f, b[j]=%f a[j]=%f  a[j]*temp[j]=%f", j,bet_old,  c[j-1], b[j], a[j], a[j]*temp[j] );
			return -1;
		}
		u[j]=(r[j]-a[j]*u[j-1])/bet;
	}
	
	// now back subsititute
	for(int j=length-2;j>=0;j--) {
		u[j] -= temp[j+1]*u[j+1];
	}
	return 0;
}


// Copied (only slightly modified) from the 2D C++ versin of the code
bool Gauss_Seidel(	double* RESTRICT a, double* RESTRICT b, double* RESTRICT xk, 
					double* RESTRICT temp, int n, int stride) {
//-------------------------------------------------------------------
//   Fills solution into xk. The other matrices are not modified
//   The function returns "false" if the matrix A is not diagonally
//   dominant
//-------------------------------------------------------------------

	double     tol(1.0e-8);    // Tolerance for absolute error
	int        MAXiter(50);    // Maximum iteration allowed

//  The Matrices all have the right dimensions
//  -------------------------------------------------------------

//  Check if the matrix A is diagonally dominant
//  -------------------------------------------------------------
	for (int i(0); i < n; ++i){
		double rowi(0.0);
		for (int j(0); j < n; ++j){
		    rowi += fabs(A(i,j));
		}
		//if (!(rowi < 2.0*A(i,i))) return false;
		//printf("[%e %e] ", fabs(A(i,i)), fabs(rowi));
		if (!(fabs(rowi) < 2.0*fabs(A(i,i)))) return false;

		// invert the diagnoal.
		temp[i] = 1.0/A(i,i);
					//xk[i] = b[i];
		// copy over xk to xold
		//xold[i] = xk[i];
	}
;//  -------------------------------------------------------------


//  Calculate and invert the diagonal elements once and for all
//    valarray<double> invDIAG(A.dim1());
//	for (int i(0); i < invDIAG.size(); ++i) invDIAG[i] = 1.0 / A(i,i);

//	valarray< complex<double> > xold(xk);
	int iteration(0);         // used to count iterations
	int conv(0);              // used to test convergence
	double delta;


//  Start the iteration loop
	while ( (iteration++ < MAXiter) && (conv <n) ) {
			conv = 0;

            //xold = xk;
            for (int i(0); i < n; ++i){
                double sigma(0.0);    // Temporary sum
                for (int j(0); j < i; ++j){
                    sigma += A(i,j)*xk[j];   
                }
                for (int j(i+1); j < n; ++j){
                    sigma += A(i,j)*xk[j];   
                }
                delta = xk[i];
                xk[i] = temp[i] * (b[i] - sigma);
                delta -= xk[i];

                if ( fabs(delta) < (tol*fabs(xk[conv] + 20.0*DBL_MIN)) ) { 
                	++conv;
                }
            }

            // Calculate Dx = x_old - x_new
            //xold -= xk;

//          If the relative error < prescribed tolerance everywhere the method has converged
//          |Dx(i)| < t*|x(i)| + eps 
//            conv = 0;
//            while ( ( conv < b.size() ) && 
//                    ( abs(xold[conv]) < (tol*abs(xk[conv] + 20.0*DBL_MIN)) ) ){ 
//                ++conv;
//            } 

            //----> Output for testing
            //--------------------------------
            //cout << "iteration = " << iteration << "    ";
            //for (int i(0); i < b.size(); ++i){
            //    cout << xk[i] << "        ";
            //}
            //cout << "\n";
            //--------------------------------

        }

        //----> Output for testing
        //--------------------------------
        // cout << "Iterations = " << iteration-1  <<"\n";
        //for (int i(0); i < b.size(); ++i) {
        //    cout << "Error |Dx| = " << abs(xold[i]) 
        //         << ",    " 
        //         << "Tolerance * |x| = " << tol*abs(xk[i]) <<"\n";
        //}
        //--------------------------------
        //printf("number of itteratins used %d anf got %d converences \n", iteration, conv);

        double its = (double) iteration;
        if(solverStatsGlobal[0] < iteration) solverStatsGlobal[0] = iteration;

        return true;
    }
#undef A

/*
=====================================================================
 E field on the Dsitibution function!
=====================================================================
*/

struct axis {
	double * values;
	double dx;
	double dp;
	int num_points;

	void init(double* vals, int num_points) {
		this->values = vals;
		this->num_points = num_points;
		dx = values[1]-values[0];
		dp = dx;
	}
};

struct Field1D {
	double * data;
	int shape[3];
	int boundry_cells[3][2];
	int dim;

	axis spatial_axis__x;
	vector<axis> spatial_axes;

	Field1D() {
		dim = 1;
		shape[0] = 0; shape[1] = 0; shape[2] = 0;
		boundry_cells[0][1] = 0;boundry_cells[0][1] = 0;
		boundry_cells[1][1] = 0;boundry_cells[1][1] = 0;
		boundry_cells[2][1] = 0;boundry_cells[2][1] = 0;
		spatial_axes.resize(3);

	};
};

struct SphericalHarmonic {
	int l;
	int m;
	double *momentum_projection;
	axis* momentum_axis;

	void init(int l, int m, axis* momentum_axis) {
		this->l = l;
		this->m = m;
		this->momentum_axis = momentum_axis;
	}
};


struct SpatialCell {
	vector<SphericalHarmonic> harmonics;
	int num_harmonics;
	int* harmonic_map;
	int num_l;
	int num_m;
	

	void init( int num_l, int num_m, axis* momentum_axis) {
		int index = 0;
		num_harmonics = 0;
		harmonic_map = new int[num_l*num_m];
		for(int l=0;l<num_l;l++) {
			for(int m=0;m<num_m;m++) {
				if(m > l) break;
				num_harmonics += 1;
			}
		}

		harmonics.resize( num_harmonics );

		for(int l=0;l<num_l;l++) {
			for(int m=0;m<num_m;m++) {
				if(m > l) break;
				harmonics[index].init( l, m, momentum_axis);
				harmonic_map[l*num_m + m] = index;
				index++;
			}
		}
		this->num_m = num_m;
		this->num_l = num_l;
	}

	SphericalHarmonic* get_harmonic(int l, int m) {
		int idx = harmonic_map[l*num_m + m];
		return &harmonics[idx];
	}
	
};

struct DistributionFunction {
	int num_total_cells[3];
	int num_plain_cells[3];
	int num_boundry_cells[3][2];

	vector<SpatialCell> cells;
	axis* momentum_axis;
	int num_l;
	int num_m;

	// this is the 'new' data layout.. where all data is in one large
	//	contigous buffer. 
	// TODO: take out old style data layout
	// TODO: investigate optimzing alignment requirements.. 
	//
	// The Dist Functions's data buffer ihas 3 indicies (in 1D)
	//		and has the shape: [ num_total_cells[0], num_l, num_p]
	//		data.shape[0] = num_total_cells[0], data.shape[1]=shape, data.shape=num_p
	//
	//	so indexing a data element data[k, j, x] looks like:
	//			data_element = data[k*num_l*num_p + j*num_p + x];
	double *data;
	
	int linear_index_for_harmonic(int l, int m) {
		return l*num_m +m;
	}

	void init(	int num_l, int num_m, axis* momentum_axis, 
				int* num_total_cells, int* num_plain_cells, int* num_boundry_cells, 
				double* data) {
		
		cells.resize( num_total_cells[0] );
		this->momentum_axis = momentum_axis;
		this->num_l = num_l;
		this->num_m = num_m;
		this->data = data;
		
		for(int i=0;i<3;i++) {
			this->num_total_cells[i] = num_total_cells[i];
			this->num_plain_cells[i] = num_plain_cells[i];
			this->num_boundry_cells[i][0] = num_boundry_cells[2*i];
			this->num_boundry_cells[i][1] = num_boundry_cells[2*i+1];
		}

		for(int i=0;i<num_total_cells[0];i++) {
			cells[i].init(num_l, num_m, momentum_axis);
		}
	}
	
};


// DEBUG vars to help test/validate
// TODO: remove these..
// helpful easy debugging aids 
int    debug__current_sim_timestep;
double debug__current_sim_time;
double debug__sim_current_dt;


struct species {

	int num_l;
	int num_m;
	axis momentum_axis;
	DistributionFunction F;
	
	void init( int num_l, int num_m, double* momentum_axis_values, int momentum_axis_num_points,
		int* F_num_total_cells, int* num_plain_cells, int* num_boundry_cells, double* F_data) {

		this->num_l = num_l;
		this->num_m = num_m;
		momentum_axis.init( momentum_axis_values, momentum_axis_num_points);
		
		F.init( num_l, num_m, &momentum_axis, F_num_total_cells, num_plain_cells, num_boundry_cells, F_data);

	// Debugging randomly shove it here since is always called.
	// TODO: take test out
	debug__current_sim_timestep = 0;
	debug__current_sim_time		= 0.0; 
	debug__sim_current_dt		= 0.0;
		
	}
	
};

// function pointer type fot the matrix solving functions....

vector<axis> SpatialAxes;
vector<species> Species;
Field1D E_field;

/*
---------------------- Fokker Plank ------------------
*/
struct fokker_plank_explicit_advance {
	numpy_array vr;
	numpy_array U4;
	numpy_array U4m1;
	numpy_array U2;
	numpy_array U2m1;
	numpy_array U1;
	numpy_array U1m1;
	numpy_array J1;
	numpy_array I4;
	numpy_array U3;
	numpy_array Qn;
	numpy_array Pn;
	numpy_array I2;

	int num_pr_cells;
	int shape_vr[3];

	// precalculations for the lower |p| cells
	double		G_constant_1;
	numpy_array G_constant_vr3;
	numpy_array G_constant_vr5;
	numpy_array G_constant_vr5_2;
	numpy_array G_constant_vr7;
	double density_np;
	int NB;
	double c_kpre;
	
	int num_subcycling_steps;
	// used in RK4
	// TODO: deallcate these when finished..
	double* F0; double* F1; double* Fh;
	
	// what species are we associated with?
	// TODO PASS THIS IN TAKE OUT HARD CODING.
	int species_index; species my_species; DistributionFunction F; double* F_data;
	
	

	void init(	double *data_vr, long long* shape_vr_long,
				double *data_U4, double *data_U4m1,
				double *data_U2, double *data_U2m1,
				double *data_U1, double *data_U1m1,
				double *data_J1, double *data_I4,
				double *data_U3, double *data_Qn,
				double *data_Pn, double* data_I2, 
				double data_G_constant_1, double* data_G_constant_vr3,
				double *data_G_constant_vr5, double* data_G_constant_vr5_2,
				double *data_G_constant_vr7, double data_density_np, int NB, double c_kpre, 
				int num_subcycling_steps, int species_index) {
		
		shape_vr[0] = (int) shape_vr_long[0];
		species_index = 0;
		my_species = Species[0]; F = my_species.F; F_data = F.data;
		
		
		vr.data = data_vr; vr.shape[0] = shape_vr[0]; vr.dim = 1;
		U4.data = data_U4; U4.shape[0] = shape_vr[0]; U4.dim = 1;
		U4m1.data = data_U4m1; U4m1.shape[0] = shape_vr[0]; U4m1.dim = 1;
		U2.data = data_U2; U2.shape[0] = shape_vr[0]; U2.dim = 1;
		U2m1.data = data_U2m1; U2m1.shape[0] = shape_vr[0]; U2m1.dim = 1;
		U1.data = data_U1; U1.shape[0] = shape_vr[0]; U1.dim = 1;
		U1m1.data = data_U1m1; U1m1.shape[0] = shape_vr[0]; U1m1.dim = 1;
		J1.data = data_J1; J1.shape[0] = shape_vr[0]; J1.dim = 1;
		I4.data = data_I4; I4.shape[0] = shape_vr[0]; I4.dim = 1;
		U3.data = data_U3; U3.shape[0] = shape_vr[0]; U3.dim = 1;
		Qn.data = data_Qn; Qn.shape[0] = shape_vr[0]; Qn.dim = 1;
		Pn.data = data_Pn; Pn.shape[0] = shape_vr[0]; Pn.dim = 1;
		I2.data = data_I2; I2.shape[0] = shape_vr[0]; I2.dim = 1;
		
		G_constant_1 = data_G_constant_1;
		G_constant_vr3.data = data_G_constant_vr3; G_constant_vr3.shape[0] = shape_vr[0]; G_constant_vr3.dim = 1;
		G_constant_vr5.data = data_G_constant_vr5; G_constant_vr5.shape[0] = shape_vr[0]; G_constant_vr5.dim = 1;
		G_constant_vr5_2.data = data_G_constant_vr5_2; G_constant_vr5_2.shape[0] = shape_vr[0]; G_constant_vr5_2.dim = 1;
		G_constant_vr7.data = data_G_constant_vr7; G_constant_vr7.shape[0] = shape_vr[0]; G_constant_vr7.dim = 1;
		
		density_np = data_density_np;
		num_pr_cells = shape_vr[0];
		this->NB = NB;
		this->c_kpre = c_kpre;
		this->num_subcycling_steps = num_subcycling_steps;
		
		int size_dist_function_buffer = F.num_plain_cells[0]*num_pr_cells*F.num_l;
		this->F0 = new double[size_dist_function_buffer];
		this->F1 = new double[size_dist_function_buffer];
		this->Fh = new double[size_dist_function_buffer];
	}

	double fokker_plank_explicit_G( int n, double* fin) {
		double f00 = (fin[0] - fin[1]*G_constant_1) / (1.0 - G_constant_1);
		double i2s = f00*G_constant_vr3.data[n]/3.0 + (fin[1]-f00)*G_constant_vr5_2.data[n];
		double i4s = f00*G_constant_vr5.data[n] + (fin[1]-f00)*G_constant_vr7.data[n];
		return fin[n]*i4s + (G_constant_vr3.data[n]*fin[n]-3.0*i2s) * J1.data[n];
	}

	double LOGee( double ne, double Te) {
		if( ne > 1e-9) {
			Te /= (3.0*ne);
			Te *= 511000.0;
			ne *= density_np;
			Te = log(Te);
			ne = log(ne);
			double lnee = 23.5 - 0.5*ne + 1.25*Te - sqrt(0.00001+0.0625*(Te-2.0)*(Te-2.0));
			if (lnee > 2.0) {
				return lnee;
			}
		}
		return 2.0;
	}

	double ZLOGei( double ne, double Te, double zeta) {
		if( ne > 1e-9 ) {
			Te /= (3.0*ne);
			Te *= 511000.0;
			ne *= density_np ;
			Te = log(Te);
			ne = log(ne);
			double lnei = 24.0 - 0.5*ne + Te;
			if( lnei > 2.0 ) return lnei*zeta;
		}
		return 2.0*zeta;
	}
	
	void fokker_plank_explicit_slope_calc(double *fin, double *fh) {

		I4.data[0] = 0.0;
		for (int n=1; n < num_pr_cells; ++n) {
			I4.data[n]  = U4.data[n]*fin[n]+U4m1.data[n]*fin[n-1]; 
			I4.data[n] += I4.data[n-1];
		}

		I2.data[0] = 0.0;
		for (int n(1); n < num_pr_cells; ++n) {
			I2.data[n]  = U2.data[n]*fin[n]+U2m1.data[n]*fin[n-1]; 
			I2.data[n] += I2.data[n-1];
		}

		J1.data[num_pr_cells-1] = 0;
		for (int n(num_pr_cells-2); n > -1; --n) {
			J1.data[n]  = U1.data[n+1]*fin[n+1]+U1m1.data[n+1]*fin[n]; 
			J1.data[n] += J1.data[n+1];
		}
		
		double Ln_ee = LOGee(4.0*PI*I2.data[num_pr_cells-1],4.0*PI*I4.data[num_pr_cells-1]) * c_kpre;
		
		for (int n(0); n < num_pr_cells; ++n) {
			I2.data[n] *= J1.data[n];				// J1(n) * I2(n)
			I2.data[n] *= -3.0;						// -3 * J1(n) * I2(n)
			I4.data[n] += U3.data[n] * J1.data[n];	// I4(n) + u_n^3 * J1(n) 
			I4.data[n] *= fin[n];					// fn * I4(n) + u_n^3 * fn * J1(n)
			I4.data[n] += I2.data[n];				// Gn = fn * I4(n) + u_n^3 * fn * J1(n) - 3 * J1(n) * I2(n)
		}

		for (int n(0); n < NB; ++n) { 
			I4.data[n] = fokker_plank_explicit_G(n,fin);
		}

		// Find -DG
		fh[0]  = (-1.0)*I4.data[0];
		for (int n(0); n < num_pr_cells-1; ++n) I4.data[n] -= I4.data[n+1];

		// Find -DG/(vDv)
		fh[0] *= 2.0/ (vr.data[0]*vr.data[0]);
		for (int n(0); n < num_pr_cells-1; ++n) I4.data[n] *= Pn.data[n];

		// Find DDG/(vDv)
		fh[0] -= I4.data[0];  
		for (int n(0); n < num_pr_cells-1; ++n) I4.data[n] -= I4.data[n+1];

		// Find DDG/(v^3*DDv)
		fh[0] *= Qn.data[0]*Ln_ee;
		for (int n(0); n < num_pr_cells-1; ++n) fh[n+1] = Qn.data[n+1]*I4.data[n]*Ln_ee;
		fh[num_pr_cells-1] *= Ln_ee;
		
		// Normalize
		//fh *=  c_kpre *  Ln_ee;
	}
	
	// 
	void RK4_for_fokker_plank( double* F_in, const double subcycle_dt, const int data_size, const int data_size__bytes) {
		// F_in is part of the offical state.. so any changes made to it will be used by rest of code..
		//	so we gotta be carefull with it.
		 
		
		memcpy( F0, F_in, data_size__bytes );
		memcpy( F1, F_in, data_size__bytes );
		
		// this function evalution must gaurentee that the input REMAINS UNCHANGED.
		//	Fh will have the slope data for RK4, F_in must be returned unchanged.
		fokker_plank_explicit_slope_calc( F1, Fh);
		
		// f1 = f1 + (h/2)*fh
		for(int i=0;i<data_size;i++) {
			Fh[i] *= 0.5*subcycle_dt;
			F1[i] += Fh[i];
			Fh[i] *= (1.0/3.0);
			F_in[i] += Fh[i];
		}
		
		//------- Step 2
		fokker_plank_explicit_slope_calc( F1, Fh);
		
		for(int i=0;i<data_size;i++) {
			F1[i] = F0[i];
			// f1 = f0 + (h/2)*fh
			Fh[i] *= (0.5*subcycle_dt);
			F1[i] += Fh[i];
			// F = F + (h/3)*Fh --> F + h/6*k1+h/2*k2
			Fh[i] *= (2.0/3.0);
			F_in[i] += Fh[i];
		}
		
		//------- Step 3
		fokker_plank_explicit_slope_calc( F1, Fh);
		for(int i=0;i<data_size;i++) {
			// f1 = f0 + h*Fh
			Fh[i] *= subcycle_dt;
			F0[i] += Fh[i];
			// F = F + (h/3)*Fh --> F+ h/6*k1+h/3*k2+h/3*k3
			Fh[i] *= (1.0/3.0);
			F_in[i] += Fh[i];
		}
		
		//------- Step 4
		fokker_plank_explicit_slope_calc( F0, Fh);
		
		for(int i=0;i<data_size;i++) {
			// F = F + (h/6)*Fh --> F+ h/6*k1+h/3*k2+h/3*k3*+1/6*k4
			Fh[i] *= (subcycle_dt/6.0);
			F_in[i] += Fh[i];
		}
		
		return;
	}
	
	// main entry point for (L=0, m=0) collsions.. results will be written
	//		directly into the input distribution function.
	void calc_f00( int n_step, double time, double dt) {
		// loop over the useful (non-guard) cells.. 
		// (Our Fokker-Plank calculation is completely local in space and is strictly
		//		an intra-spatial-cell calculation.. so no guard cell are required )
		const int num_cells__x = F.num_plain_cells[0];
		const int F_data__cell_stride = F.num_l * num_pr_cells;
		const double subcycle_dt = dt / (double) num_subcycling_steps;
		const int data_size = num_pr_cells; 
		const int data_size__bytes = num_pr_cells*sizeof(double);
		
		int starting_F_data_location = F.num_boundry_cells[0][0] * F_data__cell_stride;
		
		//printf("----------entering explicit C++ collison loop --------------- subcyle %d with dt=%f\n", num_subcycling_steps, subcycle_dt);
		for(int x=0; x < num_cells__x; x++) {
			// TODO: use higher abstracted addressiing
			double * momentum_projection = &( F_data[starting_F_data_location] );
			
			for(int n=0; n < num_subcycling_steps; n++) {
				RK4_for_fokker_plank(momentum_projection , subcycle_dt, data_size, data_size__bytes);
			}
			starting_F_data_location += F_data__cell_stride;
		}
	}
	
};


struct fokker_plank_implicit_advance {

	fokker_plank_explicit_advance* explict_code;
	
	numpy_array vr;
	numpy_array vr3;
	numpy_array ddf0;
	numpy_array df0;
	numpy_array Scattering_Term;
	numpy_array Alpha;
	numpy_array AlphaTri;
	numpy_array J1m;
	numpy_array TriI1;
	numpy_array TriI2;
	numpy_array IvDnDm1;
	numpy_array IvDnDp1;
	numpy_array Ivsq2Dn;
	numpy_array I0;
	int num_pr_cells;
	size_t alpha_matrix_size__bytes;
	
	double* solverStats;

	double* _LOGee;
	double* _ZLOGei;
	double zeta;
	double f00_factor;
	double I0_density;
	int is_tridiagonal;
	int if_implicit1D;
	int shape_vr[3];
	
	//int(*tridiag_matrix_solver)(double*);
	//int(*full_matrix_solver)(double*);
	matrix_solver_t tridiag_matrix_solver;
	matrix_solver_t full_matrix_solver;
	no_arg_function_t debug_matrix_solver;
	
	// full C++ solver....
	bool use_C_matrix_solvers;
	int* solver_indx; double* solver_temp;
	double* solver_tridiag__a; double* solver_tridiag__b; double* solver_tridiag__c; double *solver_tridiag__solution;
	double *solver_scratch;
	#ifdef USE_EIGEN
	MatrixXd _Alpha__eigen;
	#endif
	
	// what species are we associated with?
	// TODO PASS THIS IN TAKE OUT HARD CODING.
	int species_index; species my_species; DistributionFunction F; double* F_data;	
	
	double kpre;
	#define  four_pi 4.0*3.141592654
	#define eight_pi 8.0*3.141592654
	
	
	void init(	double *data_vr, long long * shape_vr_long,
				double *data_vr3,
				double *data_df0,
				double *data_ddf0, 
				double *data_Scattering_Term, 
				double *data_Alpha,
				double *data_AlphaTri,
				double *data_J1m,
				double *data_TriI1,
				double *data_TriI2,
				double* data_IvDnDm1,
				double* data_IvDnDp1,
				double* data_Ivsq2Dn,
				double* data_I0,
				double kpre, double zeta, double f00_factor,
				int is_tridiagonal, int if_implicit1D, int species_index,
				double* _LOGee, double* _ZLOGei,
				matrix_solver_t tridiag_matrix_solver1,
				matrix_solver_t full_matrix_solver1,
				no_arg_function_t debug_matrix_solver,
				double* solverStats) {

				//int(*tridiag_matrix_solver1)(double*),
				//int(*full_matrix_solver1   )(double*)  ) {
		
		
		this->solverStats = solverStats; solverStatsGlobal=solverStats;
		this->solverStats[0] = 0.0;this->solverStats[1] = 0.0;this->solverStats[2] = 0.0;this->solverStats[3] = 0.0;

		//is_tridiagonal = 1;
		//if_implicit1D = 1;
		//species_index = 0;
		shape_vr[0] = (int) shape_vr_long[0];
		my_species = Species[species_index]; F = my_species.F; F_data = F.data;
		
		vr.data = data_vr; vr.shape[0] = shape_vr[0]; vr.dim = 1;
		vr3.data = data_vr3; vr3.shape[0] = shape_vr[0]; vr3.dim = 1;
		df0.data = data_df0; df0.shape[0] = shape_vr[0]; df0.dim = 1;
		ddf0.data = data_ddf0; ddf0.shape[0] = shape_vr[0]; ddf0.dim = 1;
		Alpha.data = data_Alpha; Alpha.shape[0] = shape_vr[0]; Alpha.shape[1] = shape_vr[0]; Alpha.dim = 2;
		AlphaTri.data = data_AlphaTri; AlphaTri.shape[0] = shape_vr[0]; AlphaTri.shape[1] = shape_vr[0]; AlphaTri.dim = 2;
		Scattering_Term.data = data_Scattering_Term; Scattering_Term.shape[0] = shape_vr[0]; Scattering_Term.dim = 1; 
		
		J1m.data = data_J1m; J1m.shape[0] = shape_vr[0]; J1m.dim = 1;
		I0.data =  data_I0;   I0.shape[0] = shape_vr[0];  I0.dim = 1;
		
		TriI1.data = data_TriI1; TriI1.shape[0] = shape_vr[0]; TriI1.dim = 1;
		TriI2.data = data_TriI2; TriI2.shape[0] = shape_vr[0]; TriI2.dim = 1;
		
		IvDnDm1.data = data_IvDnDm1; IvDnDm1.shape[0] = shape_vr[0]; IvDnDm1.dim = 1;
		IvDnDp1.data = data_IvDnDp1; IvDnDp1.shape[0] = shape_vr[0]; IvDnDp1.dim = 1;
		Ivsq2Dn.data = data_Ivsq2Dn; Ivsq2Dn.shape[0] = shape_vr[0]; Ivsq2Dn.dim = 1;
		
		num_pr_cells = shape_vr[0];
		this->kpre = kpre;
		this->zeta = zeta;
		this->f00_factor = f00_factor;
		this->alpha_matrix_size__bytes = sizeof(double)*Alpha.shape[0]*Alpha.shape[1];
		
		this->is_tridiagonal = is_tridiagonal;
		this->if_implicit1D = if_implicit1D;
		this->species_index = species_index;
		
		// these are funciton pointer to supported Matrix libaries.
		//this->tridiag_matrix_solver = tridiag_matrix_solver;
		//this->full_matrix_solver = full_matrix_solver;
		this->tridiag_matrix_solver = tridiag_matrix_solver1;
		this->full_matrix_solver = full_matrix_solver1;
		this->debug_matrix_solver = debug_matrix_solver;

		
		// TOTOD: REMOVE THIS... ITS A BIT CONVOLUTED.
		this->_LOGee = _LOGee;
		this->_ZLOGei = _ZLOGei;
		
		#ifdef USE_C_MATRIX_SOLVERS
			this->use_C_matrix_solvers = true;
		#else
			this->use_C_matrix_solvers = false;
		#endif

		if(this->use_C_matrix_solvers) {
			printf("USING C-ONLY Solvers...\n\n");
			// setup space for the cheapo C solvers...
			//if(this->is_tridiagonal == 0) {
				solver_indx = new int[num_pr_cells];
				solver_temp = new double[num_pr_cells];
				solver_scratch = new double[num_pr_cells];
			//} else {
				solver_temp = new double[num_pr_cells];
				solver_tridiag__a = new double[num_pr_cells];
				solver_tridiag__b = new double[num_pr_cells];
				solver_tridiag__c = new double[num_pr_cells];
				solver_tridiag__solution = new double[num_pr_cells];

				// the following cells are not used by the solver.. fo set them to an silly, easy to spot
				// value for debugging purposes.
				solver_tridiag__a[0] = -666.0;
				solver_tridiag__c[num_pr_cells-1] = -666.0;
			//}
		}
	}
	
	void reset_coeff__cpp( double * fin, double Dt) {
		
		const int num_pr_cellsm1 = num_pr_cells -1; 
//printf("IN reset_coeff__cpp c++\n");
		// INTEGRALS
		// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
		// Temperature integral I2 = 4*pi / (v^2) * int_0^v f(u)*u^4du]
		
		numpy_array& I2 = explict_code->I2; I2.data[0] = 0.0;
		numpy_array& U1 = explict_code->U1;
		numpy_array& U1m1 = explict_code->U1m1;
		numpy_array& U2 = explict_code->U2;
		numpy_array& U2m1 = explict_code->U2m1;
		numpy_array& U4 = explict_code->U4;
		numpy_array& U4m1 = explict_code->U4m1;
		//numpy_array& I0 = explict_code->I0;
		
		for(int k=1; k<num_pr_cells;k++) {
			I2.data[k] = U4.data[k]*fin[k]+U4m1.data[k]*fin[k-1]; 
			I2.data[k] *= four_pi;
			I2.data[k] += I2.data[k-1];
//if(print_output1) {
	//printf("k=%d: u4[k]=%e fin[k]=%e,  fin[k-1]=%e, U4m1[k]=%e I2.data[k]=%e\n", k,  U4.data[k], fin[k], fin[k-1], U4m1.data[k], I2.data[k]);
//}
		}
		double I2_temperature = I2.data[num_pr_cells - 1];
		for (int k=0; k<num_pr_cells;k++) {
			I2.data[k] /=  (vr.data[k]*vr.data[k]);
		}
		
		// Density integral I0 = 4*pi*int_0^v f(u)*u^2du 
		I0.data[0] = 0.0;
		for(int k=1; k<num_pr_cells;k++) {
			I0.data[k]  = U2.data[k]*fin[k]+U2m1.data[k]*fin[k-1];
			I0.data[k] *= four_pi;
			I0.data[k] += I0.data[k-1];
		}
		I0_density = I0.data[num_pr_cells - 1];
		
		// Integral J_(-1) = 4*pi * v * int_0^v f(u)*u^4du
		J1m.data[ num_pr_cells - 1] = 0.0;
		//for(int k=(num_pr_cells-2); k<=0 ;k--) {
		for(int k=(num_pr_cells-2); k>=0 ;k--) {
		  J1m.data[k]  = U1.data[k+1]*fin[k+1]+U1m1.data[k+1]*fin[k]; 
		  J1m.data[k] += J1m.data[k+1];
		}
		for (int k=0; k<num_pr_cells;k++) {
			J1m.data[k] *= four_pi * vr.data[k];
		}
		
		// COULOMB LOGARITHMS
		// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		
		// these values are not shiped to exxternal functions anymore..
		//	And it didnt work anyway becasue Ctypes makes a new pointer object very function call
		//	TODO: So take the pointers to these Float python variables out of the init routine and dont store them
		//*_LOGee  = explict_code->LOGee( I0_density, I2_temperature );
		//*_ZLOGei = explict_code->ZLOGei( I0_density, I2_temperature, zeta );
		
/*
if(print_output1 && debug__current_sim_timestep < 10 ) {
	printf("\t REST_COEFF. We at step %d logee: %e (cpp: %e python: %e)\n", debug__current_sim_timestep, *_LOGee, *_ZLOGei, *(this->_LOGee));
	printf("\tele: I0_density=%e, I2_temperature=%e \n", I0_density, I2_temperature);
	printf("\tion: I0_density=%e, I2_temperature=%e,zeta=%e \n", I0_density, I2_temperature, zeta);
}
*/
		
		const double  LOGee_value = explict_code->LOGee( I0_density, I2_temperature );
		const double ZLOGei_value = explict_code->ZLOGei( I0_density, I2_temperature, zeta );
		//printf( "VALUE _ZLOGei %f  %f \n", ZLOGei_value, fred );

		// BASIC INTEGRALS FOR THE TRIDIAGONAL PART
		//   and SCATTERING TERM
		// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		
		
		// Take care of 0th element...
		Scattering_Term.data[0] = (ZLOGei_value * I0_density); 
		Scattering_Term.data[0] *= (kpre*Dt/vr3.data[0]);
		
		for(int i=1; i<num_pr_cells;i++) {
			double temp =   I0.data[i] + (2.0*J1m.data[i] - I2.data[i]) / 3.0;
			TriI2.data[i] = ( I2.data[i] + J1m.data[i] ) / 3.0;             		// ( I2 + J_{-1} ) / 3
			TriI1.data[i] = temp; 													// (-I2 + 2*J_{-1} + 3*I0) / 3
			Scattering_Term.data[i] = temp*LOGee_value;
			Scattering_Term.data[i] += (ZLOGei_value*I0_density);
			Scattering_Term.data[i] *= (kpre*Dt/vr3.data[i]);
		}
		
		
		const int alpha_stride = Alpha.shape[1];
		// Todo: is this zerosing necessary?]
		double factor = (-1.0) * LOGee_value * kpre * Dt;
		memset( AlphaTri.data, 0.0, alpha_matrix_size__bytes);
		
		AlphaTri.data[0] = eight_pi*fin[0]*factor;
		
		for(int i=1;i<num_pr_cellsm1;i++) {
			AlphaTri.data[i*alpha_stride+i  ] =  eight_pi*fin[i] - TriI2.data[i] * (IvDnDm1.data[i] + IvDnDp1.data[i]);
			AlphaTri.data[i*alpha_stride+i  ] *= factor;
			AlphaTri.data[i*alpha_stride+i-1] =  TriI2.data[i] * IvDnDm1.data[i] - TriI1.data[i] * Ivsq2Dn.data[i];
			AlphaTri.data[i*alpha_stride+i-1] *= factor;
			AlphaTri.data[i*alpha_stride+i+1] =  TriI2.data[i] * IvDnDp1.data[i] + TriI1.data[i] * Ivsq2Dn.data[i];
			AlphaTri.data[i*alpha_stride+i+1] *= factor;
			
		}
		
		
		// calculate the derivative
		for(int n=1;n<num_pr_cellsm1;n++) {
			df0.data[n]  = fin[n+1]-fin[n-1];
			df0.data[n] /= vr.data[n+1]-vr.data[n-1];
		}
		
		// Evaluate the second derivative
		// -df/dv_{n-1/2}
		for(int n=1;n<num_pr_cells;n++) {
			ddf0.data[n] = (fin[n-1]-fin[n]);
			ddf0.data[n]  /= (vr.data[n] - vr.data[n-1]);
		}
		
		// D(df/dv)/Dv
		for(int n=1;n<num_pr_cellsm1;n++) {
			ddf0.data[n] -= ddf0.data[n+1];
			ddf0.data[n]  /= 0.5*(vr.data[n+1]-vr.data[n-1]);
		}
		
		// Calculate zeroth cell
		double f00 =  ( fin[0] - f00_factor*fin[1] ) / ( 1.0 -  f00_factor );
		ddf0.data[0] = 2.0 * (fin[1] - f00) / (vr.data[1]*vr.data[1]);
		df0.data[0] = ddf0.data[0] * vr.data[0];
		
		// Calculate 1/(2v)*(d^2f)/(dv^2),  1/v^2*df/dv
		for(int n=0;n<num_pr_cellsm1;n++) {
			 df0.data[n] /= (vr.data[n]*vr.data[n]);
			ddf0.data[n] /= (2.0  *vr.data[n]);
		}
	}
	
	bool print_output1;
	//print_output1 = true;
	
	// TODO: take out the __LOGee.. in cpp mode we calulate in cpp and it is saved as a member variable..
	//		and its not used by the python card in any other way at this point...
	void implicit_advance(double * fin, int LL, double LOGee_from_python, double Dt, int current_x) {
		//printf("IN IMPLICIT c++\n");
		
		// if the member variable '_LOGee' is positive, then the cpp code caluclated this value
		//	in 'reset_coeff__cpp'. If not, the the python code calcualted '_LOGee' so use that value.
		// TODO: take this out.. it legacy for times when part of the FP operator was
		//		in python and part in C++.. So this wont be needed once the C++ and Fortran  FP
		//		implementation are sealed.
		double LOGee;
		if ( LOGee_from_python > 0.0) {
			LOGee = LOGee_from_python;
		} else {
			LOGee = *(this->_LOGee);
		}
		
//if(print_output1 && debug__current_sim_timestep < 10 ) {
//printf("\t IMPlicit adv step %d logee: %e (cpp: %e python: %e)\n", debug__current_sim_timestep, LOGee, LOGee_from_python, *(this->_LOGee));
//}
		
		const int num_pr_cellsm1 = num_pr_cells -1; //const int num_pr_cells1 = num_pr_cells;
		const double factor1 = (-1.0) * (*_LOGee) * kpre * Dt;
		const int alpha_stride = Alpha.shape[1];
		
		// the Alpha_Tri (tridagonal) array is reused, but the matrix solver is decructive..
		//	so we gotta copy Alpha_Tri in the buffer Alpha (Alpha is allowed to be killed by matrix solver)
		memcpy( Alpha.data, AlphaTri.data, alpha_matrix_size__bytes);
		
		#ifdef USE_EIGEN
		_Alpha__eigen
		#endif
		
		if(LL > 1) Alpha.data[0] = 0.0;
		
		
		if( is_tridiagonal == 0) {
			double A1 = (LL+1.0)*(LL+2.0) / ((2.0*LL+1.0)*(2.0*LL+3.0));
			double A2 = (-1.0) *(LL-1.0)* LL      / ((2.0*LL+1.0)*(2.0*LL-1.0));
			double B1 = (-1.0) *( 0.5 *LL*(LL+1.0) +(LL+1.0) ) / ((2.0*LL+1.0)*(2.0*LL+3.0));
			double B2 = (       (-0.5)*LL*(LL+1.0) +(LL+2.0) ) / ((2.0*LL+1.0)*(2.0*LL+3.0));
			double B3 = ( 0.5 *LL*(LL+1.0) +(LL-1.0) ) / ((2.0*LL+1.0)*(2.0*LL-1.0));
			double B4 = ( 0.5 *LL*(LL+1.0) - LL      ) / ((2.0*LL+1.0)*(2.0*LL-1.0));
			
			for (int i=0; i< num_pr_cellsm1; i++) {
				double t1 = A1*ddf0.data[i] + B1*df0.data[i];
				t1 *= factor1;
				double t2 = A1*ddf0.data[i] + B2*df0.data[i];
				t2 *= factor1;
				double t3 = A2*ddf0.data[i] + B3*df0.data[i];
				t3 *= factor1;
				double t4 = A2*ddf0.data[i] + B4*df0.data[i];
				t4 *= factor1;
				
				Alpha.data[i*alpha_stride] += t1 * ( 2.0*PI* pow(vr.data[0]/vr.data[i],LL+2)*vr.data[0]*vr.data[0]*(vr.data[1]-vr.data[0]) );
				Alpha.data[i*alpha_stride] += t3 * ( 2.0*PI*pow(vr.data[0]/vr.data[i],LL)  *vr.data[0]* vr.data[0]*(vr.data[1]-vr.data[0]) );

				for(int j=1;j<i;j++) {
					Alpha.data[i*alpha_stride+j] += t1 * ( 2.0*PI*pow(vr.data[j]/vr.data[i],LL+2)*vr.data[j]*vr.data[j]*(vr.data[j+1]-vr.data[j-1]) );
					Alpha.data[i*alpha_stride+j] += t3 * ( 2.0*PI*pow(vr.data[j]/vr.data[i],LL)  *vr.data[j]*vr.data[j]*(vr.data[j+1]-vr.data[j-1]) );

				}
				
				Alpha.data[i*alpha_stride+i] += t1 * ( 2.0*PI*vr.data[i]*vr.data[i]*(vr.data[i]-vr.data[i-1]) );
				Alpha.data[i*alpha_stride+i] += t3 * ( 2.0*PI*vr.data[i]*vr.data[i]*(vr.data[i]-vr.data[i-1]) );

				Alpha.data[i*alpha_stride+i] += t2 * ( 2.0*PI*vr.data[i]*vr.data[i]*(vr.data[i+1]-vr.data[i]) );
				Alpha.data[i*alpha_stride+i] += t4 * ( 2.0*PI*vr.data[i]*vr.data[i]*(vr.data[i+1]-vr.data[i]) );

				for(int j=i+1; j< num_pr_cellsm1; j++) {
					Alpha.data[i*alpha_stride+j] += t2 * ( 2.0*PI*pow(vr.data[j]/vr.data[i],-LL-1)*vr.data[j]*vr.data[j]*(vr.data[j+1]-vr.data[j-1]) );
					Alpha.data[i*alpha_stride+j] += t4 * ( 2.0*PI*pow(vr.data[j]/vr.data[i],-LL+1)*vr.data[j]*vr.data[j]*(vr.data[j+1]-vr.data[j-1]) ); 
				}
			}
		}
		double ll1 = (double) LL;
		ll1 *= (-0.5)*(ll1 + 1.0);
		for(int i=0; i< num_pr_cells; i++) { // was num_pr_cells1
			//double temp_calc = 1.0 - ll1 * Scattering_Term.data[i];
			//printf("\t diagonal cleanup [%d,%d]=%e  and scattering[i]=%e and additive=%e bytes copied=%d \n", i, i, Alpha.data[i*alpha_stride+i],  Scattering_Term.data[i], temp_calc, (int) alpha_matrix_size__bytes);
			Alpha.data[i*alpha_stride+i] += 1.0 - ll1 * Scattering_Term.data[i];
		}
		
		
		/*
		if(current_x==0 and LL == 1) {
			printf( "Data f1 x=0, l=1\n------------\n");
			double* crap = Scattering_Term.data;
			printf("Scatterign term\n");
			for(int xx=0;xx<num_pr_cells;xx++) {
				double m = *crap;
				printf( "%e ", m);
				crap++;
			}
			
			printf("\nAlpha \n");
			crap = Alpha.data;
				for(int xx=0;xx<3;xx++) {
					for(int yy=0;yy<3;yy++) {
						double m = *crap;
						printf( "%e ", m);
						crap++;
					}
					printf("\n");
				}
			printf("\n \n");
		}
		*/
		// now go and solve this matrix....
		// The tridiagonal version is an truncation, but is much faster
		//		to invert.
		//
		//	The full case is a full, general case inversion. This is the prefered
		//		configuration.
		
		//printf("I am c++ and i say tridiagonal is %d", is_tridiagonal);
		
		if(this->use_C_matrix_solvers) {
			int err;
			if( is_tridiagonal == 1) {
				// the C solver expects the first elemnt of 'a' to be unsused.. no need to set.
				// TODO: if anyone cares.. just put coefficents in these array rather then in full matrix
				//		form (to eliminate the copy)
				//solver_tridiag__a[0] = -666.0;
				solver_tridiag__b[0] = Alpha.data[0];
				solver_tridiag__c[0] = 0.0;
				for(int i=1;i<num_pr_cellsm1;i++) {
					// solver_tridiag__a is the lower diagonal
					solver_tridiag__a[i] = Alpha.data[i*alpha_stride+i-1];
					// solver_tridiag__b is the diagonal
					solver_tridiag__b[i] = Alpha.data[i*alpha_stride+i  ];
					// solver_tridiag__a is the upper diagonal
					solver_tridiag__c[i] = Alpha.data[i*alpha_stride+i+1];
				}
				//solver_tridiag__c[num_pr_cells-1] = -666.0;
				solver_tridiag__b[num_pr_cells-1] = Alpha.data[alpha_stride*alpha_stride - 1];
				solver_tridiag__a[num_pr_cells-1] = Alpha.data[alpha_stride*alpha_stride - 2];
				// the C solver expects the last elemnt of 'c' to be unsused.. no need to set.
				
				err = tridiag__c( solver_tridiag__a,solver_tridiag__b,solver_tridiag__c, fin, solver_tridiag__solution, solver_temp, num_pr_cells);
	#ifndef USE_C_MATRIX_SOLVERS_AND_PYTHON_AND_COMPARE
				memcpy(fin ,solver_tridiag__solution, num_pr_cells*sizeof(double) );
	#else
				// in this case, also do the tridiagonal matrix solve using the Scipy version
				int nothing1 = (*(this->tridiag_matrix_solver))( fin, current_x, LL, Dt );
				// and compare
				for(int i=0;i<num_pr_cells;i++) {
					double diff = fin[i] - solver_tridiag__solution[i];
					double percentError = diff/fin[i];
					if(diff < 0.0) diff *= -1.0; if(percentError<0.0) percentError *= -1.0;					
					if(diff > 1e-8 || percentError > 1e-8) {
						printDebuggingTridiag(fin, i);
					}
				}
	#endif
				if(err < 0) {
					printf("ERROR in C TRIDIAGONAL SOLVER!!!!! QUITTING!!!!! \n\n\n");
					exit(-1);
				}
			} else {
				
				/*
				double indata[9] = {	1.22394929835659072, 0.6153283800788669 , 0.008652750339172455,
									0.5870178580524767 , 1.34035867284704846, 0.14970862337579716 ,
									0.746821615681835  , 0.800571987570456  , 1.7782562541153384};
				double dumbRHS[3] = {0.791091183804244, 1.002017904259492, 1.7104985098039038};

				bool solverError1 = Gauss_Seidel(indata, dumbRHS, solver_temp, solver_scratch, 3,3);
				printf("solution %f   %f    %f\n", solver_temp[0],solver_temp[1],solver_temp[2] );
				*/

				// this write the solution, in place, into fin.. this is exactly what we want anyway.
				//ldu_nonfancy(Alpha.data, num_pr_cells, num_pr_cells, fin, solver_indx, solver_temp);


	#ifndef USE_C_MATRIX_SOLVERS_AND_PYTHON_AND_COMPARE
							// use tri daignoal to guess inital Gauss-Siedel.
							solver_tridiag__b[0] = Alpha.data[0];
							solver_tridiag__c[0] = 0.0;
							for(int i=1;i<num_pr_cellsm1;i++) {
								// solver_tridiag__a is the lower diagonal
								solver_tridiag__a[i] = Alpha.data[i*alpha_stride+i-1];
								// solver_tridiag__b is the diagonal
								solver_tridiag__b[i] = Alpha.data[i*alpha_stride+i  ];
								// solver_tridiag__a is the upper diagonal
								solver_tridiag__c[i] = Alpha.data[i*alpha_stride+i+1];
							}
							//solver_tridiag__c[num_pr_cells-1] = -666.0;
							solver_tridiag__b[num_pr_cells-1] = Alpha.data[alpha_stride*alpha_stride - 1];
							solver_tridiag__a[num_pr_cells-1] = Alpha.data[alpha_stride*alpha_stride - 2];
							// the C solver expects the last elemnt of 'c' to be unsused.. no need to set.
							
							err = tridiag__c( solver_tridiag__a,solver_tridiag__b,solver_tridiag__c, fin, solver_tridiag__solution, solver_temp, num_pr_cells);
							memcpy(solver_temp ,solver_tridiag__solution, num_pr_cells*sizeof(double) );
				// solution is returned in the 'solver_temp' array.
				bool solverError = Gauss_Seidel(Alpha.data, fin, solver_temp, solver_scratch, num_pr_cells,num_pr_cells);
				if(!solverError) {
					printf("Warning! Matrix is not diagonall doinmate so  the implix collision operator may not coverge to a solution.\n");
					exit(-1);
				}
				memcpy(fin ,solver_temp, num_pr_cells*sizeof(double) );			
	#else


							// use tri daignoal to guess inital Gauss-Siedel.
							solver_tridiag__b[0] = Alpha.data[0];
							solver_tridiag__c[0] = 0.0;
							for(int i=1;i<num_pr_cellsm1;i++) {
								// solver_tridiag__a is the lower diagonal
								solver_tridiag__a[i] = Alpha.data[i*alpha_stride+i-1];
								// solver_tridiag__b is the diagonal
								solver_tridiag__b[i] = Alpha.data[i*alpha_stride+i  ];
								// solver_tridiag__a is the upper diagonal
								solver_tridiag__c[i] = Alpha.data[i*alpha_stride+i+1];
							}
							//solver_tridiag__c[num_pr_cells-1] = -666.0;
							solver_tridiag__b[num_pr_cells-1] = Alpha.data[alpha_stride*alpha_stride - 1];
							solver_tridiag__a[num_pr_cells-1] = Alpha.data[alpha_stride*alpha_stride - 2];
							// the C solver expects the last elemnt of 'c' to be unsused.. no need to set.
							
							err = tridiag__c( solver_tridiag__a,solver_tridiag__b,solver_tridiag__c, fin, solver_tridiag__solution, solver_temp, num_pr_cells);
							memcpy(solver_temp ,solver_tridiag__solution, num_pr_cells*sizeof(double) );
				//double* temp = new double[num_pr_cells];
				//double* tempMatrix = new double[num_pr_cells*num_pr_cells];
				//memcpy(temp, fin, num_pr_cells*sizeof(double) );
				//memcpy(tempMatrix, Alpha.data, num_pr_cells*num_pr_cells*sizeof(double) );
				// run the C version 
				bool solverError = Gauss_Seidel(Alpha.data, fin, solver_temp, solver_scratch, num_pr_cells,num_pr_cells);
				if(!solverError) {
					printf("Warning! Matrix is not diagonall doinmate so  the implix collision operator may not coverge to a solution.\n");
					exit(-1);
				}

/*
				// in this case copy the right hand side...
				double* temp = new double[num_pr_cells];
				double* tempMatrix = new double[num_pr_cells*num_pr_cells];
				memcpy(temp, fin, num_pr_cells*sizeof(double) );
				memcpy(tempMatrix, Alpha.data, num_pr_cells*num_pr_cells*sizeof(double) );
				// run the C version 
				int solverError = ldu_nonfancy(tempMatrix, num_pr_cells, num_pr_cells, temp, solver_indx, solver_temp);
				if(solverError!=0) {
					printf("Adam your shitty linalg matrix sover sucks and died %d\n", solverError);
					exit(-1);
				} */
				// run the scipy python version
				int nothing2 = (*(this->full_matrix_solver))( fin, current_x, LL, Dt );			
				// and compare
				double maxElement = 0.0;
				int diaIdx = 0;
				for(int i=0;i<num_pr_cells;i++) {
					if( Alpha.data[diaIdx] > fabs(maxElement) ) maxElement =  Alpha.data[diaIdx];
					diaIdx += (num_pr_cells+1);
				}

				for(int i=0;i<num_pr_cells;i++) {
					//double diff = fin[i] - solver_scratch[i];
					double diff = fin[i] - solver_temp[i];
					double percentError = diff/diaIdx;
					if(diff < 0.0) diff *= -1.0; if(percentError<0.0) percentError *= -1.0;	
					if(diff > 1e-8 || percentError > 1e-8) {
						printf("The C full-matrix  solver got different results from the Scipy (python) one. Aborting (%e,%e, difference %e, percent %e).\n", fin[i], solver_temp[i],diff, percentError);
						printf("The C full-matrix solver got different results from the Scipy (python) one. Aborting. was at spatial cell %d, index %d \n", current_x, i);
						(*(this->debug_matrix_solver))(1);
						exit(-1);
					}
				}
//				delete temp;
//				delete tempMatrix;
	#endif
			}
			
		} else {
			if( is_tridiagonal == 1) {  
				int nothing1 = (*(this->tridiag_matrix_solver))( fin, current_x, LL, Dt );			//this->tridiag_matrix_solver( fin );
			} else {
				int nothing2 = (*(this->full_matrix_solver))( fin, current_x, LL, Dt ); 				//this->full_matrix_solver( fin );
			}
		}
	}
	
	void printDebuggingTridiag(double *fin, int i) {
		double diff = fin[i] - solver_tridiag__solution[i];
		double percentError = diff/fin[i];
		if(diff < 0.0) diff *= -1.0; if(percentError<0.0) percentError *= -1.0;		
		printf("The C tridiagnoal solver got different results from the Scipy (python) one. Aborting (%e,%e, difference %e, percent %e).\n", fin[i], solver_tridiag__solution[i],diff, percentError);
		printf("The python version:\n");
		printf("----------------------------------\n");
		(*(this->debug_matrix_solver))(1);
		printf("\n");
		printf("The C version:\n");
		printf("----------------------------------\n");
		printf("\n\nbelow\n:");
		for(int j=0;j<num_pr_cells;j++) {
			printf("%e,",solver_tridiag__a[j]);
		}
		printf("\n\ndiag:\n");
		for(int j=0;j<num_pr_cells;j++) {
			printf("%e,",solver_tridiag__b[j]);
		}
		printf("\n\nabove:\n");
		for(int j=0;j<num_pr_cells;j++) {
			printf("%e,",solver_tridiag__c[j]);
		}
		
		exit(-1);
	}

		// The Dist Functions's data buffer ihas 3 indicies (in 1D)
	//		and has the shape: [ num_total_cells[0], num_l, num_p]
	//		data.shape[0] = num_total_cells[0], data.shape[1]=shape, data.shape=num_p
	//
	//	so indexing a data element data[k, j, x] looks like:
	//			data_element = data[k*num_l*num_p + j*num_p + x];
	void f1_loop( double dt, int num_m ) {
		const int num_cells__x = F.num_plain_cells[0];
		const int F_data__cell_stride = F.num_l * num_pr_cells;
		const int data_size = num_pr_cells; 
		const int data_size__bytes = num_pr_cells*sizeof(double);
		//printf("IN f1_loop stride: %d\n", F_data__cell_stride);
		
		
		//int starting_F_data_location = F.num_boundry_cells[0][0] * F_data__cell_stride;
		for(int x=0; x < num_cells__x; x++) {		
			// TODO: use higher abstracted addressiing
			int starting_F_data_location = (F.num_boundry_cells[0][0] + x) * F_data__cell_stride;
			double * momentum_projection = &( F_data[starting_F_data_location] );
			//double temp_1 = momentum_projection[1];
			reset_coeff__cpp(  momentum_projection, dt );
			//double temp_12 = momentum_projection[1];
			
			
			starting_F_data_location += num_pr_cells;
			momentum_projection =  &( F_data[starting_F_data_location] );
			//double temp_2 = momentum_projection[1];
			implicit_advance( momentum_projection, 1, -1.0, dt, x);
			//double temp_3 = momentum_projection[1];
			
			//printf("For annoying 'x=%d+%d,L=1' %e, %e, %e, %e \n", x, F.num_boundry_cells[0][0], temp_1, temp_12, temp_2, temp_3);
			
			//starting_F_data_location += F_data__cell_stride;
		}
print_output1 = false;
	}
	
	void flm_loop( double dt, const int num_l, const int num_m ) {
		//printf("IN flm_loop\n");
		if(debug__current_sim_timestep < 10) {
		}
		
		const int num_cells__x = F.num_plain_cells[0];
		//const int F_data__cell_stride = F.num_l * num_pr_cells;
		int starting_F_data_location = F.num_boundry_cells[0][0] * num_l * num_pr_cells;
		
		for(int x=0; x < num_cells__x; x++) {
			
			
			// move pointer to the start of the spatial cell's data
			double * momentum_projection = &( F_data[starting_F_data_location] );
			reset_coeff__cpp(  momentum_projection, dt );
			
			// skip over L=0,L=1 to get to L=2
			// TODO: This is danergous and sloppy because it assume that num_m = 1
			// TODO: fix this.
			//starting_F_data_location += num_pr_cells;
			starting_F_data_location += num_pr_cells;
			printf("num_lnum_l %i \n", num_l);
			for(int L=1; L < num_l; L++) {
			
				// DEBUGING.....
//if(x==1 && L == 2) print_output1 = true;
//else print_output1 = false;
				
				momentum_projection = &( F_data[starting_F_data_location] );
				//double temp_1 = momentum_projection[0];
				implicit_advance(momentum_projection, L, -1.0, dt, x );
				//double temp_2 = momentum_projection[0];
				starting_F_data_location += num_pr_cells;
				//printf("In c++ x=%d, l=%d   %e,%e\n", x, L, temp_1, temp_2);
				/*
				if max_m < l:
					for m in xrange(0, max_m):
						self.implicit_advance( Yin.get_harmonic_at_location(x + self.NB, l, m).momentum_projection, l)					
				else:
					for m in xrange(0, l):
						self.implicit_advance( Yin.get_harmonic_at_location(x + self.NB, l, m).momentum_projection, l)
				*/
			}
		}
	}
	
	void calc_flm( int n_step, double time, double dt ) {
#ifdef USE_C_MATRIX_SOLVERS_AND_PYTHON_AND_COMPARE
		printf("-------------NOTE: Code is compiled with the USE_C_MATRIX_SOLVERS_AND_PYTHON_AND_COMPARE flag. This is significantly slower and is only for validation testing of numerical matrix solving routines.\n");
#endif

		// slot zero is the maximum itteration used this time step
		solverStatsGlobal[0] = 0.0;

		const int num_cells__x = F.num_plain_cells[0];
		const int num_p_cells = num_pr_cells;
		const int num_Ls = F.num_l;
		
		//printf("num_Ls %i \n", num_Ls);
		//printf("num_p_cells %i \n", num_p_cells);
		//printf("num_cells__x %i \n", num_cells__x);
		//printf("F.num_boundry_cells[0][0] %i \n", F.num_boundry_cells[0][0]);
		int starting_F_data_location = F.num_boundry_cells[0][0] * num_Ls * num_p_cells;

		//printf("starting_F_data_location %i \n", starting_F_data_location);

		for(int x=0; x < num_cells__x; x++) {
			double * momentum_projection = &( F_data[starting_F_data_location] );
			reset_coeff__cpp(  momentum_projection, dt );
			starting_F_data_location += num_p_cells;
			
			for(int L=1; L < num_Ls; L++) {
				momentum_projection = &( F_data[starting_F_data_location] );
				implicit_advance(momentum_projection, L, -1.0, dt, x );
				starting_F_data_location += num_p_cells;
			}
		}

		// keep track og the maximum itterations ever across all timesteps
		if(solverStatsGlobal[1] < solverStatsGlobal[0]) solverStatsGlobal[1] = solverStatsGlobal[0];
	}

	// does the same thing as the 'calc_flm' function above, but only on L=1 harmonic
	//		used in the implicit E-Field solver
	void advance_1( int n_step, double time, double dt ) {
	#ifdef USE_C_MATRIX_SOLVERS_AND_PYTHON_AND_COMPARE
			printf("-------------NOTE: Code is compiled with the USE_C_MATRIX_SOLVERS_AND_PYTHON_AND_COMPARE flag. This is significantly slower and is only for validation testing of numerical matrix solving routines.\n");
	#endif

			const int num_cells__x = F.num_plain_cells[0];
			const int num_p_cells = num_pr_cells;
			const int num_Ls = F.num_l;
			
			//printf("num_Ls %i \n", num_Ls);
			//printf("num_p_cells %i \n", num_p_cells);
			//printf("num_cells__x %i \n", num_cells__x);
			//printf("F.num_boundry_cells[0][0] %i \n", F.num_boundry_cells[0][0]);
			int spatialCellDataStartingLocation = F.num_boundry_cells[0][0] * num_Ls * num_p_cells; 
			int starting_F_data_location;// = F.num_boundry_cells[0][0] * num_Ls * num_p_cells;

			//printf("starting_F_data_location %i \n", starting_F_data_location);

			for(int x=0; x < num_cells__x; x++) {

				double * momentum_projection = &( F_data[spatialCellDataStartingLocation] );
				reset_coeff__cpp(  momentum_projection, dt );
				starting_F_data_location = spatialCellDataStartingLocation + num_p_cells;
				
				momentum_projection = &( F_data[starting_F_data_location] );
				implicit_advance(momentum_projection, 1, -1.0, dt, x );

				spatialCellDataStartingLocation +=  num_Ls * num_p_cells;
			}
		}

	// does the same thing as the 'calc_flm' function above, but only on L>1 harmonics
	//		used in the implicit E-Field solver
	void advance_flm( int n_step, double time, double dt ) {
#ifdef USE_C_MATRIX_SOLVERS_AND_PYTHON_AND_COMPARE
		printf("-------------NOTE: Code is compiled with the USE_C_MATRIX_SOLVERS_AND_PYTHON_AND_COMPARE flag. This is significantly slower and is only for validation testing of numerical matrix solving routines.\n");
#endif

		const int num_cells__x = F.num_plain_cells[0];
		const int num_p_cells = num_pr_cells;
		const int num_Ls = F.num_l;
		
		//printf("num_Ls %i \n", num_Ls);
		//printf("num_p_cells %i \n", num_p_cells);
		//printf("num_cells__x %i \n", num_cells__x);
		//printf("F.num_boundry_cells[0][0] %i \n", F.num_boundry_cells[0][0]);
		int starting_F_data_location = F.num_boundry_cells[0][0] * num_Ls * num_p_cells;

		//printf("starting_F_data_location %i \n", starting_F_data_location);

		for(int x=0; x < num_cells__x; x++) {
			double * momentum_projection = &( F_data[starting_F_data_location] );
			reset_coeff__cpp(  momentum_projection, dt );
			// move from L=0 to L=1
			starting_F_data_location += num_p_cells;
			// move from L=1 to L=2 (i.e. we are skipping L=1)
			starting_F_data_location += num_p_cells;

			for(int L=2; L < num_Ls; L++) {
				momentum_projection = &( F_data[starting_F_data_location] );
				implicit_advance(momentum_projection, L, -1.0, dt, x );
				starting_F_data_location += num_p_cells;
			}
		}
	}
		

	
	// Main entry function for higher harmonic collisons calculations.... results will be
	//		written directly into the input distribution function.
	void calc_flm_older( int n_step, double time, double dt ) {
		// this does the integration for (l=1,m=0) and (l=1,m=1)
		//euler_backward_solver_for_fp_advance_1(  Yin, dt);
		//printf("IN CALC_FLM\n");
		
print_output1 = false;
		
		#define num_m 1
		
		if(if_implicit1D) {
			
			// if 1D, only m=0 modes are relevent so just process (L=1,m=0)
			f1_loop(dt, num_m);
			// if 1D, only m=0 modes are relevent so process modes L>=2,m=0
			flm_loop(dt, F.num_l, 1);
		} else {
			// do calculation for L=1 modes, both m values [ (L=1,m=0) and (L=1,m=1) ]
			f1_loop(dt, num_m);
			// do the modes with L >= 2
			flm_loop(dt,  F.num_l, num_m);
		}
	}
		#undef num_m
};
/*
---------------------------------------------------
*/
struct HarmonicCentricView {
	double* data;
	axis *momentum_axis;
	int num_total_cells[3];

	HarmonicCentricView() {
		num_total_cells[0] = 0; num_total_cells[1] = 0; num_total_cells[2] = 0;
	}
};

struct effect_of_E_on_distribution {

	int num_l;
	int num_m;


	numpy_array pr;
	
	numpy_array A1;
	numpy_array A2;
	numpy_array B1;
	numpy_array B2;
	numpy_array C1;
	numpy_array C2;
	numpy_array C3;
	numpy_array C4;
	numpy_array Hp0;

	HarmonicCentricView G;
	HarmonicCentricView H;

	axis momentum_axis;
	axis invpr;

	void init(	species* species, int num_l, int num_m, double *data_A1, double *data_A2, double *data_B1, double *data_B2, 
				double *data_C1, double *data_C2, double* data_C3, double *data_C4,
				double *data_Hp0, double *data_G, double *data_H, double* invpr_values) {
		
		
		this->momentum_axis.init( species->momentum_axis.values,  species->momentum_axis.num_points);
		this->invpr.init( invpr_values,  species->momentum_axis.num_points);
		this->num_l = species->num_l;
		this->num_m = species->num_m;

		A1.data = data_A1; A1.shape[0] = num_l; A1.shape[1] = num_m; A1.dim = 2;
		A2.data = data_A2; A2.shape[0] = num_l; A2.shape[1] = num_m; A2.dim = 2;
		B1.data = data_B1; B1.shape[0] = num_l; B1.dim = 1;
		B2.data = data_B2; B2.shape[0] = num_l; B2.dim = 1;
		C1.data = data_C1; C1.shape[0] = num_l; C1.dim = 1;
		C2.data = data_C2; C2.shape[0] = num_l; C2.shape[1] = num_m; C2.dim = 2;
		C3.data = data_C3; C3.shape[0] = num_l; C3.dim = 1;
		C4.data = data_C4; C4.shape[0] = num_l; C4.shape[1] = num_m; C4.dim = 2;

		Hp0.data = data_Hp0; Hp0.shape[0] = num_l; Hp0.dim = 1;


		A1.data = data_A1; A1.shape[0] = num_l; A1.shape[1] = num_m; A1.dim = 2;
		A1.data = data_A1; A1.shape[0] = num_l; A1.shape[1] = num_m; A1.dim = 2;
		A1.data = data_A1; A1.shape[0] = num_l; A1.shape[1] = num_m; A1.dim = 2;
		A1.data = data_A1; A1.shape[0] = num_l; A1.shape[1] = num_m; A1.dim = 2;


		G.data = data_G;
		H.data = data_H;

	}
};



struct oshun_calcs {
	fokker_plank_explicit_advance _fokker_plank_explicit_advance;
	fokker_plank_implicit_advance _fokker_plank_implicit_advance;
	effect_of_E_on_distribution	  _effect_of_E_on_distribution;
};

oshun_calcs OSHUN_CALCS;


extern "C" {
	extern	NUMPY1_API void general_setup(	int num_species, double* E_field_data, long long* E_field_shape, int E_field_dim, int* E_field_boundry_cells, 
											double* E_field_spatial_axis_values, int E_field_spatial_axis_num_points) {
	
	//extern	NUMPY1_API void general_setup(	int num_species, double* E_field_data, int* E_field_shape, int E_field_dim, int* E_field_boundry_cells, 
	//										double* E_field_spatial_axis_values, int E_field_spatial_axis_num_points) {
		Species.resize( num_species );

		E_field.data = E_field_data;
		E_field.dim = E_field_dim;
		for(int i=0;i < E_field_dim;i++) {
			E_field.shape[i] = (int) E_field_shape[i];		// NEW for 64 bit windows
		}
		E_field.boundry_cells[0][0] = E_field_boundry_cells[0]; E_field.boundry_cells[0][1] = E_field_boundry_cells[1];
		E_field.boundry_cells[1][0] = E_field_boundry_cells[2]; E_field.boundry_cells[1][1] = E_field_boundry_cells[3];
		E_field.boundry_cells[2][0] = E_field_boundry_cells[4]; E_field.boundry_cells[2][1] = E_field_boundry_cells[5];
		
		SpatialAxes.resize( E_field_dim );
		for(int i=0;i < E_field_dim;i++) {
			SpatialAxes[i].values = E_field_spatial_axis_values;
			SpatialAxes[i].num_points = E_field_spatial_axis_num_points;
			SpatialAxes[i].dx = E_field_spatial_axis_values[1] - E_field_spatial_axis_values[0];
			SpatialAxes[i].dp = E_field_spatial_axis_values[1] - E_field_spatial_axis_values[0];

			E_field.spatial_axes[0].values = E_field_spatial_axis_values;
			E_field.spatial_axes[0].num_points = E_field_spatial_axis_num_points;
			E_field.spatial_axes[0].dx = E_field_spatial_axis_values[1] - E_field_spatial_axis_values[0];
			E_field.spatial_axes[0].dp = E_field_spatial_axis_values[1] - E_field_spatial_axis_values[0];
		}
	}

	extern	NUMPY1_API void general_species_setup(	int species_idx, int num_l, int num_m, 
													double* mometum_axis_values, int momentum_axis_num_points,
													int* F_num_total_cells, int* num_plain_cells, int* num_boundry_cells,
													double* F_data) {
		Species[species_idx].num_l = num_l;
		Species[species_idx].num_l = num_l;
		Species[species_idx].init(	num_l, num_m, mometum_axis_values, momentum_axis_num_points, 
									F_num_total_cells, num_plain_cells, num_boundry_cells, F_data );
		
		// TODO enable mutli-species.
		OSHUN_CALCS._fokker_plank_implicit_advance.explict_code = &(OSHUN_CALCS._fokker_plank_explicit_advance);

	}

	extern NUMPY1_API void transfer_species_F_buffers( int species_idx, int harmonic_indx, int ix, double* data) {
		SphericalHarmonic * sh = &(Species[species_idx].F.cells[ix].harmonics[harmonic_indx]);
		sh->momentum_projection = data;
	}
/*
---------------------------- Effect of E ---------------------------
*/
	extern NUMPY1_API void init_effect_of_E_on_distribution(	int species_index, int num_l, int num_m, 
																double *data_A1, double *data_A2, 
																double *data_B1, double *data_B2, 
																double *data_C1, double *data_C2, 
																double* data_C3, double *data_C4,
																double *data_Hp0, double *data_G, 
																double *data_H, double* invpr_values) {
		DEBUG("\t Entering init_effect_of_E_on_distribution \n");
		OSHUN_CALCS._effect_of_E_on_distribution.init(	&(Species[species_index]), num_l, num_m,
														data_A1, data_A2,
														data_B1, data_B2,
														data_C1, data_C2,
														data_C3, data_C4,
														data_Hp0, data_G,
														data_H, invpr_values);
		DEBUG("\t  Leaving init_effect_of_E_on_distribution \n");

	}
	
	extern NUMPY1_API void effect_of_E_on_distribution_calc( ) {
		
	}  
		
/* 
------------------------------- Collsions------------------------------
*/

	
	extern  NUMPY1_API void init__fokker_plank__implicit(	double *data_vr, long long* shape_vr,
															double *data_vr3,
															double *data_df0,
															double *data_ddf0,
															double *data_Scattering_Term,
															double *data_Alpha,
															double *data_AlphaTri,
															double *data_J1m,
															double *data_TriI1,
															double *data_TriI2,
															double *data_IvDnDm1,
															double *data_IvDnDp1,
															double *data_Ivsq2Dn,
															double *data_I0,
															double kpre, double zeta, double f00_factor,
															int is_tridiagonal, int if_implicit1D, int species_index,
															double* _LOGee, double* _ZLOGei,
															matrix_solver_t tridiag_solver, matrix_solver_t full_solver, no_arg_function_t debug_matrix_solver,
															double* solverStats) {
		DEBUG("\t Entering init__fokker_plank__implicit \n");
		OSHUN_CALCS._fokker_plank_implicit_advance.init(	data_vr, shape_vr, 
															data_vr3,
															data_df0,
															data_ddf0,
															data_Scattering_Term,
															data_Alpha,
															data_AlphaTri,
															data_J1m,
															data_TriI1,
															data_TriI2,
															data_IvDnDm1,
															data_IvDnDp1,
															data_Ivsq2Dn,
															data_I0,
															kpre, zeta, f00_factor,
															is_tridiagonal, if_implicit1D, species_index,
															_LOGee, _ZLOGei,
															tridiag_solver,full_solver, debug_matrix_solver,
															solverStats);
		DEBUG("\t  Leaving init__fokker_plank__implicit \n");
	}

	extern  NUMPY1_API void fokker_plank__implicit__reset_coeff( double* fin, double dt	) {
		DEBUG("\t Entering fokker_plank__implicit__reset_coeff \n");
		OSHUN_CALCS._fokker_plank_implicit_advance.reset_coeff__cpp( fin, dt);
		DEBUG("\t  Leaving fokker_plank__implicit__reset_coeff \n");
	}
	
	
	extern  NUMPY1_API void fokker_plank__implicit__advance(double * fin, int LL, double _LOGee, double Dt, int current_x_index) {
		DEBUG("\t Entering fokker_plank__implicit__advance \n");
		OSHUN_CALCS._fokker_plank_implicit_advance.implicit_advance(fin, LL, _LOGee, Dt, current_x_index);
		DEBUG("\t  Leaving fokker_plank__implicit__advance \n");
	}
	
	extern  NUMPY1_API void fokker_plank__implicit__calc_flm( int n_step, double time, double dt) {
		DEBUG("\t Entering fokker_plank__implicit__calc_flm \n");
		OSHUN_CALCS._fokker_plank_implicit_advance.calc_flm( n_step, time, dt );
		DEBUG("\t  Leaving fokker_plank__implicit__calc_flm \n");
	}

	// needed for implict E calculations.. does what the calc_clm version deos but only on the l=1 harmonic
	extern  NUMPY1_API void fokker_plank__implicit__advance_1( int n_step, double time, double dt) {
		DEBUG("\t Entering fokker_plank__implicit__calc_flm \n");
		OSHUN_CALCS._fokker_plank_implicit_advance.advance_1( n_step, time, dt );
		DEBUG("\t  Leaving fokker_plank__implicit__calc_flm \n");
	}
	// needed for implict E calculations.. does what the calc_clm version deos but only on the l=1 harmonic
	extern  NUMPY1_API void fokker_plank__implicit__advance_flm( int n_step, double time, double dt) {
		DEBUG("\t Entering fokker_plank__implicit__calc_flm \n");
		OSHUN_CALCS._fokker_plank_implicit_advance.advance_flm( n_step, time, dt );
		DEBUG("\t  Leaving fokker_plank__implicit__calc_flm \n");
	}
	extern  NUMPY1_API void init__fokker_plank__explicit(	double *data_vr, long long* shape_vr, 
															double *data_U4, double *data_U4m1,
															double *data_U2, double *data_U2m1,
															double *data_U1, double *data_U1m1,
															double *data_J1, double *data_I4,
															double *data_U3, double *data_Qn,
															double *data_Pn, double* data_I2, 
															double data_G_constant_1, double* data_G_constant_vr3,
															double *data_G_constant_vr5, double* data_G_constant_vr5_2,
															double *data_G_constant_vr7, double data_density_np,
															int NB, double c_kpre, int num_subcycling_steps, int species_index) {
		DEBUG("\t Entering init__fokker_plank__explicit \n");
		OSHUN_CALCS._fokker_plank_explicit_advance.init(	data_vr, shape_vr, data_U4, data_U4m1,
															data_U2, data_U2m1,data_U1, data_U1m1,
															data_J1, data_I4, data_U3, data_Qn,
															data_Pn, data_I2, data_G_constant_1,
															data_G_constant_vr3, data_G_constant_vr5, 
															data_G_constant_vr5_2, data_G_constant_vr7, 
															data_density_np, NB, c_kpre, num_subcycling_steps, species_index );
		DEBUG("\t  Leaving init__fokker_plank__explicit \n");
	}

	extern  NUMPY1_API void fokker_plank_explicit_slope_calc(double *fin, double *fh) {
		DEBUG("\t Entering fokker_plank_explicit_slope_calc \n");
		OSHUN_CALCS._fokker_plank_explicit_advance.fokker_plank_explicit_slope_calc( fin, fh );
		DEBUG("\t  Leaving fokker_plank_explicit_slope_calc \n");
	}
	
	extern  NUMPY1_API void fokker_plank_explicit__calc_f00( int n_step, double time, double dt) {
		DEBUG("\t Entering fokker_plank_explicit__calc_f00 \n");
		OSHUN_CALCS._fokker_plank_explicit_advance.calc_f00( n_step, time, dt );
		DEBUG("\t  Leaving fokker_plank_explicit__calc_f00 \n");
	}

	
	
	// abject tester funcitona......
	extern  NUMPY1_API void testa_cpp_passing( double* _time) {
		printf( "I see this as A value %f \n", (*_time) );
		(*_time) = 666.0;
	}
	extern  NUMPY1_API void testb_cpp_passing( double* _time2) {
		printf( "I see this as B value %f \n", (*_time2) );
		(*_time2) += 100.0;
	}
	
	extern NUMPY1_API void advance_cpp_sim_clock( int time_step, double time, double dt) {
		//DEBUG("%d , %e, %e \n", time_step , time, dt);
		
		debug__current_sim_timestep = time_step;
		debug__current_sim_time = time;
		debug__sim_current_dt = dt;
	}

	extern NUMPY1_API unsigned int get_pointer_size_in_bytes() {
		return sizeof(size_t);
	}
	extern NUMPY1_API unsigned int get_int_size_in_bytes() {
		return sizeof(int);
	}
	extern NUMPY1_API unsigned int get_float_size_in_bytes() {
		return sizeof(float);
	}
	extern NUMPY1_API unsigned int get_double_size_in_bytes() {
		return sizeof(double);
	}
	extern NUMPY1_API unsigned int get_long_size_in_bytes() {
		return sizeof(long);
	}
	extern NUMPY1_API unsigned int get_longlong_size_in_bytes() {
		return sizeof(long long);
	}

}





extern "C" {

}
