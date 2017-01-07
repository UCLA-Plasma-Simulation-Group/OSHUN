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
#ifdef _WIN32 || _WIN64
	#include "stdafx.h"
#endif

#include "numpy1.h"
#include <stdio.h>
#include <cstring>			// for memcpy
#include <vector>
#include <math.h>

using namespace std;
#define PI 3.14159265358979323846

// This is an example of an exported variable
NUMPY1_API int nnumpy1=0;

// This is an example of an exported function.
NUMPY1_API int fnnumpy1(void)
{
	return 42;
}


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
	int dim_temp2[1];
	int dim_temp3[1];
	int dim_temp4[1];
	int dim_temp5[1];

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
			printf("Ddecutrcution and meyhem!!! \n")
			delete data
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

	// precalculations for the lower |p| cells
	double		G_constant_1;
	numpy_array G_constant_vr3;
	numpy_array G_constant_vr5;
	numpy_array G_constant_vr5_2;
	numpy_array G_constant_vr7;
	double density_np;
	int NB;
	double c_kpre;

	void init(	double *data_vr, int* shape_vr,
				double *data_U4, double *data_U4m1,
				double *data_U2, double *data_U2m1,
				double *data_U1, double *data_U1m1,
				double *data_J1, double *data_I4,
				double *data_U3, double *data_Qn,
				double *data_Pn, double* data_I2, 
				double data_G_constant_1, double* data_G_constant_vr3,
				double *data_G_constant_vr5, double* data_G_constant_vr5_2,
				double *data_G_constant_vr7, double data_density_np, int NB, double c_kpre) {

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
		for (int n=0; n < num_pr_cells; ++n) {
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
	
	double _LOGee;
	double _ZLOGei;
	double _zeta;
	double f00_factor;
	double I0_density;
	
	
	double kpre;
	#define  four_pi 4.0*3.141592654
	#define eight_pi 8.0*3.141592654
	
	
	void init(	double *data_vr, int* shape_vr,
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
				double kpre, double _zeta, double f00_factor) {
		vr.data = data_vr; vr.shape[0] = shape_vr[0]; vr.dim = 1;
		vr3.data = data_vr3; vr3.shape[0] = shape_vr[0]; vr3.dim = 1;
		df0.data = data_df0; df0.shape[0] = shape_vr[0]; df0.dim = 1;
		ddf0.data = data_ddf0; ddf0.shape[0] = shape_vr[0]; ddf0.dim = 1;
		Alpha.data = data_Alpha; Alpha.shape[0] = shape_vr[0]; Alpha.shape[1] = shape_vr[0]; Alpha.dim = 2;
		AlphaTri.data = data_AlphaTri; AlphaTri.shape[0] = shape_vr[0]; AlphaTri.shape[1] = shape_vr[0]; AlphaTri.dim = 2;
		Scattering_Term.data = data_Scattering_Term; Scattering_Term.shape[0] = shape_vr[0]; Scattering_Term.dim; 
		
		J1m.data = data_J1m; J1m.shape[0] = shape_vr[0]; J1m.dim = 1;
		I0.data =  data_I0;   I0.shape[0] = shape_vr[0];  I0.dim = 1;
		
		TriI1.data = data_TriI1; TriI1.shape[0] = shape_vr[0]; TriI1.dim = 1;
		TriI2.data = data_TriI2; TriI2.shape[0] = shape_vr[0]; TriI2.dim = 1;
		
		IvDnDm1.data = data_IvDnDm1; IvDnDm1.shape[0] = shape_vr[0]; IvDnDm1.dim = 1;
		IvDnDp1.data = data_IvDnDp1; IvDnDp1.shape[0] = shape_vr[0]; IvDnDp1.dim = 1;
		Ivsq2Dn.data = data_Ivsq2Dn; Ivsq2Dn.shape[0] = shape_vr[0]; Ivsq2Dn.dim = 1;
		
		num_pr_cells = shape_vr[0];
		this->kpre = kpre;
		this->_zeta = _zeta;
		this->f00_factor = f00_factor;
		this->alpha_matrix_size__bytes = sizeof(double)*Alpha.shape[0];
		
		// TOTOD: REMOVE THIS... ITS A BIT CONVOLUTED.
		this->_LOGee = -1.0;
		this->_ZLOGei = -1.0;

	}
	
	void reset_coeff__cpp( double * fin, double Dt) {
		
		const int num_pr_cellsm1 = num_pr_cells -1; 
		
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
			I2.data[k] += I2.data[k-1];
			I2.data[k] *= four_pi;
		}
		double I2_temperature = I2.data[num_pr_cells - 1];
		for (int k=0; k<num_pr_cells;k++) {
			I2.data[k] /=  (vr.data[k]*vr.data[k]);
		}
		
		// Density integral I0 = 4*pi*int_0^v f(u)*u^2du 
		I0.data[0] = 0.0;
		for(int k=1; k<num_pr_cells;k++) {
			I0.data[k]  = U2.data[k]*fin[k]+U2m1.data[k]*fin[k-1];
			I0.data[k] += I0.data[k-1];
			I0.data[k] *= four_pi;
		}
		I0_density = I0.data[num_pr_cells - 1];
		
		// Integral J_(-1) = 4*pi * v * int_0^v f(u)*u^4du
		J1m.data[ num_pr_cells - 1] = 0.0;
		for(int k=(num_pr_cells-2); k<=0 ;k--) {
		  J1m.data[k]  = U1.data[k+1]*fin[k+1]+U1m1.data[k+1]*fin[k]; 
		  J1m.data[k] += J1m.data[k+1];
		}
		for (int k=0; k<num_pr_cells;k++) {
			J1m.data[k] *= four_pi * vr.data[k];
		}
		
		// COULOMB LOGARITHMS
		// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		_LOGee  = explict_code->LOGee( I0_density, I2_temperature );
		_ZLOGei = explict_code->ZLOGei( I0_density, I2_temperature, _zeta );
		
		// BASIC INTEGRALS FOR THE TRIDIAGONAL PART
		//   and SCATTERING TERM
		// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		
		
		// Take care of 0th element...
		Scattering_Term.data[0] = (_ZLOGei * I0_density); 
		Scattering_Term.data[0] *= (kpre*Dt/vr3.data[0]);
		
		for(int i=1; i<num_pr_cells;i++) {
			double temp =   I0.data[i] + (2.0*J1m.data[i] - I2.data[i]) / 3.0;
			TriI2.data[i] = ( I2.data[i] + J1m.data[i] ) / 3.0;             		// ( I2 + J_{-1} ) / 3
			TriI1.data[i] = temp; 													// (-I2 + 2*J_{-1} + 3*I0) / 3
			Scattering_Term.data[i] = temp*_LOGee;
			Scattering_Term.data[i] += (_ZLOGei*I0_density);
			Scattering_Term.data[i] *= (kpre*Dt/vr3.data[i]);
		}
		
		
		const int alpha_stride = Alpha.shape[1];
		// Todo: is this zerosing necessary?]
		double factor = (-1.0) * _LOGee * kpre * Dt;
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
	
	// TODO: take out the __LOGee.. in cpp mode we calulate in cpp and it is saved as a member variable..
	//		and its not used by the python card in any other way at this point...
	void implicit_advance(double * fin, int LL, double LOGee_from_python, double Dt, int is_tridiagonal) {
		
		// if the member variable '_LOGee' is positive, then the cpp code caluclated this value
		//	in 'reset_coeff__cpp'. If not, the the python code calcualted '_LOGee' so use that value.
		double LOGee;
		if ( LOGee_from_python > 0.0) {
			LOGee = LOGee_from_python;
		} else {
			LOGee = this->_LOGee;
		}
		
		const int num_pr_cellsm1 = num_pr_cells -1; //const int num_pr_cells1 = num_pr_cells;
		const double factor1 = (-1.0) * LOGee * kpre * Dt;
		const int alpha_stride = Alpha.shape[1];
		
		// the Alpha_Tri (tridagonal) array is reused, but the matrix solver is decructive..
		//	so we gotta copy Alpha_Tri in the buffer Alpha (Alpha is allowed to be killed by matrix solver)
		memcpy( Alpha.data, AlphaTri.data, alpha_matrix_size__bytes);
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
			Alpha.data[i*alpha_stride+i] += 1.0 - ll1 * Scattering_Term.data[i];
		}
	}
	
	// main entry point for (L=0, m=0) collsions.. results will be written
	//		directly into the input distribution function.
	void calc_f00( Yin, species, state, CFG, n_step, time ) {
		// loop over the useful (non-guard) cells.. 
		// (Our Fokker-Plank calculation is completely local in space and is strictly
		//		an intra-spatial-cell calculation.. so no guard cell are required )
		for(int x=0; x < 
	}
	
	// Main entry function for higher harmonic collisons calculations.... results will be
	//		written directly into the input distribution function.
	void calc_flm() {
		// this does the integration for (l=1,m=0) and (l=1,m=1)
		euler_backward_solver_for_fp_advance_1(  Yin, dt);
		// this does the remaining higher harmonics.. (i.e. harmonic with L >= 2
		euler_backward_solver_for_fp_advance_lm(Yin, dt);
	}

};


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
	//	so indexing an element data[k, j, x] looks like:
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
	}
};


vector<axis> SpatialAxes;
vector<species> Species;
Field1D E_field;

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
	extern	NUMPY1_API void general_setup(	int num_species, double* E_field_data, int* E_field_shape, int E_field_dim, int* E_field_boundry_cells, 
											double* E_field_spatial_axis_values, int E_field_spatial_axis_num_points) {
		Species.resize( num_species );

		E_field.data = E_field_data;
		E_field.dim = E_field_dim;
		for(int i=0;i < E_field_dim;i++) {
			E_field.shape[i] = E_field_shape[i];
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
		OSHUN_CALCS._effect_of_E_on_distribution.init(	&(Species[species_index]), num_l, num_m,
														data_A1, data_A2,
														data_B1, data_B2,
														data_C1, data_C2,
														data_C3, data_C4,
														data_Hp0, data_G,
														data_H, invpr_values);


	}
	
	extern NUMPY1_API void effect_of_E_on_distribution_calc( ) {
		
	}  
		
/* 
------------------------------- Collsions------------------------------
*/
	extern  NUMPY1_API void init__fokker_plank__implicit(	double *data_vr, int* shape_vr,
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
															double kpre, double _zeta, double f00_factor) {
		
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
															kpre, _zeta, f00_factor);
	}

	extern  NUMPY1_API void fokker_plank__implicit__reset_coeff( double* fin, double dt	) {
		OSHUN_CALCS._fokker_plank_implicit_advance.reset_coeff__cpp( fin, dt);
	}
	
	
	extern  NUMPY1_API void fokker_plank__implicit__advance(double * fin, int LL, double _LOGee, double Dt, int is_tridiagonal) {
		OSHUN_CALCS._fokker_plank_implicit_advance.implicit_advance(fin, LL, _LOGee, Dt, is_tridiagonal);
	}


	extern  NUMPY1_API void init__fokker_plank__explicit(	double *data_vr, int* shape_vr, 
															double *data_U4, double *data_U4m1,
															double *data_U2, double *data_U2m1,
															double *data_U1, double *data_U1m1,
															double *data_J1, double *data_I4,
															double *data_U3, double *data_Qn,
															double *data_Pn, double* data_I2, 
															double data_G_constant_1, double* data_G_constant_vr3,
															double *data_G_constant_vr5, double* data_G_constant_vr5_2,
															double *data_G_constant_vr7, double data_density_np,
															int NB, double c_kpre ) {
		
		OSHUN_CALCS._fokker_plank_explicit_advance.init(	data_vr, shape_vr, data_U4, data_U4m1,
															data_U2, data_U2m1,data_U1, data_U1m1,
															data_J1, data_I4, data_U3, data_Qn,
															data_Pn, data_I2, data_G_constant_1,
															data_G_constant_vr3, data_G_constant_vr5, 
															data_G_constant_vr5_2, data_G_constant_vr7, 
															data_density_np, NB, c_kpre );
	}

	extern  NUMPY1_API void fokker_plank_explicit_slope_calc(double *fin, double *fh) {
		OSHUN_CALCS._fokker_plank_explicit_advance.fokker_plank_explicit_slope_calc( fin, fh );
	}
}





extern "C" {

}
