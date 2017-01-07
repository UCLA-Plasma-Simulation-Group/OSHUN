// kernals
/*
	Version 4:
		Cleanup up code... made it so that there were no register spillings or stack frmae bytes genrated as reported by PTAX.. This is the baseline version
		for the furthur cleanups.

	Version 5:
		Try to minimize all mem access... calculated most everything in the from scratch in the kernal.. Only question is how well looks work for the
		p bin values... should they be repeated for each thread? should there be 1 copy? will it be cached?
*/
#include <iostream>
#include <vector>
#include <valarray>
#include <complex>
//#include <stdint.h>

//#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <float.h>
#include <cuda.h>
#include <cuda_runtime_api.h>
//#include <thrust/scan.h>
#include "matrices.h"

// Interface 
//#include "decl-input.h"
//#include "decl-state.h"
//#include "decl-fokkerplanck.h"

//#include "oshun_cuda.h"
#include "oshun_cuda.cuh"
#define FP_TYPE double
#define DATA double *
#define uint unsigned int
#define M_PI 3.14159265358979323846


extern "C" void set_FP_constants(int _leng, FP_TYPE _c_kpre, uint _NB, FP_TYPE _density_np);
extern "C" void simple_test_kernal();
extern "C" void add_1_to_all();
extern "C" void cudasafe( cudaError error, const char *message);
extern "C" void scan_test();
extern "C" void scan_test_4_at_a_time();
extern "C" void eval_test(	double * fc, double * dest, double *vr, 
							double * U1, double *U1m1, double * U2, double *U2m1,
							double * U3,
							double * U4, double *U4m1,
							double * Pn, double *Qn,
							uint size, int debug);

extern "C" void test_scans(double * crap);

extern "C" void eval_rk4_v6 (	double * fc, double * dest, double *vr, 
							//double * U1, double *U1m1, double * U2, double *U2m1,
							//double * U3,
							//double * U4, double *U4m1,
							//double * Pn, double *Qn,
							double * __restrict__ cell_data__precomp1,   // these should be the size of the number of cells.*pr
							double * __restrict__ cell_data__precomp2,
							double * __restrict__ cell_data__precomp3,
							double * __restrict__ cell_data__precomp4,
							double * __restrict__ cell_data__precomp6,
							uint size, int padded_size, int debug, double h, int numh, int num_cells_x, int num_cells_y,
							const double * __restrict__ U4,
									const double * __restrict__ U4m1,
									const double * __restrict__ U2,
									const double * __restrict__ U2m1);
								
//extern GpuProperties gpu_properties;

#define id threadIdx.x
#define bi blockIdx.x == 0 &&

#define cutilCheckMsg(msg)           __cutilGetLastError (msg, __FILE__, __LINE__)

inline void __cutilGetLastError( const char *errorMessage, const char *file, const int line )
{
    cudaError_t err = cudaGetLastError();
    if( cudaSuccess != err) {
        printf("%s(%i) : cutilCheckMsg() CUTIL CUDA error : %s : (%d) %s.\n", file, line, errorMessage, (int)err, 
			cudaGetErrorString( err ) );
        exit(-1);
    }
} 

// constants for FP evaluatoin...
__constant__ int leng;
__constant__ FP_TYPE c_kpre;
__constant__ uint NB;
__constant__ double density_np;
__device__ double LOG_LAMBDA;

void set_FP_constants(int _leng, FP_TYPE _c_kpre, uint _NB, FP_TYPE _density_np) {

	cudaMemcpyToSymbol("leng", &_leng, 4, 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol("c_kpre", &_c_kpre, 8, 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol("NB", &_NB, 4, 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol("density_np", &_density_np, 8, 0, cudaMemcpyHostToDevice);
}



inline __device__ double LOGee_cuda__v2(double ne, double Te) {
		
        double lnee;
        
//      if the density is positive
        if (ne > 0.000000001) {
            Te /= (3.0*ne);
            Te *= 511000; // Temperature in eV
            ne *= density_np;

            Te = (double) __logf((float) Te); 
            ne = (double) __logf((float) ne);
            lnee = 23.5 - 0.5*ne + 1.25*Te - sqrt(0.00001+0.0625*(Te-2.0)*(Te-2.0));

            if (lnee > 2.0) return lnee;
        }
        // Default minimum "2"

        return 2.0; 
}

inline __device__ double LOGee_cuda(double ne, double Te) {
//-------------------------------------------------------------------
//   Calculate the Coulomb logarithm for electron-electron collisions
//-------------------------------------------------------------------
//      Note: that the results here assume the distribution functions
//      are nonrelativistic, as is the case for the rest of the F-P 
//      part of the code.

        double lnee;
        
//      if the density is positive
        if (ne > 0.000000001) {
            Te /= (3.0*ne);
            Te *= 511000; // Temperature in eV
            ne *= density_np;

            Te = log(Te); 
            ne = log(ne);
            lnee = 23.5 - 0.5*ne + 1.25*Te - sqrt(0.00001+0.0625*(Te-2.0)*(Te-2.0));

            if (lnee > 2.0) return lnee;
        }
        // Default minimum "2"

        return 2.0; 
}


// TODO: this can be loads better
__device__  double G_cuda(const int n, const DATA fin,
					const DATA vr, const FP_TYPE J1) {
	double i2s, i4s;
	double f00( (fin[0] - fin[1]*(vr[0]*vr[0])/(vr[1]*vr[1]))/ (1.0 - (vr[0]*vr[0])/(vr[1]*vr[1])) );
	//printf("\tf00 %e\n", f00);
	i2s = f00*pow(vr[n],3)/3.0 + (fin[1]-f00)*pow(vr[n],5)/(vr[1]*vr[1])*0.2;
	//printf("\ti2s %e\n", i2s);
	i4s = f00*pow(vr[n],5)*0.2 + (fin[1]-f00)*pow(vr[n],7)/(vr[1]*vr[1]*7.0);
	//printf("\ti4s %e\n", i4s);
	return fin[n]*i4s + (pow(vr[n],3)*fin[n]-3.0*i2s) * J1;
}


////////////////////////////////////////////////////////////////////////////////
// Scan kernels double
////////////////////////////////////////////////////////////////////////////////
// orginall 256
#define THREADBLOCK_SIZE 512
#define LOG2_WARP_SIZE 5U
#define WARP_SIZE (1U << LOG2_WARP_SIZE)

// REVERSE SCAN KERNALS
//-------------------------------------------------------------------------------------------------
	//Almost the same as naive scan1Inclusive, but doesn't need __syncthreads()
	//assuming size <= WARP_SIZE
	inline __device__ double warpScanInclusive_rev(double idata, volatile double *s_Data, uint size) {
		// this genrates a a sequence of [0,size][2*(size+1),(3*size-1)]
		
		//uint pos = 2 * threadIdx.x - (threadIdx.x & (size - 1));
		
		// diff
		uint pos = (2*threadIdx.x - 3*(threadIdx.x & (size - 1))) + size - 1;
		s_Data[pos] = 0.0;
		
		//pos += (2*(size-threadIdx.x) -1); // DUMB TODO: ADAM FACTOR THIS BETTER BUT I AM TIRED //size;
		//pos = 2*size - 1 - pos;
		pos += size;
		s_Data[pos] = idata;
		
		//printf("Thread %d (size %d) mapped to %d (%d)\n", threadIdx.x, size, pos, pos-size);
		
		for(uint offset = 1; offset < size; offset <<= 1) {
			s_Data[pos] += s_Data[pos - offset];
		}
		
		return s_Data[pos];
	}

	inline __device__ double warpScanExclusive_rev(double idata, volatile double *s_Data, uint size) {
		return warpScanInclusive_rev(idata, s_Data, size) - idata;
	}

	inline __device__ double scan1Inclusive_rev(double idata, volatile double *s_Data, uint size) {
		
		if(size > WARP_SIZE) {
			//Bottom-level inclusive warp scan
			double warpResult = warpScanInclusive_rev(idata, s_Data, WARP_SIZE);
			
			//Save top elements of each warp for exclusive warp scan
			//sync to wait for warp scans to complete (because s_Data is being overwritten)
			__syncthreads();
			
			//if( (threadIdx.x & (WARP_SIZE - 1)) == (WARP_SIZE - 1) ) {
			if( ((threadIdx.x-1) & (WARP_SIZE - 1)) == (WARP_SIZE - 1) )  // so ghetto
				s_Data[threadIdx.x >> LOG2_WARP_SIZE] = warpResult;
			
			// wait for warp scans to complete
			__syncthreads();
			if( threadIdx.x < (THREADBLOCK_SIZE / WARP_SIZE) ) {
				//grab top warp elements
				double val = s_Data[threadIdx.x];
				//calculate exclsive scan and write back to shared memory
				s_Data[threadIdx.x] = warpScanExclusive_rev(val, s_Data, size >> LOG2_WARP_SIZE);
			}

			//return updated warp scans with exclusive scan results
			__syncthreads();
			return warpResult + s_Data[threadIdx.x >> LOG2_WARP_SIZE];
		} else { 
			return warpScanInclusive_rev(idata, s_Data, size);
		}
	}

	inline __device__ double scan1Exclusive_rev(double idata, volatile double *s_Data, uint size){
		return scan1Inclusive_rev(idata, s_Data, size) - idata;
	}
// Simple Forward Scan Kernal
//-------------------------------------------------------------------------------------------------
	// simple 128
	inline __device__ double warpScanInclusiveSimple(double idata, double *s_Data, uint size) {
		
		s_Data[threadIdx.x] = idata;
		
		if(threadIdx.x > 0) s_Data[threadIdx.x] += s_Data[threadIdx.x - 1];
		__syncthreads();
		if(threadIdx.x > 1) s_Data[threadIdx.x] += s_Data[threadIdx.x - 2];
		__syncthreads();
		if(threadIdx.x > 3) s_Data[threadIdx.x] += s_Data[threadIdx.x - 4];
		__syncthreads();
		if(threadIdx.x > 7) s_Data[threadIdx.x] += s_Data[threadIdx.x - 8];
		__syncthreads();
		if(threadIdx.x > 15) s_Data[threadIdx.x] += s_Data[threadIdx.x - 16];
		__syncthreads();
		if(threadIdx.x > 31) s_Data[threadIdx.x] += s_Data[threadIdx.x - 32];
		__syncthreads();
		if(threadIdx.x > 63) s_Data[threadIdx.x] += s_Data[threadIdx.x - 64];
		__syncthreads();
		
		//if( blockIdx.x == 0 ) {
		//	printf("|%d:%d: %e| \n", blockIdx.x, threadIdx.x, s_Data[threadIdx.x]);
		//}
		return s_Data[threadIdx.x];
	}
	
	// assume 128 for now... and WARP_SIZE=32
	inline __device__ double warpScanInclusiveSimpleWarp(double idata, volatile double *s_Data , uint size) {   //, volatile double *s_data__warp_results, uint size) {
		
		s_Data[threadIdx.x] = idata;
		__syncthreads(); // suspect
		
		//uint warp_index = (threadIdx.x & (size - 1));
		uint warp_index = (threadIdx.x & 31);
		
		if(warp_index > 0) s_Data[threadIdx.x] += s_Data[threadIdx.x - 1];
		//__syncthreads();
		if(warp_index > 1) s_Data[threadIdx.x] += s_Data[threadIdx.x - 2];
		//__syncthreads();
		if(warp_index > 3) s_Data[threadIdx.x] += s_Data[threadIdx.x - 4];
		//__syncthreads();
		if(warp_index > 7) s_Data[threadIdx.x] += s_Data[threadIdx.x - 8];
		//__syncthreads();
		if(warp_index > 15) s_Data[threadIdx.x] += s_Data[threadIdx.x - 16];
		__syncthreads();
		
		double intra_warp_result = s_Data[threadIdx.x];
		
		//return intra_warp_result;
		
		if(warp_index == 31) {
			s_Data[threadIdx.x >> LOG2_WARP_SIZE] = intra_warp_result;
		}
		__syncthreads();
		
		if(threadIdx.x < 5) {
			double my_val = s_Data[threadIdx.x];
			
			if(threadIdx.x > 0) s_Data[threadIdx.x] += s_Data[threadIdx.x - 1];
			if(threadIdx.x > 1) s_Data[threadIdx.x] += s_Data[threadIdx.x - 2];
			s_Data[threadIdx.x] -= my_val;
		}
		__syncthreads();
		
		
		return s_Data[threadIdx.x >> LOG2_WARP_SIZE] + intra_warp_result; 
	}

		// assume 128 for now... and WARP_SIZE=32
	inline __device__ double warpScanInclusiveSimpleWarp_rev(double idata, volatile double *s_Data , uint size) {   //, volatile double *s_data__warp_results, uint size) {
		
		__syncthreads();
		s_Data[threadIdx.x] = idata;
		__syncthreads(); // suspect
		
		
		
		//uint warp_index = (threadIdx.x & (size - 1));
		uint warp_index = (threadIdx.x & 31);
		__syncthreads();
		if(warp_index < 31) s_Data[threadIdx.x] += s_Data[threadIdx.x + 1];
		__syncthreads();
		if(warp_index < 30) s_Data[threadIdx.x] += s_Data[threadIdx.x + 2];
		__syncthreads();
		if(warp_index < 28) s_Data[threadIdx.x] += s_Data[threadIdx.x + 4];
		__syncthreads();
		if(warp_index < 24) s_Data[threadIdx.x] += s_Data[threadIdx.x + 8];
		__syncthreads();
		if(warp_index < 16) s_Data[threadIdx.x] += s_Data[threadIdx.x + 16];
		__syncthreads();
		
		double intra_warp_result = s_Data[threadIdx.x];
		__syncthreads();
		
		if(warp_index == 0) {
			s_Data[threadIdx.x >> LOG2_WARP_SIZE] = intra_warp_result;
		}
		__syncthreads();
		
		if(threadIdx.x < 4) {
			double my_val = s_Data[threadIdx.x];
			
			if(threadIdx.x < 3) s_Data[threadIdx.x] += s_Data[threadIdx.x + 1];
			if(threadIdx.x < 2) s_Data[threadIdx.x] += s_Data[threadIdx.x + 2];
			s_Data[threadIdx.x] -= my_val;
			
			//printf("%d %e \n", threadIdx.x, s_Data[threadIdx.x]);
		}
		__syncthreads();
		
		return s_Data[threadIdx.x >> LOG2_WARP_SIZE] + intra_warp_result; 
	}

// FORWARD SCAN KERNALS
//-------------------------------------------------------------------------------------------------
	// Almost the same as naive scan1Inclusive, but doesn't need __syncthreads()
	// assuming: size <= WARP_SIZE
	inline __device__ double warpScanInclusive(double idata, volatile double *s_Data, uint size){
		// this genrates a a sequence of [0,size][2*(size+1),(3*size-1)]
		
		/*
		 If 'size' is a power of 2 such that 'size' = 2^p then ('size'-1) is a bit pattern of all 1's for the p least significant bits. and all zeros for all
		    higher bits (for example if the '=size of an interger is 2^32 and 'size' is 2^8 then ('size'-1) = 00000000000000000000000011111111
		    so '(threadIdx.x & (size - 1)' will effectivly cause all bits in 'threadIdx.x ' above p to be zeroed and will pass all lower p bits of 'threadIdx.x ' 
		    to be passed though unaltered. This is another way to write '(threadIdx.x mod 'size')' which gives a sequence of [0....(size-1)] repeated.
		   (we have to so this because Mod is pathalogically slow on Cuda at this point... cuda is terrible at integer aritmetic.)
		
		 So '2 * threadIdx.x - (threadIdx.x & (size - 1))' does the following mapping (assuming size==WARP_SIZE but threadIdx.x can be > WARP_SIZE):
		              0  <= threadIdx.x <  1*WARP_SIZE  maps to  0          ..( 1*WARP_SIZE - 1 )     for WARP_SIZE 32:   0   .. 31
		      WARP_SIZE  <= threadIdx.x <  2*WARP_SIZE  maps to  2*WARP_SIZE..( 3*WARP_SIZE - 1 )     for WARP_SIZE 32:   64  .. 95
		    2*WARP_SIZE  <= threadIdx.x <  3*WARP_SIZE  maps to  4*WARP_SIZE..( 5*WARP_SIZE - 1 )     for WARP_SIZE 32:   128 .. 159
		    3*WARP_SIZE  <= threadIdx.x <  4*WARP_SIZE  maps to  6*WARP_SIZE..( 7*WARP_SIZE - 1 )     for WARP_SIZE 32:   192 .. 223
		    4*WARP_SIZE  <= threadIdx.x <  5*WARP_SIZE  maps to  8*WARP_SIZE..( 9*WARP_SIZE - 1 )     for WARP_SIZE 32:   256 .. 287
		    ....
		   15*WARP_SIZE  <= threadIdx.x < 16*WARP_SIZE  maps to 30*WARP_SIZE..(31*WARP_SIZE - 1 )     for WARP_SIZE 32:   960 .. 991
		*/
		uint pos = 2 * threadIdx.x - (threadIdx.x & (size - 1));
		s_Data[pos] = 0;
		pos += size;
		s_Data[pos] = idata;
		/*
			so the above statements create a structure that looks like (for warp of 32)
			
			
			buffer index: |---000 to 031---|---032 to 063---|---064 to 095---|---096 to 127---|---128 to 159---|---160 to 191---|
			      values: |        0       |vals 000 to 031 |        0       |vals 032 to 063 |        0       |vals 064 to 095 |
			
			where "vals 000 to 031" means the values passed in by the threads with 000 < threadIdx.x < 031 written in consecutive memory locations.
			
		*/
		
		// in this routine size <= WARP_SIZE...
		// the following statements are a bit magical... because they are GAURENTEED to executed simultaneously. It's kinda freeky.
		// Assume size == 32:  [X] = current value of index location X
		//                     |x| = value that was at index x at start of loop
		//                     {x} = value passed ib by thread with threadIdx.x = x 
		//                     |x..y| = sum all values inclusive of x any y that were there at start of loop
		//                     {x..y} = sum all values inclusive of values passed in thread with threadIdx.x={x to y}
		/*
			1st time though for loop (offset == 1):
				index[32] = |32| + |31| == {0} + 0
				index[33] = |33| + |32| == {1} + {0}
				index[34] = |34| + [33] == {2} + {1}
				index[35] = |35| + [34] == {3} + {2}
				index[36] = |36| + [35] == {4} + {3}
				index[62] = |62| + [61] == {30} + {29}
				index[63] = |63| + |62| == {31} + {30}
			2nd time though for loop (offset == 2):
				index[32] = [32] + [30] = (|32| + |31|) + (0) == {0} + 0 + 0 
				index[33] = [33] + [31] = (|33| + |32|) + (0) == {0} + {1} + 0 
				index[34] = [34] + [32] = (|34| + |33|) + (|32| + |31|) == {2} + {1} + {0} + 0
				index[35] = [35] + [33] = (|35| + [34]) + (|33| + |32|) == {3} + {2} + {1} + {1}
				index[36] = [36] + [34] = (|36| + [35]) + (|34| + [33]) == {4} + {3} + {2} + {1}
				index[62] = [62] + [60] = (|62| + [61]) + (|60| + [59]) == {30} + {29} + {28} + {27}
				index[63] = [63] + [61] = (|63| + |62|) + (|62| + [61]) == {31} + {30} + {29} + {28}
			3rd time through for loop (offset == 4):
				index[32] = [32] + [28] = [(|32| + |31|) + (0) ] + (0) == [{0} + 0 + 0 + 0] + [0 + 0 + 0 + 0]
				index[33] = [33] + [29] = [(|33| + |32|) + (0) ] + (0) == [{1} + {0} + 0 + 0] + [0 + 0 + 0 + 0]
				index[34] = [34] + [30] = [(|34| + |33|) + (|32| + |31|)] + [(0)] = [{2}+{1}+{0} + 0] + [0 + 0 + 0 + 0]
				index[35] = [35] + [31] = [(|35| + [34]) + (|33| + |32|)] + [(0)] = [{3}+{2}+{1}+{0}] + [0 + 0 + 0 + 0]
				index[36] = [36] + [32] = [(|36| + [35]) + (|34| + [33])] + [(|32| + |31|) + (0)]  = [{4}+{3}+{2}+{1}]+[{0} + 0 + 0 + 0]
				index[62] = [62] + [58] = [(|62| + [61]) + (|60| + [59])] + [(|58| + [57]) + (|56| + [55])] = [{30} + {29} + {28} + {27}] + [{26} + {25} + {24} + {23}]
				index[63] = [63] + [59] = [(|63| + |62|) + (|62| + [61])] + [(|59| + [58]) + (|57| + [56])] = [{31} + {30} + {29} + {28}] + [{27} + {26} + {25} + {24}]
			4th time through for loop (offset == 8):
				index[32] = [32] + [24] = [32] + 0 = {0}
				index[33] = [33] + [25] = [33] + 0 = {0..1}
				index[34] = [34] + [26] = [34] + 0 = {0..2}
				index[35] = [35] + [27] = [35] + 0 = {0..3}
				index[36] = [36] + [28] = [36] + 0 = {0..4}
				index[62] = [62] + [54] = [|62..55|] + [|54..47|] = {30..23} + {22..15} = {30..15}
				index[63] = [63] + [55] = [|63..56|] + [|55..48|] = {31..24} + {23..16} = {31..16}
			5th time through for loop (offset == 16):
				index[32] = {0}
				index[33] = {0..1}
				index[34] = {0..2}
				index[35] = {0..3}
				index[36] = {0..4}
				index[62] = [62] + [46] = [|62..47|] + [|46..30|] = {30..15} + {14..0} = {30..0}
				index[63] = [63] + [47] = [|63..48|] + [|47..31|] = {31..16} + {15..0} = {31..0}
			Done.
		*/
		for(uint offset = 1; offset < size; offset <<= 1)
			s_Data[pos] += s_Data[pos - offset];
		
		return s_Data[pos];
	}
	
	inline __device__ double warpScanExclusive(double idata, volatile double *s_Data, uint size){
		return warpScanInclusive(idata, s_Data, size) - idata;
	}
	
	inline __device__ double scan1Inclusive(double idata, volatile double *s_Data, uint size) {
		
		if(size > WARP_SIZE) {
			//Bottom-level inclusive warp scan
			double warpResult = warpScanInclusive(idata, s_Data, WARP_SIZE);
			
			// first wait for all treads to finish
			// since s_data is being written.. 
			// and we need it in it's final state before moving on.
			__syncthreads();
			
			// --------------------------------------------------------
			//Save top elements of each warp for exclusive warp scan
			//
			// this if statement is only true for the 'top' threadIdx.x in each warp
			//    so WARP_SIZE=32, the if is true for threads with threadIdx.x = [31,63, 95, 127, 159, 191, 223, 255, 287, 319, 351, 383, 415, 447, 479, 511]
			// and writes the value to index in s_data of                      = [0 , 1,  2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15]
			if( (threadIdx.x & (WARP_SIZE - 1)) == (WARP_SIZE - 1) )
				s_Data[threadIdx.x >> LOG2_WARP_SIZE] = warpResult;
			
			//wait for warp scans to complete
			__syncthreads();
			
			// Done saving the 'top' warp results to the s_data buffer
			// --------------------------------------------------------
			
			
			//--------------------------------------------------------
			// Add the 'top' results to get final result.
			//
			// for all threads with threadinx.x = 0 to 15  [one thread for every 'top' result written in previous result].
			// Do another 'warp scan' but with values pass in that are the sum of a Warp.
			if( threadIdx.x < (THREADBLOCK_SIZE / WARP_SIZE) ) {
				//grab top warp elements
				double val = s_Data[threadIdx.x];
				//printf("FUCK %d: %e \n", threadIdx.x,  val);
				//calculate exclsive scan and write back to shared memory
				s_Data[threadIdx.x] = warpScanExclusive(val, s_Data, size >> LOG2_WARP_SIZE);
			}
			
			/* so now the first 16 elements (THREADBLOCK_SIZE / WARP_SIZE elements in general... 16 for 512/32) of s_data have:
				s_data[0] =  warps(0..0) =  {0..31}
				s_data[1] =  warps(0..1) =  {0..63}
				s_data[2] =  warps(0..2) =  {0..95}
				s_data[3] =  warps(0..3) =  {0..127}
				s_data[4] =  warps(0..4) =  {0..159}
				s_data[14] = warps(0..14) = {0..479}
				s_data[15] = warps(0..15) = {0..511}
			*/
			// wait for all write to s_data to be done bewfore moving on.
			__syncthreads();
			
			// Done adding top results...
			//--------------------------------------------------------
			
			
			// Now all threads partisipate again...
			// The final result passed back is the "preceeding warp sum" plus the sum within warp that thread is in.!
			
			
			return warpResult + s_Data[threadIdx.x >> LOG2_WARP_SIZE];
			
		} else {
			
			return warpScanInclusive(idata, s_Data, size);
			
		}
	}
	
	inline __device__ double scan1Exclusive(double idata, volatile double *s_Data, uint size) {
		return scan1Inclusive(idata, s_Data, size) - idata;
	}
//-------------------------------------------------------------------------------------------------


inline __device__ double4 scan4Inclusive(double4 idata4, volatile double *s_Data, uint size) {
	
	//Level-0 inclusive scan
	idata4.y += idata4.x;
	idata4.z += idata4.y;
	idata4.w += idata4.z;
	
	//Level-1 exclusive scan
	double oval = scan1Exclusive(idata4.w, s_Data, size / 4);
	
	idata4.x += oval;
	idata4.y += oval;
	idata4.z += oval;
	idata4.w += oval;
	
	return idata4;
}

inline __device__ double4 scan4Exclusive(double4 idata4, volatile double *s_Data, uint size){
	double4 odata4 = scan4Inclusive(idata4, s_Data, size);
	odata4.x -= idata4.x;
	odata4.y -= idata4.y;
	odata4.z -= idata4.z;
	odata4.w -= idata4.w;
	return odata4;
}

__global__ void scanInclusiveShared(
    double4 *d_Dst,
    double4 *d_Src,
    uint size
){
    __shared__ double s_Data[2 * THREADBLOCK_SIZE];

    uint pos = blockIdx.x * blockDim.x + threadIdx.x;

    //Load data
    double4 idata4 = d_Src[pos];

    //Calculate exclusive scan
    //uint4 odata4 = scan4Exclusive(idata4, s_Data, size);
	double4 odata = scan4Inclusive(idata4, s_Data, size);

    //Write back
    d_Dst[pos] = odata;
}

__global__ void scanInclusiveShared_one_at_a_time(	double *d_Dst,
													double *d_Src,
													uint size){
	__shared__ double s_Data[2 * THREADBLOCK_SIZE];
	uint pos = blockIdx.x * blockDim.x + threadIdx.x;
	double data = d_Src[pos];
	double final_val = scan1Inclusive(data, s_Data, size);

	//data = 1.4;
	double final_val2 = scan1Inclusive_rev(data, s_Data, size);
	d_Dst[pos] = final_val2;

}

__global__ void scanInclusiveShared_one_at_a_time_rev(	double *d_Dst,
													double *d_Src,
													uint size){
	__shared__ double s_Data[2 * THREADBLOCK_SIZE];
	uint pos = blockIdx.x * blockDim.x + threadIdx.x;
	
	double data = d_Src[pos];
	double final_val = scan1Inclusive_rev(data, s_Data, size);
	d_Dst[pos] = final_val;

}


void scan_test() {
	int size = 512;
	int size_logical = 512;

	// cheat for now and bump up mem size so eveything is even.
	int blocks_needed = size / (THREADBLOCK_SIZE);
	if((size % (THREADBLOCK_SIZE)) != 0) {
		blocks_needed += 1;
	}
	size_logical = blocks_needed *THREADBLOCK_SIZE;

	double * host_buffer = new double[size];
	double * gpu_buffer_src = 0;
	double * gpu_buffer_dest = 0;

	for(int i=0;i<size;i++) {
		host_buffer[i] = 3.0;
	}

	cout << "Hi from scan test (one element at a time and REV)" << endl;
	cudasafe(cudaMalloc((void**)&(gpu_buffer_src),size * sizeof(double)), 			"while mallocing on the gpu: __FILE__ __LINE__");
	cudasafe(cudaMalloc((void**)&(gpu_buffer_dest),size * sizeof(double)), 			"while mallocing on the gpu: __FILE__ __LINE__");

	// copy data to gpu array
	cudaMemcpy(gpu_buffer_src, host_buffer,  sizeof(double)*size, cudaMemcpyHostToDevice); 

	for(int i=0;i<size;i++) {
		host_buffer[i] = 0.0;
	}

	// this is linear 1-d grid array.
	dim3 dimBlock(size);
	dim3 dimGrid(blocks_needed);
	
	// call the kernal...
	// it requires an array that is power of 2.
	cout << "Calling kddernal with " << size << " and " << blocks_needed << endl;
	//scanInclusiveShared_one_at_a_time_rev<<<dimGrid, dimBlock>>>( gpu_buffer_dest, gpu_buffer_src, size);
	scanInclusiveShared_one_at_a_time<<<dimGrid, dimBlock>>>( gpu_buffer_dest, gpu_buffer_src, size);
	cutilCheckMsg("Fuck up");

	// copy results back to host's fortran array.
	//cudaMemcpy(parts.species[i].x, parts.species[i].cudaX,  sizeof(double)*parts.species[i].num_par_max, cudaMemcpyDeviceToHost); 
	cudaMemcpy(host_buffer, gpu_buffer_dest, sizeof(double)*size, cudaMemcpyDeviceToHost); 
	for(int i=0;i<size;i++) {
		if (i%16==0) cout << endl;
		cout << host_buffer[i] << " ";
	}
	cout << endl;
}

void scan_test_4_at_a_time() {
	int size = 64;
	int size_logical = 64;

	// cheat for now and bump up mem size so eveything is even.
	int blocks_needed = size / (4*THREADBLOCK_SIZE);
	if((size % (4*THREADBLOCK_SIZE)) != 0) {
		blocks_needed += 1;
	}
	size_logical = blocks_needed * 4*THREADBLOCK_SIZE;

	double * host_buffer = new double[size];
	double * gpu_buffer_src = 0;
	double * gpu_buffer_dest = 0;

	for(int i=0;i<size;i++) {
		host_buffer[i] = 1.0;
	}

	cout << "Hi from scan test!" << endl;
	cudasafe(cudaMalloc((void**)&(gpu_buffer_src),size * sizeof(double)), 			"while mallocing on the gpu: __FILE__ __LINE__");
	cudasafe(cudaMalloc((void**)&(gpu_buffer_dest),size * sizeof(double)), 			"while mallocing on the gpu: __FILE__ __LINE__");

	// copy data to gpu array
	cudaMemcpy(gpu_buffer_src, host_buffer,  sizeof(double)*size, cudaMemcpyHostToDevice); 

	for(int i=0;i<size;i++) {
		host_buffer[i] = 0.0;
	}

	// this is linear 1-d grid array.
	dim3 dimBlock(size);
	dim3 dimGrid(blocks_needed);
	
	// call the kernal...
	cout << "Calling kernal with " << size << " and " << blocks_needed << endl;
	scanInclusiveShared<<<dimGrid, dimBlock>>>( (double4 *) gpu_buffer_dest, (double4 *) gpu_buffer_src, size);
	cutilCheckMsg("Fuck up");

	// copy results back to host's fortran array.
	//cudaMemcpy(parts.species[i].x, parts.species[i].cudaX,  sizeof(double)*parts.species[i].num_par_max, cudaMemcpyDeviceToHost); 
	cudaMemcpy(host_buffer, gpu_buffer_dest, sizeof(double)*size, cudaMemcpyDeviceToHost); 
	for(int i=0;i<size;i++) {
		if (i%20==0) cout << endl;
		cout << host_buffer[i] << " ";
	}
	cout << endl;
}
//__constant__ int NB;
//__constant__ double density_np;


/*
#######################################################################################################
#######################################################################################################
	6th attempt..... try to save some shared mem
#######################################################################################################
#######################################################################################################
*/

// This is the unsplit evaluation step.
 __device__ inline void eval_f_slopey_v6(double * fc, double * dest, double * logglam, const double *vr,
							//double const &_vr, double  const &_vrm1, double  const &_vrp1, 
							const uint size, const int padded_size, const int debug, volatile double *s_Data, const uint pos, const int eval_ind) {

	// these should be a 'synch' called before this routine is invoked.
	// the synch is needed bcuase data needs to be shared between the threads via the shared mem..
	double _vr = vr[id];
	double _vrm1 = 0;
	double _vrp1 = 0;
	if(id > 0) _vrm1 = vr[id-1];
	if(id < size-1) _vrp1 = vr[id+1];



	// load our data value in from fc (the shared mem chached input data)
	__syncthreads();
	
	double data = fc[id];
	double datam1 = 0.0;					// TODO: this needs to be fixed.I just wanna allocate the buffer to be a bit bigger do pos-1 won't seg fault at 0.
	if(id > 0) datam1 = fc[id-1];			

	//printf("%e %e %e %e %e\n",_vr, _vrm1,_vrp1,data, datam1)  ;
	double t = .5*_vr*_vr*_vr*_vr;
	double tm1 = .5*_vrm1*_vrm1*_vrm1*_vrm1;//*(_vr - _vrm1);
	double reduce_me = t*data + tm1*datam1;
	reduce_me *= (_vr - _vrm1);
	if(id == 0) reduce_me = 0.0;

	double final_I4 = scan1Inclusive(reduce_me, s_Data, padded_size);
	

																					//if(bi eval_num id == 0) printf("\n\tI4[n] inital: \n\t\t");
																					//if(bi eval_num id < 5) printf("%d:%e ", id, final_I4);
	t = .5*_vr*_vr;
	tm1 = .5*_vrm1*_vrm1;
	reduce_me = t*data + tm1*datam1;
	reduce_me *= (_vr - _vrm1); 
	if(id == 0) reduce_me = 0.0;
	double final_I2 = scan1Inclusive(reduce_me, s_Data, padded_size);
																					//if(bi eval_num id == 0) printf("\n\tI2[n] inital: \n\t\t");
																					//if(bi eval_num id < 5) printf("%d:%e ", id, final_I2);
	
	t = .5*_vrp1;
	tm1 = .5*_vr;
	reduce_me = t*fc[id+1] + tm1*data;
	

	reduce_me *= (_vrp1 - _vr);
	if(id == 31) reduce_me = 0.0;
																					//if(bi eval_num id == 0) printf("\n\tJ1 input: \n\t\t");
																					//if(bi eval_num id < 5) printf("%d:%e ", id, reduce_me);


	double final_J1 = scan1Inclusive_rev(reduce_me, s_Data, padded_size);
	__syncthreads();
	
																					//if(	bi eval_num id == 0) printf("\n\tJ1[n] input: \n\t\t");
																					//if(bi eval_num id < 5) printf("%d:%e ", id, final_J1);

	if( id  == (size-1)) {
		*logglam = LOGee_cuda(4.0*M_PI*final_I2, 4.0*M_PI*final_I4);
																					//if(	bi eval_num) printf("\n\tLoglam\n\t\t%e", (*logglam));
	}
					
	final_I2 *= final_J1;
	
	final_I2 *= -3.0;
	final_I4 += _vr*_vr*_vr * final_J1;
	final_I4 *= data;
	final_I4 += final_I2;


	
	__syncthreads();
																					//if(	bi eval_num id == 0) printf("\n\tG inital: \n\t\t");
																					//if(bi eval_num id < 5) printf("%d:%e ",id,  final_I4);


	if(threadIdx.x < NB) {
		final_I4 = G_cuda(threadIdx.x, fc, vr, final_J1);
	}
																					//if(	bi eval_num id == 0) printf("\n\tG after bound: \n\t\t");
																					//if(bi eval_num id < 5) printf("%d:%e ", id, final_I4);

	// add another divergencew to take care of first cell...
	if(threadIdx.x == 0) {
		// store answer in u
		reduce_me = -1.0 * final_I4;
		reduce_me *= 2.0 / (vr[0]*vr[0]);
		dest[0] = reduce_me;
	}


	__syncthreads();
	s_Data[id] = final_I4;
	tm1 = s_Data[id+1];		// um1 is now I4[n+1]

	final_I4 -= tm1;			// I4[n] -= I4[n+1];

																					//if(	bi eval_num id == 0) printf("\n\tFirst 14 diff \n\t\t");
																					//if(bi eval_num id < 5) printf("%d:%e ",id,  final_I4);

	//u = Pn[id];
	
	t = ((_vrp1-_vr)*(_vrp1+_vr))/2.0; ///2.0/((_vrp1-_vr)*(_vrp1+_vr))
	
	//__syncthreads();
	if(t > 0)
		final_I4 /= t;				// I4[n] *= Pn[n];
																					//if(bi eval_num id == 0 && debug==1) printf("\n\tAfter Pn[n]:\n\t\t");
																					//if(bi eval_num id < 5 && debug==1) printf("%d:%e ", id, final_I4);
	s_Data[id] = final_I4;
	tm1 = s_Data[id+1];		// um1 is now I4[n+1]

	// this is terrible.. everyone will have to wait till this is done.
	if(threadIdx.x == 0) {
		reduce_me -= final_I4;
		reduce_me *= 2.0/(_vr*_vr*_vrp1);
		dest[0] = reduce_me;
	}
	__syncthreads();


	final_I4 -= tm1;			// I4[n] -= I4[n+1];
																					//if(bi eval_num id == 0 && debug==1) printf("\n\tSecond I4 diff\n\t\t");
																					//if(bi eval_num id < 5 && debug==1) printf("%d:%e ", id, final_I4);
	//u = Qn[id+1];

//1.0 / (vr[i+1]*vr[i+1]*(vr[i+2]-vr[i])/2.0);



	//t = ((_vrp1 - _vrm1)*_vr*_vr)/2.0;    //2.0/((_vrp1 - _vrm1)*_vr*_vr);
	t = (_vrp1*_vrp1*(vr[id+2]-_vr))/ 2.0;
	if(t > 0)
		final_I4 /= t;				// final_I4 = Qn[n+1]*I4[n];
	__syncthreads();

							 														//if(bi eval_num id == 0 && debug==1) printf("\n\tBefore Norm\n\t\t");
																					//if(bi eval_num id < 5 && debug==1) printf("%d:%e ", id, final_I4);		
	// TODO: prolly should put log ee in a register otr shared mem.
	// dp the normalization...

	final_I4 *= (*logglam);
	final_I4 *= c_kpre;

	dest[id+1] = final_I4;

	if(threadIdx.x == 0) {
		dest[0] *= (*logglam) *c_kpre;
	}
							 														//if(bi eval_num id == 0 && debug==1) printf("\n\After\n\t\t");
																					//if(bi eval_num id < 5 && debug==1) printf("%d:%e ", id, final_I4);		

}


__device__  double G_cuda_v6(const int n, 
							const double data0, const double data1, const double datan, 
							double vr0,double vr1, double vrn,
							double J1) {
	double i2s, i4s;
	// right now, everyone reads f00 and passes it in a register.. I think shared broadcasting may be faster.
	// read (vr0*vr0)/(vr1*vr1) from constant mem (which also boradcasts)
	
	double f00( (data0 - data1*(vr0*vr0)/(vr1*vr1))/ (1.0 - (vr0*vr0)/(vr1*vr1)) );
	i2s = f00*pow(vrn,3)/3.0 + (data1-f00)*pow(vrn,5)/(vr1*vr1)*0.2;
	i4s = f00*pow(vrn,5)*0.2 + (data1-f00)*pow(vrn,7)/(vr1*vr1*7.0);
	return datan*i4s + (pow(vrn,3)*datan-3.0*i2s) * J1;
	/*
	double i2s, i4s;
	double f00( (fin[0] - fin[1]*(vr[0]*vr[0])/(vr[1]*vr[1]))/ (1.0 - (vr[0]*vr[0])/(vr[1]*vr[1])) );
	//printf("\tf00 %e\n", f00);
	i2s = f00*pow(vr[n],3)/3.0 + (fin[1]-f00)*pow(vr[n],5)/(vr[1]*vr[1])*0.2;
	//printf("\ti2s %e\n", i2s);
	i4s = f00*pow(vr[n],5)*0.2 + (fin[1]-f00)*pow(vr[n],7)/(vr[1]*vr[1]*7.0);
	//printf("\ti4s %e\n", i4s);
	return fin[n]*i4s + (pow(vr[n],3)*fin[n]-3.0*i2s) * J1;*/

}



//#define DBROOT blockIdx.x == 0 && threadIdx.x==0 && debug == 1

template <unsigned int rk4_step>
 __global__ void rk4__ee_v6_sortee_PHYSCO(	double * __restrict__ fin, 
											double * __restrict__ trial_input_vector_to_use_in_evalation, 
											double * __restrict__ y_new, 
											const double * __restrict__ vr, 
											double * __restrict__ cell_data__loglam, 
											double * __restrict__ cell_data__precomp1,
											double * __restrict__ cell_data__precomp2,
											double * __restrict__ cell_data__precomp3,
											uint size, int padded_size, int debug, double h, int numh, double rk4_trial_factor, double rk4_factor) {	
	
	// declare needed shared mem.
	//__shared__ double swap_space[THREADBLOCK_SIZE];	
	__shared__ double swap_space[128];	
	__shared__ double zeroth;
	
																								LOG_ARRAY_INT("Reviving in physvo at RK4 step rk4_step with pos: ", debug, "\t");
	/* ------------------------------------------------- */
	/* --- restart and resume place boilerplate--------- */
	uint pos = blockIdx.x* blockDim.x + threadIdx.x;
	
	double final_I4 = cell_data__precomp1[pos];
	double _vr = vr[id];
	double _vrp1 = vr[id+1];	// we are ok with _vr because it havily zero padded wht only 1 copy in whole..
	double _vrp2 = vr[id+2];
	
	double final_I4_p1 = 0.0;
	if(id < (size-1)) {
		final_I4_p1 = cell_data__precomp1[pos+1];
	}
	/* --- END restart and resume place boilerplate----- */
	/* ------------------------------------------------- */
																								LOG_ARRAY("final I4 (results from prev. steps loaded from global mem):", final_I4, "\t");
																						
	
	// take care of the zero compoenent...
	if(threadIdx.x == 0) {
		zeroth = -1.0 * final_I4;
		zeroth *= 2.0 / (_vr*_vr);
	}
	
	final_I4 -= final_I4_p1;
																								LOG_ARRAY1("Result of first final_I4 diff:", zeroth, final_I4, "\t");
	
	
	double t = ((_vrp1-_vr)*(_vrp1+_vr))/2.0;  //(this is p[n])
	if(t > 0)
		final_I4 /= t;
																								LOG_ARRAY("After p[n]: ", final_I4, "\t");
	
	// This is Find DDG/(vDv)
	/*t= ((_vrp2-_vrp1)*(_vrp2+_vrp1))/2.0;  //(this is p[n+1])
	if(t > 0)
		final_I4 /= t;	
	__syncthreads();
	*/
	
	// Now take a forward derivitive (  I4[n] -= I4[n+1] )
	swap_space[id] = final_I4;
	__syncthreads();
	final_I4 -= swap_space[id+1];
																								LOG_ARRAY("After second I4 diff: ", final_I4, "\t");
	
	
	// final_I4 = Qn[n+1]*I4[n];
	t = (_vrp1*_vrp1*(_vrp2-_vr))/ 2.0;
	if(t > 0)
		final_I4 /= t;
	__syncthreads();
																								LOG_ARRAY("tBefore Noralization ", final_I4, "\t");
	
	final_I4 *= c_kpre;  // finishing the normalizing process..
																								LOG_ARRAY1("tFINAL EVAL ", rk4_trial_factor, final_I4, "\t");
	
	
	// So final_I4 contains the results from an evaulation of the F function of rk4.
	// But there is a slight wrinkle... it holds the result not of this 
	// position in the array, but the n+1 space in the array.
	// KILL THIS... it's just for debugging.
	t =  final_I4*rk4_trial_factor;
																								LOG_ARRAY1("eval*factor h/2: ", rk4_trial_factor, t, "\t");
	
	
	if(id < (size-1)) {
		/* use the answer we just computed and combine it with the original input vector
		   to use as the inital evulation input vector for next rk4 substep. */
		t = fin[pos+1];
		//t+= final_I4*rk4_trial_factor;
		//trial_input_vector_to_use_in_evalation[pos] = t
		trial_input_vector_to_use_in_evalation[pos+1] = t + final_I4*rk4_trial_factor;
		
		/* Also, we need to add this step into the overall answer all rk4 steps.
		   In the case of the 1st evaulation step, we need to prime the y_new vector.. so ass the result from this
		   step to the original input vector ams write to 'ynew'. 
		   Otherwise we simply add the result of this step to 'ynew'*/
		if(rk4_step == 1) {
			y_new[pos+1] = t + (final_I4*rk4_factor);
		} else if(rk4_step == 4) {
			t = y_new[pos+1];
			t += (final_I4*rk4_factor);
			fin[pos+1] = t;
		} else {
			t = y_new[pos+1];
			t += (final_I4*rk4_factor);
			y_new[pos+1] = t;
		}
	}
	__syncthreads();
	
	
	// I can do better then this..... maybe move it up to the last thread?
	// This is occuring because the answer in a thread is really the final answer for the [n+1] position..
	// so the zeroth thread gets left out... to really preserve symmetry, this calculation should really go into
	//   the last thread.. since that thread is not doing any usefull work (it's result get thrown out).
	if(threadIdx.x == 0) {
		zeroth -= swap_space[0]; // final_I4;
		zeroth *= 2.0/(_vr*_vr*_vrp1) * c_kpre;
		
																								LOG_ROOT("\tFinal 0th component this eval:        %e\n\t\t", zeroth);
		t = fin[pos];
		trial_input_vector_to_use_in_evalation[pos] = t + zeroth*rk4_trial_factor;
		if(rk4_step == 1) {
			y_new[pos] = fin[pos] + (zeroth*rk4_factor);
		} else if(rk4_step == 4) {
			t = y_new[pos];
			t += (zeroth*rk4_factor);
			fin[pos] = t;
		} else {
			t = y_new[pos];
																								LOG_ROOT("\tSaved 0th val before timestep update: %.8e \n\t\t", t);
			t += (zeroth*rk4_factor);
																								LOG_ROOT("\tSaved 0th val after timestep update:  %.8e \n\t\t", t);
			y_new[pos] = t;
		}
	}
	
	//__syncthreads();
	// KILL THIS... it's just for debugging.
	//t = trial_input_vector_to_use_in_evalation[pos];
	//																							LOG_ARRAY("The Starting values for next RK4 itteration sucker", t, "\t");
	//
	// If this is the 4th (and final) eval step, then copy results
	// to the results (input) buffer. [This is just copying.. so we can 
	// assume that every thread[n] has data for position[n]]
	// NOTE: this 'if(rk4_step==4)' is templated... so has no run-time cost.
	//if(rk4_step==4) {
	//	__syncthreads();
	//	t = y_new[pos];
	//	fin[pos] = t;
	//																							LOG_ROOT("-----------------------------------\n\t----------- Done with an RK4 Eval -----------\n-----------------------------------\n", t);
	//}
																								LOG_ARRAY("The current values of work-in-progress RK4 solution\n", t, "\t");
 }

 
/*
//    split version1
//
*/

// Shared mem for rk4 stuff..
// explicit:gl
//	 	1 global read (double)f
//	f 	1 global write (double)
//
#define SynchThreads __syncthreads();

__device__ inline void ass_fuck1( const int padded_size,  const uint pos, double * shared_buffer, const double vr, const double vrm1, const double ddata, const double ddatam1, double * __restrict__ out) {
	double t = .5*vr*vr;
	double tm1 = .5*vrm1*vrm1;
	double reduce_me = t*ddata + tm1*ddatam1;
	reduce_me *= (vr - vrm1); 
	if(id == 0) reduce_me = 0.0;
	
	out[pos] = scan1Inclusive(reduce_me, shared_buffer, padded_size);
}
__device__ inline void ass_fuck2( const int padded_size, const uint pos, double * shared_buffer, const double vr, const double vrm1, const double ddata, const double ddatam1, double * __restrict__ out) {
	double t = .5*vr*vr*vr*vr;
	double tm1 = .5*vrm1*vrm1*vrm1*vrm1;
	double reduce_me = t*ddata + tm1*ddatam1;
	reduce_me *= (vr - vrm1); 
	if(id == 0) reduce_me = 0.0;
	
	out[pos] = scan1Inclusive(reduce_me, shared_buffer, padded_size);
}

template <unsigned int test_type>
__global__ void test_forward_prefix_scans(double * buffer) {
	if (test_type == 0) {
		__shared__ double s_Data[132];
		double shit = 1.0;
		//printf("FUCK 0 \n");
		double crap = warpScanInclusiveSimpleWarp(shit, s_Data, 128) ;
		buffer[threadIdx.x] = crap;
	} else if (test_type == 1) {
		__shared__ double s_Data[132];
		double shit = 1.0;
		//printf("FUCK 1 \n"); 
		double crap = warpScanInclusiveSimple(shit, s_Data, 128) ;
		buffer[threadIdx.x] = crap;
	} else if (test_type == 2) {
		__shared__ double s_Data[132];
		double shit = 1.0;
		//printf("FUCK 2 \n");
		double crap = warpScanInclusiveSimpleWarp_rev(shit, s_Data, 128) ;
		buffer[threadIdx.x] = crap;
		//printf("%d %e \n", threadIdx.x, crap);
	} else {
		
		__shared__ double s_Data[2 * THREADBLOCK_SIZE];
		double shit = 1.0;
		//printf("FUCK 3 \n");
		double crap = scan1Inclusive(shit, s_Data, 128);
		buffer[threadIdx.x] = crap; 
	}
}


__global__ void rk4__ee_v6_sortee_FAGGOT_split3__1(	double * __restrict__ fin, double * __restrict__ dest, const double * __restrict__ vr, 
								    double * __restrict__ cell_data__loglam, 
									double * __restrict__ cell_data__precomp1,
									double * __restrict__ cell_data__precomp2,
									double * __restrict__ cell_data__precomp3,
									double * __restrict__ cell_data__precomp4,
									uint size, int padded_size, int debug, double h, int numh,
									const double * __restrict__ U4,
									const double * __restrict__ U4m1,
									const double * __restrict__ U2,
									const double * __restrict__ U2m1) {	
	
	//__shared__ double s_Data[2 * THREADBLOCK_SIZE];
	__shared__ double s_Data[128];
	//__shared__ double logLambda;
	
	// Note: when comparing values in this funtion with the CPU e-e collsion code, keep in mind these values will be different
	//  by a factor of LogLambda (in the cpu code LogLombda is multiplied in the last step.. in the GPU code, it is mulitped eailier in the process)
	
	/* ------------------------------------------------- */
	/* --- restart and resume place boilerplate--------- */
	uint pos = blockIdx.x* blockDim.x + threadIdx.x;
	
	double _vr = vr[id];										/* get momentum basis value for this cell. */
	double _vrm1 = 0;											/* get vr for point one behind. */
	if(id > 0) _vrm1 = vr[id-1];
	//double _vrp1 = 0;											/* get vr for point one ahead. */
	//if(id < size-1) _vrp1 = vr[id+1];
	
	double data = fin[pos];										/* get data point. */
	double datam1 = 0.0;										/* get data point for one position behind. */
	if(id > 0) datam1 = fin[pos-1];
	//double datap1 = 0.0;										/* get data point for one position ahead. */
	//if(id < size-1) datap1 = fin[pos+1];
	
	/* --- END restart and resume place boilerplate----- */
	/* ------------------------------------------------- */
																								LOG_ARRAY("Incoming Data to New Eval Step", data, "\t"); 
	
	/* now set up the calculation */
	
	/*
		// try functions...
		ass_fuck1(padded_size, pos, s_Data, _vr, _vrm1, data, datam1, cell_data__precomp4);
		//ass_fuck2(padded_size,pos, s_Data, _vr, _vrm1, data, datam1, cell_data__precomp1);
		double t = .5*_vr*_vr*_vr*_vr;
		double tm1 = .5*_vrm1*_vrm1*_vrm1*_vrm1;
		double reduce_me = t*data + tm1*datam1;
		reduce_me *= (_vr - _vrm1); 
		if(id == 0) reduce_me = 0.0;
		cell_data__precomp1[pos] = scan1Inclusive(reduce_me, s_Data, padded_size);
	*/
	
	/*
	// does the 1,2 clac at once using prcomed values (uses 19 registers)
	double t = U2[threadIdx.x];
	double tm1 = U2m1[threadIdx.x];
	double reduce_me = t*data + tm1*datam1;
	if(id == 0) reduce_me = 0.0;
	cell_data__precomp4[pos] = scan1Inclusive(reduce_me, s_Data, padded_size);
	
	t = U4[threadIdx.x];
	tm1 = U4m1[threadIdx.x];
	reduce_me = t*data + tm1*datam1;
	if(id == 0) reduce_me = 0.0;
	cell_data__precomp1[pos] = scan1Inclusive(reduce_me, s_Data, padded_size);
	*/
	
	
	/*
	// does the 1,2 calc at once. tries to resue variables (uses 28 fucking registers)
	double t = .5*_vr*_vr;
	double tm1 = .5*_vrm1*_vrm1;
	double reduce_me = t*data + tm1*datam1;
	reduce_me *= (_vr - _vrm1); 
	if(id == 0) reduce_me = 0.0;
	
	cell_data__precomp4[pos] = scan1Inclusive(reduce_me, s_Data, padded_size);
	
	
	t *= _vr*_vr;
	tm1 *= _vrm1*_vrm1;
	reduce_me = t*data + tm1*datam1;
	reduce_me *= (_vr - _vrm1);
	if(id == 0) reduce_me = 0.0;
	cell_data__precomp1[pos] = scan1Inclusive(reduce_me, s_Data, padded_size);
	*/
	
	/*
	// does the 1,2 calc at once. tries to resue variables another way(uses 28 fucking registers)
	double t = .5*_vr*_vr;
	double tm1 = .5*_vrm1*_vrm1;
	double reduce_me = t*data + tm1*datam1;
	reduce_me *= (_vr - _vrm1); 
	if(id == 0) reduce_me = 0.0;
	//cell_data__precomp4[pos] = scan1Inclusive(reduce_me, s_Data, padded_size);
	cell_data__precomp4[pos] = warpScanInclusiveSimpleWarp(reduce_me, s_Data, padded_size);
	*/
	
	double t = .5*_vr*_vr*_vr*_vr;
	double tm1 =.5*_vrm1*_vrm1*_vrm1*_vrm1;
	double reduce_me = t*data + tm1*datam1;
	reduce_me *= (_vr - _vrm1);
	if(id == 0) reduce_me = 0.0;
	//cell_data__precomp1[pos] = scan1Inclusive(reduce_me, s_Data, padded_size);
	cell_data__precomp1[pos] = warpScanInclusiveSimpleWarp(reduce_me, s_Data, padded_size);
	
	
	//t = scan1Inclusive(reduce_me, s_Data, padded_size);
	//cell_data__precomp4[pos] = 43.9;
	
	
	//cell_data__precomp1[pos] = 42.0;
	
	
	/*
	//  does the 1,2 calc at once. renames variables. (uses 28 fucking registers)
	
	double t = .5*_vr*_vr;
	double tm1 = .5*_vrm1*_vrm1;
	double reduce_me = t*data + tm1*datam1;
	reduce_me *= (_vr - _vrm1); 
	if(id == 0) reduce_me = 0.0;
	cell_data__precomp4[pos] = scan1Inclusive(reduce_me, s_Data, padded_size);
	
	double fuck_nvidas_shitty_compilier_t = .5*_vr*_vr*_vr*_vr;
	double fuck_nvidas_shitty_compilier_tm1 = .5*_vrm1*_vrm1*_vrm1*_vrm1;
	double fuck_nvidas_shitty_compilier_reduce_me = fuck_nvidas_shitty_compilier_t*data + fuck_nvidas_shitty_compilier_tm1*datam1;
	fuck_nvidas_shitty_compilier_reduce_me *= (_vr - _vrm1);
	if(id == 0) fuck_nvidas_shitty_compilier_reduce_me = 0.0;
	cell_data__precomp1[pos] = scan1Inclusive(fuck_nvidas_shitty_compilier_reduce_me, s_Data, padded_size);
	*/
	
	/*
	// just does 1 calc then leaves..
	// now set up the calculation 
	double t = .5*_vr*_vr*_vr*_vr;
	double tm1 = .5*_vrm1*_vrm1*_vrm1*_vrm1;
	double reduce_me = t*data + tm1*datam1;
	reduce_me *= (_vr - _vrm1);
	if(id == 0) reduce_me = 0.0;
	
	// do the 'integral' (sum) 
	double final_I4 = scan1Inclusive(reduce_me, s_Data, padded_size);
																								LOG_ARRAY("I4[n] inital", final_I4, "\t"); 
	
	
	//save result in global (persistent) mem for next subkernal 
	cell_data__precomp1[pos] = final_I4;
																								LOG_ROOT("------------END OF SPLIT3_1--------------", debug);
	*/
	
	/*
	//newer, low shared mem version of 'just one calc then leave'
	double t = .5*_vr*_vr*_vr*_vr;
	double tm1 =.5*_vrm1*_vrm1*_vrm1*_vrm1;
	double reduce_me = t*data + tm1*datam1;
	reduce_me *= (_vr - _vrm1);
	if(id == 0) reduce_me = 0.0;
	
	//double final_I4 = scan1Inclusive(reduce_me, s_Data, padded_size);
	double final_I4 = warpScanInclusiveSimpleWarp(reduce_me, s_Data, padded_size);
																								LOG_ARRAY("I4[n] inital", final_I4, "\t"); 
	
	
	//save result in global (persistent) mem for next subkernal 
	cell_data__precomp1[pos] = final_I4;
																								LOG_ROOT("------------END OF SPLIT3_1--------------", debug);
	*/
	
}

__global__ void rk4__ee_v6_sortee_FAGGOT_split3__2(	double * __restrict__ fin, double * __restrict__ dest, const double * __restrict__ vr, 
								    double * __restrict__ cell_data__loglam, 
									double * __restrict__ cell_data__precomp1,
									double * __restrict__ cell_data__precomp2,
									double * __restrict__ cell_data__precomp3,
									double * __restrict__ cell_data__precomp4,
									uint size, int padded_size, int debug, double h, int numh) {	
	
	//__shared__ double s_Data[2 * THREADBLOCK_SIZE];
	__shared__ double s_Data[128];
	
	/* ------------------------------------------------- */
	/* --- restart and resume place boilerplate--------- */
	uint pos = blockIdx.x* blockDim.x + threadIdx.x;
	
	double _vr = vr[id];										/* get momentum basis value for this cell. */
	double _vrm1 = 0;											/* get vr for point one behind. */
	if(id > 0) _vrm1 = vr[id-1];
	//double _vrp1 = 0;											/* get vr for point one ahead. */
	//if(id < size-1) _vrp1 = vr[id+1];
	
	double data = fin[pos];										/* get data point. */
	double datam1 = 0.0;										/* get data point for one position behind. */
	if(id > 0) datam1 = fin[pos-1];
	//double datap1 = 0.0;										/* get data point for one position ahead. */
	//if(id < size-1) datap1 = fin[pos+1];
	
	/* --- END restart and resume place boilerplate----- */
	/* ------------------------------------------------- */
	
	
	/* now set up the calculation */
	double t = .5*_vr*_vr;
	double  tm1 = .5*_vrm1*_vrm1;
	double  reduce_me = t*data + tm1*datam1;
	reduce_me *= (_vr - _vrm1); 
	if(id == 0) reduce_me = 0.0;
	
	/* do the 'integral' (sum) */
	double final_I2 = 	warpScanInclusiveSimpleWarp(reduce_me, s_Data, padded_size);
	//double final_I2 = scan1Inclusive(reduce_me, s_Data, padded_size);
	
	/* save result in global (persistent) mem for next subkernal */
	cell_data__precomp4[pos] = final_I2;
																								LOG_ARRAY("I2[n] inital", final_I2, "\t"); 
																								LOG_ROOT("------------END OF SPLIT3_2--------------", debug);
	
}

__global__ void rk4__ee_v6_sortee_FAGGOT_split3__2_5(	double * __restrict__ fin, double * __restrict__ dest, const double * __restrict__ vr, 
								    double * __restrict__ cell_data__loglam, 
									double * __restrict__ cell_data__precomp1,
									double * __restrict__ cell_data__precomp2,
									double * __restrict__ cell_data__precomp3,
									double * __restrict__ cell_data__precomp4,
									double * __restrict__ cell_data__precomp5,
									uint size, int padded_size, int debug, double h, int numh) {
	
	//__shared__ double s_Data[2 * THREADBLOCK_SIZE];
	__shared__ double s_Data[128];
	
	/* ------------------------------------------------- */
	/* --- restart and resume place boilerplate--------- */
	uint pos = blockIdx.x* blockDim.x + threadIdx.x;
	
	double _vr = vr[id];									/* get momentum basis value for this cell. */
	double _vrp1 = 0;										/* get vr for point one ahead. */
	if(id < size-1) _vrp1 = vr[id+1];
	
	double data = fin[pos];									/* get data point. */
	double datap1 = 0.0;									/* get data point for one position ahead. */
	if(id < size-1) datap1 = fin[pos+1];
	
	/* --- END restart and resume place boilerplate----- */
	/* ------------------------------------------------- */
	
	/* now set up the calculation */
	double t = .5*_vrp1;
	double tm1 = .5*_vr;
	double reduce_me = t*datap1 + tm1*data;
	reduce_me *= (_vrp1 - _vr);
	if(id >= (size-1)) reduce_me = 0.0;
	
																								LOG_ARRAY("J1 input", reduce_me, "\t"); 
	/* do the 'integral' (sum) note: this is reversed. */
	double final_J1 = warpScanInclusiveSimpleWarp_rev(reduce_me, s_Data, padded_size);
	//double final_J1 = scan1Inclusive_rev(reduce_me, s_Data, padded_size);
	
	/* save result in global (persistent) mem for next subkernal */
	cell_data__precomp5[pos] = final_J1;
	
}

__global__ void rk4__ee_v6_sortee_FAGGOT_split3__3(	double * __restrict__ fin, double * __restrict__ dest, const double * __restrict__ vr, 
								    double * __restrict__ cell_data__loglam, 
									double * __restrict__ cell_data__precomp1,
									double * __restrict__ cell_data__precomp2,
									double * __restrict__ cell_data__precomp3,
									double * __restrict__ cell_data__precomp4,
									double * __restrict__ cell_data__precomp5,
									uint size, int padded_size, int debug, double h, int numh) {
	
	__shared__ double logLambda;
	//__shared__ double s_Data[2 * THREADBLOCK_SIZE];
	
	/* ------------------------------------------------- */
	/* --- restart and resume place boilerplate--------- */
	uint pos = blockIdx.x* blockDim.x + threadIdx.x;
	
	double _vr = vr[id];									/* get momentum basis value for this cell. */
	//double _vrm1 = 0;										/* get vr for point one behind. */
	//if(id > 0) _vrm1 = vr[id-1];
	//double _vrp1 = 0;										/* get vr for point one ahead. */
	//if(id < size-1) _vrp1 = vr[id+1];
	
	double data = fin[pos];									/* get data point. */
	//double datap1 = 0.0;									/* get data point for one position ahead. */
	//if(id < size-1) datap1 = fin[pos+1];
	//double datam1 = 0.0;									/* get data point for one position behind. */
	//if(id > 0) datam1 = fin[pos-1];
	
	/* --- END restart and resume place boilerplate----- */
	/* ------------------------------------------------- */
	
	/*
	double t = .5*_vrp1;
	double tm1 = .5*_vr;
	double reduce_me = t*datap1 + tm1*data;
	
	
	reduce_me *= (_vrp1 - _vr);
	if(id >= (size-1)) reduce_me = 0.0;
																								LOG_ARRAY("J1 input", reduce_me, "\t"); 
	
	double final_J1 = scan1Inclusive_rev(reduce_me, s_Data, padded_size);
	__syncthreads();
																								LOG_ARRAY("J1[n] processed", final_J1, "\t");
	*/
	
	/* recall the values calculated by the previous sub kernals */
	double final_I4 = cell_data__precomp1[pos];
	double final_I2 = cell_data__precomp4[pos];
	double final_J1 = cell_data__precomp5[pos];
	
																								LOG_ARRAY("Data from Split-Step-1 coming in Split-Step-3: ", final_I4, "\t");
																								LOG_ARRAY("Data from Split-Step-2 coming in Split-Step-3: ", final_I2, "\t");
	
	if(id == size-1) {
        //double logLambda;
		logLambda = 2;
        double ne = 4.0*M_PI*final_I2; double Te = 4.0*M_PI*final_I4;
		// if the density is positive
        if (ne > 0.000000001) {
            Te /= (3.0*ne);
            Te *= 511000; // Temperature in eV
            ne *= density_np;

            Te = log(Te); 
            ne = log(ne);
            logLambda = 23.5 - 0.5*ne + 1.25*Te - sqrt(0.00001+0.0625*(Te-2.0)*(Te-2.0));

            if (logLambda < 2.0) logLambda = 2;
        }
		//logLambda = LOGee_cuda(4.0*M_PI*final_I2,4.0*M_PI*final_I4);
	}
																								LOG_ROOT("\tLog Lambda\n\t\t%e", 	logLambda);
	__syncthreads();
	
	
	final_I4 = data*final_I4 + _vr*_vr*_vr * final_J1 * data - 3.0*final_J1*final_I2;
	
																								LOG_ARRAY("G intital", final_I4, "\t");
	/* fixup the lower cells (filter them) to fix up issues iwth shpherical harmonincs */
	if (id < NB) { 
		final_I4 = G_cuda_v6(	id, fin[blockIdx.x*blockDim.x], 
								fin[blockIdx.x*blockDim.x + 1], 
								data, vr[0], vr[1], _vr, final_J1);
	}
																								LOG_ARRAY("G after lower boundry conditions applied", final_I4, "\t");
	
	final_I4 *= logLambda;
	
	/* save result in global (persistent) mem for next subkernal */
	cell_data__precomp1[pos] = final_I4;
																								LOG_ROOT("------------END OF SPLIT3_3--------------", debug);
	
}

__global__ void rk4__ee_v6_sortee_FAGGOT(	double * __restrict__ fin, 
											double * __restrict__ trial_input_vector_to_use_in_evalation, 
											const double * __restrict__ vr, 
											double * __restrict__ cell_data__loglam, 
											double * __restrict__ cell_data__precomp1,
											double * __restrict__ cell_data__precomp2,
											double * __restrict__ cell_data__precomp3,
											uint size, int padded_size, 
											int debug, double h, int numh) {	
	
	//__shared__ double s_Data[2 * THREADBLOCK_SIZE];
	__shared__ double s_Data[128];
	__shared__ double logLambda;
	// Note: when comparing values in this funtion with the CPU e-e collsion code, keep in mind these values will be different
	//  by a factor of LogLambda (in the cpu code LogLombda is multiplied in the last step.. in the GPU code, it is mulitped eailier in the process)
	uint pos = blockIdx.x* blockDim.x + threadIdx.x;
	
	double _vr = vr[id];
	double _vrm1 = 0;
	double _vrp1 = 0;
	if(id > 0) _vrm1 = vr[id-1];
	if(id < size-1) _vrp1 = vr[id+1];
	
	__syncthreads();
	//if(rk4_eval_step==1) {
		double data = fin[pos];
		double datam1 = 0.0;					// TODO: this needs to be fixed.I just wanna allocate the buffer to be a bit bigger do pos-1 won't seg fault at 0.
		if(id > 0) datam1 = fin[pos-1];			
		double datap1 = 0.0;
		if(id < size-1) datap1 = fin[pos+1];
	//} else {
	//	data = trial_input_vector_to_use_in_evalation[id];
	//	datam1 = 0.0;					// TODO: this needs to be fixed.I just wanna allocate the buffer to be a bit bigger do pos-1 won't seg fault at 0.
	//	if(id > 0) datam1 = trial_input_vector_to_use_in_evalation[id-1];			
	//	datap1 = 0.0;
	//	if(id < size-1) datap1 = trial_input_vector_to_use_in_evalation[id+1];		
	//}
																					LOG_ARRAY("Incoming Data to New Eval Step", data, "\t"); 
	double t = .5*_vr*_vr*_vr*_vr;
	double tm1 = .5*_vrm1*_vrm1*_vrm1*_vrm1;//*(_vr - _vrm1);
	double reduce_me = t*data + tm1*datam1;
	reduce_me *= (_vr - _vrm1);
	if(id == 0) reduce_me = 0.0;
	
	//double final_I4 = scan1Inclusive(reduce_me, s_Data, padded_size);
	//double final_I4 = warpScanInclusiveSimpleWarp(reduce_me, s_Data, padded_size);
	double final_I4 = warpScanInclusiveSimpleWarp(reduce_me, s_Data, padded_size);
	
	
																					LOG_ARRAY("I4[n] inital", final_I4, "\t"); 
	t = .5*_vr*_vr;
	tm1 = .5*_vrm1*_vrm1;
	reduce_me = t*data + tm1*datam1;
	reduce_me *= (_vr - _vrm1); 
	if(id == 0) reduce_me = 0.0;
	__syncthreads(); // suspect.
	
	//double final_I2 = scan1Inclusive(reduce_me, s_Data, padded_size);
	double final_I2 = warpScanInclusiveSimpleWarp(reduce_me, s_Data, padded_size);
	
																					LOG_ARRAY("I2[n] inital", final_I2, "\t"); 

	if(id == size-1) {
		logLambda = LOGee_cuda(4.0*M_PI*final_I2,4.0*M_PI*final_I4);//*c_kpre;	
		//printf("\tLog lambda (args=%e,%e): %e", 4.0*M_PI*final_I2, 4.0*M_PI*final_I4, logLambda);
		//logLambda *= c_kpre;
	}

	//if(DBROOT) {
	//	printf("\npadding %e %e %e %e \n", fin[blockIdx.x* blockDim.x + 28],fin[blockIdx.x* blockDim.x + 29],fin[blockIdx.x* blockDim.x + 30],fin[blockIdx.x* blockDim.x + 31]);
	//}
	
	t = .5*_vrp1;
	tm1 = .5*_vr;
	reduce_me = t*datap1 + tm1*data;
	
	
	reduce_me *= (_vrp1 - _vr);
	//if(id == (size-1)) reduce_me = 0.0;
	// to be very safe for now.. relax later...
	if(id >= (size-1)) reduce_me = 0.0;
																					LOG_ARRAY("J1 input", reduce_me, "\t"); 


	//double final_J11 = scan1Inclusive_rev(reduce_me, s_Data, padded_size);
	double final_J1 = warpScanInclusiveSimpleWarp_rev(reduce_me, s_Data, padded_size);
																					LOG_ARRAY("J1[n] processed", final_J1, "\t");
/* 
	//
	// way 1 stop computation here ( Used 31 registers, 4096+0 bytes smem) these are nice numbers.
	//

	cell_data__precomp1[blockIdx.x] = final_I4;
	cell_data__precomp2[blockIdx.x] = final_I2;
	cell_data__precomp3[blockIdx.x] = final_J1;
*/
	//if( id  == (size-1)) {
	//	double logglam = LOGee_cuda(4.0*M_PI*final_I2, 4.0*M_PI*final_I4);				if(DBROOT) printf("\n\tLoglam\n\t\t%e", (*logglam));
	//}

	//
	// Also try doing a bit more computation over here..
	//
	final_I4 = data*final_I4 + _vr*_vr*_vr * final_J1 * data - 3.0*final_J1*final_I2;		LOG_ARRAY("G intital", final_I4, "\t");
																							
	
	if (id < NB) {
		final_I4 = G_cuda_v6(id, fin[blockIdx.x*blockDim.x], fin[blockIdx.x*blockDim.x + 1], data, vr[0], vr[1], _vr, final_J1);
	}
																							LOG_ROOT("\tLog Lambda\n\t\t%e", 	logLambda);																			

	///     Evaluate G assuming a parabolic f(v << vt)
	//for (int n(0); n < NB; ++n) { 
	//	I4[n] = G(n,fin);
	//}
																							LOG_ARRAY("G after lower boundry conditions applied", final_I4, "\t");
																							
	// this may not work... but pre normalize...
	
    final_I4 *= logLambda;


	//cell_data__precomp1[blockIdx.x] = final_I4;
	cell_data__precomp1[pos] = final_I4;
	__syncthreads();
}

/*
//    split eval into 2 functions!
//
*/
__global__ void rk4__ee_v6_sortee_FAGGOT_split2__1__v2(	double * __restrict__ fin,
														double * __restrict__ dest,
														const double * __restrict__ vr, 
														double * __restrict__ cell_data__loglam, 
														double * __restrict__ cell_data__precomp1,
														double * __restrict__ cell_data__precomp2,
														double * __restrict__ cell_data__precomp3,
														uint size, int padded_size, int debug, double h, int numh) {
	
	__shared__ double s_Data[2 * THREADBLOCK_SIZE];
	
	/* ------------------------------------------------- */
	/* --- restart and resume place boilerplate--------- */
	uint pos = blockIdx.x* blockDim.x + threadIdx.x;
	
	double _vr = vr[id];										/* get momentum basis value for this cell. */
	double _vrm1 = 0;											/* get vr for point one behind. */
	if(id > 0) _vrm1 = vr[id-1];
	//double _vrp1 = 0;											/* get vr for point one ahead. */
	//if(id < size-1) _vrp1 = vr[id+1];
	
	double data = fin[pos];										/* get data point. */
	double datam1 = 0.0;										/* get data point for one position behind. */
	if(id > 0) datam1 = fin[pos-1];
	//double datap1 = 0.0;										/* get data point for one position ahead. */
	//if(id < size-1) datap1 = fin[pos+1];
	
	/* --- END restart and resume place boilerplate----- */
	/* ------------------------------------------------- */
	
																								LOG_ARRAY("Incoming Data to New Eval Step", data, "\t"); 
	double t = .5*_vr*_vr*_vr*_vr;
	double tm1 = .5*_vrm1*_vrm1*_vrm1*_vrm1;//*(_vr - _vrm1);
	double reduce_me = t*data + tm1*datam1;
	reduce_me *= (_vr - _vrm1);
	if(id == 0) reduce_me = 0.0;
	
	double final_I4 = scan1Inclusive(reduce_me, s_Data, padded_size);
																								LOG_ARRAY("I4[n] inital", final_I4, "\t"); 
	
	/* save result in global (persistent) mem for next subkernal */
	cell_data__precomp1[pos] = final_I4;
																								LOG_ROOT("------------END OF SPLIT1--------------", debug);
	
}

__global__ void rk4__ee_v6_sortee_FAGGOT_split2__2__v2(	double * __restrict__ fin, 
														double * __restrict__ dest,
														const double * __restrict__ vr, 
														double * __restrict__ cell_data__loglam, 
														double * __restrict__ cell_data__precomp1,
														double * __restrict__ cell_data__precomp2,
														double * __restrict__ cell_data__precomp3,
														uint size, int padded_size, int debug, double h, int numh) {	
	// --- restart and resume place boilerplate---------
	__shared__ double s_Data[2 * THREADBLOCK_SIZE];
	__shared__ double logLambda;
	
	/* ------------------------------------------------- */
	/* --- restart and resume place boilerplate--------- */
	uint pos = blockIdx.x* blockDim.x + threadIdx.x;
	
	double _vr = vr[id];										/* get momentum basis value for this cell. */
	double _vrm1 = 0;											/* get vr for point one behind. */
	if(id > 0) _vrm1 = vr[id-1];
	double _vrp1 = 0;											/* get vr for point one ahead. */
	if(id < size-1) _vrp1 = vr[id+1];
	
	double data = fin[pos];										/* get data point. */
	double datam1 = 0.0;										/* get data point for one position behind. */
	if(id > 0) datam1 = fin[pos-1];
	double datap1 = 0.0;										/* get data point for one position ahead. */
	if(id < size-1) datap1 = fin[pos+1];
	
	/* --- END restart and resume place boilerplate----- */
	/* ------------------------------------------------- */
	
	/* load previously saved result */
	double final_I4 = cell_data__precomp1[pos];
																								LOG_ARRAY("Data from Split-Step-1 coming in Split-Step-2: ", final_I4, "\t");
	
	
	/* now set up the calculation */
	double t = .5*_vr*_vr;
	double  tm1 = .5*_vrm1*_vrm1;
	double  reduce_me = t*data + tm1*datam1;
	reduce_me *= (_vr - _vrm1); 
	if(id == 0) reduce_me = 0.0;
	
	/* do the 'integral' (sum) */
	double final_I2 = scan1Inclusive(reduce_me, s_Data, padded_size);
																								LOG_ARRAY("I2[n] inital", final_I2, "\t"); 
	
	if(id == size-1) {
		logLambda = LOGee_cuda(4.0*M_PI*final_I2,4.0*M_PI*final_I4);
	}
	
	
	/* now set up the calculation */
	t = .5*_vrp1;
	tm1 = .5*_vr;
	reduce_me = t*datap1 + tm1*data;
	reduce_me *= (_vrp1 - _vr);
	// to be very safe for now.. relax later...
	if(id >= (size-1)) reduce_me = 0.0;
																								LOG_ARRAY("J1 input", reduce_me, "\t"); 
	/* do the 'integral' (sum) note: this is reversed. */
	double final_J1 = scan1Inclusive_rev(reduce_me, s_Data, padded_size);
	__syncthreads();
																								LOG_ARRAY("J1[n] calculated", final_J1, "\t");
	
	
	/* calculate the result we after by combining all the previous integral values! */
	final_I4 = data*final_I4 + _vr*_vr*_vr * final_J1 * data - 3.0*final_J1*final_I2;
																								LOG_ARRAY("G intital", final_I4, "\t");
	
	/* fixup the lower cells (filter them) to fix up issues iwth shpherical harmonincs */
	if (id < NB) {
		final_I4 = G_cuda_v6(	id, fin[blockIdx.x*blockDim.x],
								fin[blockIdx.x*blockDim.x + 1], 
								data, vr[0], vr[1], _vr, final_J1);
	}
	
																								LOG_ROOT("\tLog Lambda\n\t\t%e", 	logLambda);																			
																								LOG_ARRAY("G after lower boundry conditions applied", final_I4, "\t");
	final_I4 *= logLambda;
	
	/* save result in global (persistent) mem for next subkernal */
	cell_data__precomp1[pos] = final_I4;
																								LOG_ROOT("----------------END OF SPLIT2---------------", debug);
}

void test_scans(double * trash_buffer) {
	
	dim3 dimGrid_test(1);
	dim3 dimBlock_test(128);
	
	/*
	{ ScopeTimer ___t("scan_test_0");
		for(int i = 0; i < 1000;i++) {
			test_forward_prefix_scans<0><<<dimGrid_test, dimBlock_test>>>(trash_buffer);
		}
	}
	
	{ ScopeTimer ___t("scan_test_1");
		for(int i = 0; i < 1000;i++) {
			test_forward_prefix_scans<1><<<dimGrid_test, dimBlock_test>>>(trash_buffer);
		}
	}
	
	{ ScopeTimer ___t("scan_test_2");
		for(int i = 0; i < 1000;i++) {
			test_forward_prefix_scans<3><<<dimGrid_test, dimBlock_test>>>(trash_buffer);
		}
	}
	
	
	double timer0 = TimerOracle::summon()->get_result("scan_test_0")->dt;
	double timer1 = TimerOracle::summon()->get_result("scan_test_1")->dt;
	double timer2 = TimerOracle::summon()->get_result("scan_test_2")->dt;
	
	cout << "Method 0: " << timer0 << " ms" << endl;
	cout << "Method 1: " << timer1 << " ms" << endl;
	cout << "Method 2: " << timer2 << " ms" << endl;
	*/
	
	test_forward_prefix_scans<2><<<dimGrid_test, dimBlock_test>>>(trash_buffer);
}

// Good entry point... works..
void eval_rk4_v6 (	double * fc, double * dest, double *vr, 
							//double * U1, double *U1m1, double * U2, double *U2m1,
							//double * U3,
							//double * U4, double *U4m1,
							//double * Pn, double *Qn,
							double * __restrict__ cell_data__precomp1,   // these should be the size of the number of cells.*pr
							double * __restrict__ cell_data__precomp2,
							double * __restrict__ cell_data__precomp3,
							double * __restrict__ cell_data__precomp4,
							double * __restrict__ cell_data__precomp5,
							uint size, int padded_size, int debug, double h, int numh, int num_cells_x, int num_cells_y,
																const double * __restrict__ U4,
									const double * __restrict__ U4m1,
									const double * __restrict__ U2,
									const double * __restrict__ U2m1) {
	
	// I gotta clean up threse names... the ultimate answer is passed back to Host via the first argument (fc)
	//									this buffer islso the buffer that has the incoming data.
	//
	//									for this reason, it is not altered until the last part of an RK4 step.
	//
	
	// cheat for now and bump up mem size so eveything is even.
	int blocks_needed = size / (1*THREADBLOCK_SIZE);
	if((size % (1*THREADBLOCK_SIZE)) != 0) {
		blocks_needed += 1;
	}
	//
	
	// this is linear 1-d grid array.
	dim3 dimBlock(padded_size);
	dim3 dimGrid(num_cells_x*num_cells_y);
	
	double * yn = cell_data__precomp3;
	double * y_new = cell_data__precomp2;
	//numh = 1;  // for debugging.. so this loop runs only once.
	
	//cout << "Hi from eval test!" << padded_size << endl;
	
	for(int i =0;i<numh;i++) {
		if(false) { 
			//cudaStream_t stream[3]; 
			//for (int i = 0; i < 3; ++i) 
			//	cudaStreamCreate(&stream[i]);
			//cout << "Trying super streams! " << endl << blocks_needed << ":"<< padded_size << endl;
			if(true) {
				// 
				// first evaulation of the function (1/4th of an RK4 step) should read from (have 1st argument) fc (current function state)
				// after 1st evaluation, the next state to use is copied into the 'yn' variable.. so all other evals should be reading from this.
				// ynew hold the answer as it is axxumated fron eaxhg step (it's a sum)
				//
				// this stramge stucture exists because RK4 does 4  funtion evaulations.. at 'locations' that are relative to
				// the input... so the input (fc) needsd to exist unchanged (can't be overwritten) until very last part of an RK4 step.
				// so the anser must be accumated in another buffer and the the trial evaluation point in a third.
				//
		// RK4 eval 1
				rk4__ee_v6_sortee_FAGGOT_split3__1<<< dimGrid, dimBlock>>>(	fc, yn, vr,  //fc, dest, vr, 
														cell_data__precomp1,
														cell_data__precomp1,
														cell_data__precomp2,
														cell_data__precomp3,
														cell_data__precomp4,
														size, padded_size, 1, h, numh, U4, U4m1, U2, U2m1);
				
				rk4__ee_v6_sortee_FAGGOT_split3__2<<< dimGrid, dimBlock>>>(	fc, yn, vr,  //fc, dest, vr, 
														cell_data__precomp1,
														cell_data__precomp1,
														cell_data__precomp2,
														cell_data__precomp3,
														cell_data__precomp4,
														size, padded_size, 1, h, numh);
				rk4__ee_v6_sortee_FAGGOT_split3__2_5<<< dimGrid, dimBlock>>>(	fc, yn, vr,  //fc, dest, vr, 
														cell_data__precomp1,
														cell_data__precomp1,
														cell_data__precomp2,
														cell_data__precomp3,
														cell_data__precomp4,
														cell_data__precomp5,
														size, padded_size, 1, h, numh);
				rk4__ee_v6_sortee_FAGGOT_split3__3<<< dimGrid, dimBlock>>>(	fc, yn, vr,  //fc, dest, vr, 
														cell_data__precomp1,
														cell_data__precomp1,
														cell_data__precomp2,
														cell_data__precomp3,
														cell_data__precomp4,
														cell_data__precomp5,
														size, padded_size, 1, h, numh);
				//cudaDeviceSynchronize();
				rk4__ee_v6_sortee_PHYSCO<1><<< dimGrid, dimBlock >>>(	fc, yn, y_new, vr, cell_data__precomp1, 
												cell_data__precomp1,
												cell_data__precomp2,
												cell_data__precomp3,
												size, padded_size, 1, h, numh, .5*h, h/6.0);
				
				// notice that the frist 2 args of 'faggor' are reverese relative to the previous invokations.
		// RK4 eval 2
				rk4__ee_v6_sortee_FAGGOT_split3__1<<< dimGrid, dimBlock>>>( yn, fc, vr, //fc, dest, vr, 
														cell_data__precomp1,
														cell_data__precomp1,
														cell_data__precomp2,
														cell_data__precomp3,
														cell_data__precomp4,
														size, padded_size, 1, h, numh, U4, U4m1, U2, U2m1);
				
				rk4__ee_v6_sortee_FAGGOT_split3__2<<< dimGrid, dimBlock>>>(	yn, fc, vr, //fc, dest, vr, 
														cell_data__precomp1,
														cell_data__precomp1,
														cell_data__precomp2,
														cell_data__precomp3,
														cell_data__precomp4,
														size, padded_size, 2, h, numh);
				rk4__ee_v6_sortee_FAGGOT_split3__2_5<<< dimGrid, dimBlock>>>(	fc, yn, vr,  //fc, dest, vr, 
														cell_data__precomp1,
														cell_data__precomp1,
														cell_data__precomp2,
														cell_data__precomp3,
														cell_data__precomp4,
														cell_data__precomp5,
														size, padded_size, 2, h, numh);
				rk4__ee_v6_sortee_FAGGOT_split3__3<<< dimGrid, dimBlock>>>(	yn, fc, vr, //fc, dest, vr, 
														cell_data__precomp1,
														cell_data__precomp1,
														cell_data__precomp2,
														cell_data__precomp3,
														cell_data__precomp4,
														cell_data__precomp5,
														size, padded_size, 2, h, numh);
				//cudaDeviceSynchronize();
				rk4__ee_v6_sortee_PHYSCO<2><<< dimGrid, dimBlock >>>(	fc, yn, y_new, vr, cell_data__precomp1, 
												cell_data__precomp1,
												cell_data__precomp2,
												cell_data__precomp3,
												size, padded_size, 2, h, numh, .5*h, h/3.0);
				
				
		// RK4 eval 3
				rk4__ee_v6_sortee_FAGGOT_split3__1<<< dimGrid, dimBlock>>>(	yn, fc, vr, //fc, dest, vr, 
														cell_data__precomp1,
														cell_data__precomp1,
														cell_data__precomp2,
														cell_data__precomp3,
														cell_data__precomp4,
														size, padded_size, 1, h, numh, U4, U4m1, U2, U2m1);
				
				rk4__ee_v6_sortee_FAGGOT_split3__2<<< dimGrid, dimBlock>>>(	yn, fc, vr, //fc, dest, vr, 
														cell_data__precomp1,
														cell_data__precomp1,
														cell_data__precomp2,
														cell_data__precomp3,
														cell_data__precomp4,
														size, padded_size, 3, h, numh);
				rk4__ee_v6_sortee_FAGGOT_split3__2_5<<< dimGrid, dimBlock>>>(	fc, yn, vr,  //fc, dest, vr, 
														cell_data__precomp1,
														cell_data__precomp1,
														cell_data__precomp2,
														cell_data__precomp3,
														cell_data__precomp4,
														cell_data__precomp5,
														size, padded_size, 3, h, numh);
				rk4__ee_v6_sortee_FAGGOT_split3__3<<< dimGrid, dimBlock>>>(	yn, fc, vr, //fc, dest, vr, 
														cell_data__precomp1,
														cell_data__precomp1,
														cell_data__precomp2,
														cell_data__precomp3,
														cell_data__precomp4,
														cell_data__precomp5,
														size, padded_size, 3, h, numh);
				//cudaDeviceSynchronize();
				rk4__ee_v6_sortee_PHYSCO<3><<< dimGrid, dimBlock >>>(	fc, yn, y_new, vr, cell_data__precomp1, 
												cell_data__precomp1,
												cell_data__precomp2,
												cell_data__precomp3,
												size, padded_size, 3, h, numh, h, h/3.0);
				
				
		// RK4 eval 4
				rk4__ee_v6_sortee_FAGGOT_split3__1<<< dimGrid, dimBlock>>>(	yn, fc, vr, //fc, dest, vr, 
														cell_data__precomp1,
														cell_data__precomp1,
														cell_data__precomp2,
														cell_data__precomp3,
														cell_data__precomp4,
														size, padded_size, 1, h, numh, U4, U4m1, U2, U2m1);
				
				rk4__ee_v6_sortee_FAGGOT_split3__2<<< dimGrid, dimBlock>>>(	yn, fc, vr, //fc, dest, vr, 
														cell_data__precomp1,
														cell_data__precomp1,
														cell_data__precomp2,
														cell_data__precomp3,
														cell_data__precomp4,
														size, padded_size, 4, h, numh);
				rk4__ee_v6_sortee_FAGGOT_split3__2_5<<< dimGrid, dimBlock>>>(	fc, yn, vr,  //fc, dest, vr, 
														cell_data__precomp1,
														cell_data__precomp1,
														cell_data__precomp2,
														cell_data__precomp3,
														cell_data__precomp4,
														cell_data__precomp5,
														size, padded_size, 4, h, numh);
				rk4__ee_v6_sortee_FAGGOT_split3__3<<< dimGrid, dimBlock>>>(	yn, fc, vr, //fc, dest, vr, 
														cell_data__precomp1,
														cell_data__precomp1,
														cell_data__precomp2,
														cell_data__precomp3,
														cell_data__precomp4,
														cell_data__precomp5,
														size, padded_size, 4, h, numh);
				//cudaDeviceSynchronize();
				rk4__ee_v6_sortee_PHYSCO<4><<< dimGrid, dimBlock >>>(	fc, yn, y_new, vr, cell_data__precomp1, 
												cell_data__precomp1,
												cell_data__precomp2,
												cell_data__precomp3,
												size, padded_size, 4, h, numh, h, h/6.0);
				
				
			} else {
				//cout << "Using split 2 version" << endl;
				
				// aternate version with 2 split threads...
		// RK4 eval 1
				rk4__ee_v6_sortee_FAGGOT_split2__1__v2<<< dimGrid, dimBlock>>>(	fc, yn, vr,  //fc, dest, vr, 
														cell_data__precomp1,
														cell_data__precomp1,
														cell_data__precomp2,
														cell_data__precomp3,
														size, padded_size, 1, h, numh);
				
				rk4__ee_v6_sortee_FAGGOT_split2__2__v2<<< dimGrid, dimBlock>>>(	fc, yn, vr, 
														cell_data__precomp1,
														cell_data__precomp1,
														cell_data__precomp2,
														cell_data__precomp3,
														size, padded_size, 1, h, numh);
				
				//cudaDeviceSynchronize();
				
				rk4__ee_v6_sortee_PHYSCO<1><<< dimGrid, dimBlock >>>(	fc, yn, y_new, vr, cell_data__precomp1, 
												cell_data__precomp1,
												cell_data__precomp2,
												cell_data__precomp3,
												size, padded_size, 1, h, numh, .5*h, h/6.0);
				
		// RK4 eval 2
				rk4__ee_v6_sortee_FAGGOT_split2__1__v2<<< dimGrid, dimBlock>>>(	yn, fc, vr,
														cell_data__precomp1,
														cell_data__precomp1,
														cell_data__precomp2,
														cell_data__precomp3,
														size, padded_size, 2, h, numh);
				
				rk4__ee_v6_sortee_FAGGOT_split2__2__v2<<< dimGrid, dimBlock>>>(	yn, fc, vr, 
														cell_data__precomp1,
														cell_data__precomp1,
														cell_data__precomp2,
														cell_data__precomp3,
														size, padded_size, 2, h, numh);
				
				//cudaDeviceSynchronize();
				
				rk4__ee_v6_sortee_PHYSCO<2><<< dimGrid, dimBlock >>>(	fc, yn, y_new, vr, cell_data__precomp1, 
												cell_data__precomp1,
												cell_data__precomp2,
												cell_data__precomp3,
												size, padded_size, 2, h, numh, .5*h, h/3.0);
		// RK4 eval 3
				rk4__ee_v6_sortee_FAGGOT_split2__1__v2<<< dimGrid, dimBlock>>>(	yn, fc, vr,
														cell_data__precomp1,
														cell_data__precomp1,
														cell_data__precomp2,
														cell_data__precomp3,
														size, padded_size, 3, h, numh);
				
				rk4__ee_v6_sortee_FAGGOT_split2__2__v2<<< dimGrid, dimBlock>>>(	yn, fc, vr, 
														cell_data__precomp1,
														cell_data__precomp1,
														cell_data__precomp2,
														cell_data__precomp3,
														size, padded_size, 3, h, numh);
				
				//cudaDeviceSynchronize();
				
				rk4__ee_v6_sortee_PHYSCO<3><<< dimGrid, dimBlock >>>(	fc, yn, y_new, vr, cell_data__precomp1, 
												cell_data__precomp1,
												cell_data__precomp2,
												cell_data__precomp3,
												size, padded_size, 3, h, numh, h, h/3.0);
		// RK4 eval 4
				rk4__ee_v6_sortee_FAGGOT_split2__1__v2<<< dimGrid, dimBlock>>>(	yn, fc, vr,
														cell_data__precomp1,
														cell_data__precomp1,
														cell_data__precomp2,
														cell_data__precomp3,
														size, padded_size, 4, h, numh);
				
				rk4__ee_v6_sortee_FAGGOT_split2__2__v2<<< dimGrid, dimBlock>>>(	yn, fc, vr,
														cell_data__precomp1,
				 										cell_data__precomp1,
														cell_data__precomp2,
														cell_data__precomp3,
														size, padded_size, 4, h, numh);
				
				//cudaDeviceSynchronize();
				
				rk4__ee_v6_sortee_PHYSCO<4><<< dimGrid, dimBlock >>>(	fc, yn, y_new, vr, cell_data__precomp1, 
												cell_data__precomp1,
												cell_data__precomp2,
												cell_data__precomp3,
												size, padded_size, 4, h, numh, h, h/6.0); 
				
			}
			//cudaStreamDestroy(stream[0]);
			//cudaStreamDestroy(stream[1]);
			//cudaStreamDestroy(stream[2]);
		} else {
		// cudaDeviceSynchronize();	
			
																												{ //ScopeTimer ___t("CUDA e*e collsion kernal launch time per RK4 advance");
		// RK4 eval 1
			rk4__ee_v6_sortee_FAGGOT	<<< dimGrid, dimBlock >>>(	fc, yn, vr, 
																	cell_data__precomp1,
																	cell_data__precomp1,
																	cell_data__precomp2,
																	cell_data__precomp3,
																	size, padded_size, 1, h, numh); 
			
			rk4__ee_v6_sortee_PHYSCO<1> <<< dimGrid, dimBlock >>>(	fc, yn, y_new, vr, 
																	cell_data__precomp1, 
																	cell_data__precomp1,
																	cell_data__precomp2,
																	cell_data__precomp3,
																	size, padded_size, 1, h, numh, 0.5*h, h/6.0); 
			
		// RK4 eval 2
			rk4__ee_v6_sortee_FAGGOT	<<< dimGrid, dimBlock >>>(	yn, fc, vr, 		// Note the reversed order of the 1st 2 arguments fc, yn
																	cell_data__precomp1,
																	cell_data__precomp1,
																	cell_data__precomp2,
																	cell_data__precomp3,
																	size, padded_size, 2, h, numh); 
			
			rk4__ee_v6_sortee_PHYSCO<2> <<< dimGrid, dimBlock >>>(	fc, yn, y_new, vr, 
																	cell_data__precomp1, 
																	cell_data__precomp1,
																	cell_data__precomp2,
																	cell_data__precomp3,
																	size, padded_size, 2, h, numh, .5*h, h/3.0); 
																					// 2
			rk4__ee_v6_sortee_FAGGOT	<<< dimGrid, dimBlock >>>(	yn, fc, vr, 		// Note the reversed order of the 1st 2 arguments fc, yn
																	cell_data__precomp1,
																	cell_data__precomp1,
																	cell_data__precomp2,
																	cell_data__precomp3,
																	size, padded_size, 3, h, numh); 
			
		// RK4 eval 3
			rk4__ee_v6_sortee_PHYSCO<3>	<<< dimGrid, dimBlock >>>(	fc, yn, y_new, vr, 
																	cell_data__precomp1, 
																	cell_data__precomp1,
																	cell_data__precomp2,
																	cell_data__precomp3,
																	size, padded_size, 3, h, numh, h, h/3.0);
																					// 3
			rk4__ee_v6_sortee_FAGGOT<<< dimGrid, dimBlock >>>(		yn, fc, vr, 		// Note the reversed order of the 1st 2 arguments fc, yn
																	cell_data__precomp1,
																	cell_data__precomp1,
																	cell_data__precomp2,
																	cell_data__precomp3,
																	size, padded_size, 4, h, numh); 
			
		// RK4 eval 4
			rk4__ee_v6_sortee_PHYSCO<4><<< dimGrid, dimBlock >>>(	fc, yn, y_new, vr, 
																	cell_data__precomp1, 
																	cell_data__precomp1,
																	cell_data__precomp2,
																	cell_data__precomp3,
  																	size, padded_size, 4, h, numh, h, h/6.0); 
																					// 4
																													}
		}
	}
	
	cutilCheckMsg("An error occured while launching the kernal. Exiting.");
	//cout << "\tLaunched " << (numh*8) << " kernals for a " << numh << " consecutive RK4 step process.\n";
}
