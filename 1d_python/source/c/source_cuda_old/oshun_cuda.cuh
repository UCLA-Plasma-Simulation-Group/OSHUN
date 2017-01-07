#ifndef DECLERATION_EXPLICIT_EC_FP_CUDA_CUH
#define DECLERATION_EXPLICIT_EC_FP_CUDA_CUH

	//#define FRED 1
	#ifdef FRED
		// These are here because we cannot have C Macros inside other Macros (nesting)
		//   so these provide a quick way to change which functions get debug output.
		inline __device__ bool match_bl()		{return (blockIdx.x==0);	}
		inline __device__ bool match_root()		{return (threadIdx.x==0);	}
		inline __device__ bool match_th_less()	{return (threadIdx.x < 5);	}
		inline __device__ bool match_db_flag1(int flag1) {return true; } //return (flag1==1);	}

		// Quick functions to control logging for denugging... These need to be Macros because We want to to able to
		//		compile them out of the code when not in use (especially because the printf function is an expensive, slow, 
		//		and dispurtive function to call on the GPU... we really don't want any to be compiled in production code.
		#define DB blockIdx.x == 0 && threadIdx.x<5 && debug > -1//== 1
		#define DBROOT blockIdx.x == 0 && threadIdx.x==0 && debug>-1// == 1

		// printf("\nvar " #VAR ": " #MESSAGE "\n");				/
		#define LOG_(X)													if(DB) { printf(X); }

		//if(blockIdx.x == 6 &&                && debug == 1) {  printf(X); }
		#define LOG_ROOT_TEXT(X)										if(match_bl()&&match_root()&&match_db_flag1(debug)) {  printf(X); }
		#define LOG_ROOT(FORMAT_STRING, VAL)							if(match_bl()&&match_root()&&match_db_flag1(debug)) {  printf("\n" FORMAT_STRING, VAL); }
		#define LOG_ARRAY(MESSAGE, VAR, INDENT)							if(match_bl()&&match_root()&&match_db_flag1(debug)) { printf("\n" INDENT MESSAGE "(var " #VAR "):" "\n" INDENT "\t"); }			\
																		if(match_bl()&&match_th_less()&&match_db_flag1(debug)) {																		\
																			printf("(%d:%d) %.9e ", blockIdx.x, threadIdx.x, VAR);																			\
																		}
		#define LOG_ARRAY_INT(MESSAGE, VAR, INDENT)						if(match_bl()&&match_root()&&match_db_flag1(debug)) { printf("\n" INDENT MESSAGE "(var " #VAR "):" "\n" INDENT "\t"); }			\
																		if(match_bl()&&match_th_less()&&match_db_flag1(debug)) {																		\
																			printf("(%d:%d) %d ", blockIdx.x, threadIdx.x, VAR);																			\
																		}
		#define LOG_ARRAY1(MESSAGE, MESSAGE_VAR, VAR, INDENT)			if(match_bl()&&match_root()&&match_db_flag1(debug)) { printf("\n" INDENT MESSAGE "(" #MESSAGE_VAR "=%e)" "\n" INDENT "\t", MESSAGE_VAR ); }		\
																		if(match_bl()&&match_th_less()&&match_db_flag1(debug)) {																						\
																			printf("(%d:%d) %.9e ", blockIdx.x, threadIdx.x, VAR);																							\
																		}
		#define LOG_ARRAY1_INT(MESSAGE, MESSAGE_VAR, VAR,INDENT)		if(match_bl()&&match_root()&&match_db_flag1(debug)) { printf("\n" INDENT MESSAGE "(" #MESSAGE_VAR "=%d)" "\n", INDENT "\t", MESSAGE_VAR ); }	\
																		if(match_bl()&&match_th_less()&&match_db_flag1(debug)) {																						\
																			printf("(%d:%d) %d ", blockIdx.x, threadIdx.x, VAR);																						\
																		}
	#else

		// or can use  (void)sizeof
		#define DB														(void)0
		//#define DBROOT													(void)0
		#define LOG_(X)													(void)0
		#define LOG_ROOT_TEXT(X)										(void)0
		#define LOG_ROOT(FORMAT_STRING, VAL)							(void)0
		#define LOG_ARRAY(MESSAGE, VAR, INDENT)							(void)0
		#define LOG_ARRAY_INT(MESSAGE, VAR, INDENT)						(void)0
		#define LOG_ARRAY1(MESSAGE, MESSAGE_VAR, VAR, INDENT)			(void)0
		#define LOG_ARRAY1_INT(MESSAGE, MESSAGE_VAR, VAR,INDENT)		(void)0

	#endif


#endif 