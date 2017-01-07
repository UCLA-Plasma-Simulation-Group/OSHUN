// The following ifdef block is the standard way of creating macros which make exporting 
// from a DLL simpler. All files within this DLL are compiled with the NUMPY1_EXPORTS
// symbol defined on the command line. this symbol should not be defined on any project
// that uses this DLL. This way any other project whose source files include this file see 
// NUMPY1_API functions as being imported from a DLL, whereas this DLL sees symbols
// defined with this macro as being exported.
#ifdef NUMPY1_EXPORTS
#define NUMPY1_API __declspec(dllexport)
#else
#define NUMPY1_API 
#endif


extern NUMPY1_API int nnumpy1;

NUMPY1_API int fnnumpy1(void);
