#include "pch.h"
#include <vector>

#ifdef DIFFURSOLVERDLL_EXPORTS
#define DIFFURSOLVERDLL_API __declspec(dllexport)
#else
#define DIFFURSOLVERDLL_API __declspec(dllimport)
#endif

#define ui unsigned int


enum Methods
{
    Euler,
    RK2,
    RK4
};


extern "C" {
    
    typedef double(__cdecl* F1)(double x, double y);

    typedef double(__cdecl* F1System)(double x, double* y);

    typedef double(__cdecl* StepF)(double x0, double y0, double h, F1 f);

    typedef double(__cdecl* F)(double x);

    DIFFURSOLVERDLL_API int SolveDiffUr(double x0, double y0, double x1, double h, F1 f, double& res, Methods method);

    DIFFURSOLVERDLL_API int  SolveDiffUrAutoH(double x0, double y0, double x1, double h, double eps, F1 f, double& res, Methods method);

    DIFFURSOLVERDLL_API int  SolveDiffUrAutoHArr(double x0, double y0, double x1, double h, double eps, F1 f, double* resX, double* resY,int& resSize, Methods method);

    DIFFURSOLVERDLL_API int  SolveDiffUrAutoMethod(double x0, double y0, double x1, double h, double eps, F1 f, double& res, Methods method);

    DIFFURSOLVERDLL_API int  SolveDiffUrAutoMethodArr(double x0, double y0, double x1, double h, double eps, F1 f, double* res,int& sizeRes, Methods method);

    DIFFURSOLVERDLL_API int  SolveDiffUrAutoMethodArr2(double x0, double y0, double x1, double h, double eps, F1 f, double* res, int& sizeRes, Methods* methodsArr, Methods method);


    DIFFURSOLVERDLL_API int  SolveDiffUrArr(double x0, double y0, double x1, double h, F1 f, double* res, int& sizeRes, Methods method);

    DIFFURSOLVERDLL_API int  SolveDiffUrSystem(double x0, double* y0, double x1, ui n, double h, F1System* f, double* res, Methods method);

    DIFFURSOLVERDLL_API int  SolveDiffUrSystemArr(double x0, double* y0, double x1, ui n, double h, F1System* f, double** res, Methods method);

}