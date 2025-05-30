#include "pch.h"
#include <vector>

#ifdef DIFFURSOLVERDLL_EXPORTS
#define DIFFURSOLVERDLL_API __declspec(dllexport)
#else
#define DIFFURSOLVERDLL_API __declspec(dllimport)
#endif

#define ui unsigned int
#define vecF std::vector<F1System>
#define vecRealF std::vector<F>


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

    DIFFURSOLVERDLL_API int  SolveDiffUrAutoHArr(double x0, double y0, double x1, double h, double eps, F1 f, std::vector<double>& resX, std::vector<double>& resY, Methods method);

    DIFFURSOLVERDLL_API int  SolveDiffUrArr(double x0, double y0, double x1, double h, F1 f, std::vector< double>& res, Methods method);

    DIFFURSOLVERDLL_API int  SolveDiffUrSystem(double x0, double* y0, double x1, ui n, double h, vecF f, double* res, Methods method);

    DIFFURSOLVERDLL_API int  SolveDiffUrSystemArr(double x0, double* y0, double x1, ui n, double h, vecF f, double** res, Methods method);

}