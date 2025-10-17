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

    typedef double(__cdecl* F1System)(double x,const double* y);

    typedef double(__cdecl* StepF)(double x0, double y0, double h, F1 f);

    typedef double(__cdecl* F)(double x);

    struct Task {
        double x0;
        double x1;
        double y0;
        ui n;
        double h;
        F1* f;
        Methods method;

    };

    struct SystemTask {
        double t0;
        double t1;
        double* y0;
        ui n;
        double h;
        F1System* f;
        Methods method;

    };

    DIFFURSOLVERDLL_API int SolveDiffUr(const Task& task, double& res);

    DIFFURSOLVERDLL_API int  SolveDiffUrArr(const Task& task, double* res);

    DIFFURSOLVERDLL_API int  SolveDiffUrAutoH(const Task& task, double eps, double& res);

    DIFFURSOLVERDLL_API int  SolveDiffUrAutoHArr(const Task& task, double eps, double* resX, double* resY,int& resSize);

    DIFFURSOLVERDLL_API int  SolveDiffUrAutoMethod(const Task& task, double eps, double& res);

    DIFFURSOLVERDLL_API int  SolveDiffUrAutoMethodArr(const Task& task, double eps, double* res);

    DIFFURSOLVERDLL_API int  SolveDiffUrAutoMethodArr2(const Task& task, double eps, double* res, Methods* methodsArr);

    DIFFURSOLVERDLL_API int  SolveDiffUrSystem(const SystemTask &task, double* res);

    DIFFURSOLVERDLL_API int  SolveDiffUrSystemArr(const SystemTask& task, double* res);

    

}