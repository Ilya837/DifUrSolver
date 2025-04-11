#pragma once
#include <iostream>
#include <vector>

#define ui unsigned int
#define vecF std::vector<double (*)(double x, double* y)>
#define vecRealF std::vector<double (*)(double x)>

enum Methods
{
	Euler,
	RK2,
	RK4
};

class DifUrSolver
{
public:

	int static SolveDiffUr(double x0, double y0,double x1,double h,	double (*f)(double x, double y), double& res, Methods method);

	int static SolveDiffUrArr(double x0, double y0, double x1, double h, double(*f)(double x, double y), double* res, Methods method);

	int static SolveDiffUrSystm(double x0, double* y0, double x1,ui n, double h,vecF f, double* res, Methods method);

	int static SolveDiffUrSystmArr(double x0, double* y0, double x1, ui n, double h, vecF f, double** res, Methods method);

	void Test(double x0, double y0, double x1, double h,
		double(*f)(double x, double y), double(*realY)(double x), Methods method, std::string header = "");

	void static TestSystem(double x0, double* y0, double x1, ui n, double h,
		vecF f, vecRealF RealY, Methods method, std::string header);

private:

	double static EulerStep(double x0, double y0, double h, double(*f)(double x, double y));

	double static RK2Step(double x0, double y0, double h, double(*f)(double x, double y));

	double static RK4Step(double x0, double y0, double h, double(*f)(double x, double y));

	void static EulerSystemStep(double x0, double* y0, double* yres, ui n, double h, vecF f);

	void static RK2SystemStep(double x0, double* y0, double* yres, double* ytmp, ui n, double h, vecF f);

	void static RK4SystemStep(double x0, double* y0, double* yres, double* ytmp, double* ytmp2, ui n, double h, vecF f);

};

