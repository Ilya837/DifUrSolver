#include "DifUrSolver.h"
#include <iostream>
#include <iomanip>
#include <vector>

#define ui unsigned int
#define vecF std::vector<double (*)(double x, double* y)>
#define vecRealF std::vector<double (*)(double x)>

double DifUrSolver::EulerStep(double x0, double y0, double h, double(*f)(double x, double y))
{
	return f(x0, y0) * h;
}

double DifUrSolver::RK2Step(double x0, double y0, double h, double(*f)(double x, double y))
{

	//return h * f(x0 + h/2, y0 + h * f(x0,y0)/2);

	double tmp = f(x0, y0);
	return h * (tmp + f(x0 + h, y0 + h * tmp)) / 2;
}

double DifUrSolver::RK4Step(double x0, double y0, double h, double(*f)(double x, double y))
{
	double sum = 0;
	double n = f(x0, y0);
	sum += n;

	n = f(x0 + h / 2, y0 + h * n / 2);
	sum += 2 * n;

	n = f(x0 + h / 2, y0 + h * n / 2);
	sum += 2 * n;

	sum += f(x0 + h, y0 + h * n);

	return h * sum / 6;
}

void DifUrSolver::EulerSystemStep(double x0, double* y0, double* yres, ui n, double h, vecF f)
{
	for (ui i = 0; i < n; i++) {
		yres[i] = h * f.at(i)(x0, y0);
	}
}

void DifUrSolver::RK2SystemStep(double x0, double* y0, double* yres, double* ytmp, ui n, double h, vecF f)
{

	for (ui i = 0; i < n; i++) {
		ytmp[i] = f.at(i)(x0, y0);
	}

	for (ui i = 0; i < n; i++) {
		y0[i] += h * ytmp[i];
	}

	for (ui i = 0; i < n; i++) {
		yres[i] = h * (ytmp[i] + f.at(i)(x0 + h, y0)) / 2;
	}

	for (ui i = 0; i < n; i++) {
		y0[i] -= h * ytmp[i];
	}
}

void DifUrSolver::RK4SystemStep(double x0, double* y0, double* yres, double* ytmp, double* ytmp2, ui n, double h, vecF f) {

	for (ui i = 0; i < n; i++) {
		yres[i] = 0;
	}

	for (ui i = 0; i < n; i++) {
		yres[i] += ytmp[i] = f.at(i)(x0, y0);
	}

	for (ui i = 0; i < n; i++) {
		ytmp2[i] = y0[i] + h * ytmp[i] / 2;
	}

	for (ui i = 0; i < n; i++) {
		yres[i] += 2 * (ytmp[i] = f.at(i)(x0 + h / 2, ytmp2));
	}

	for (ui i = 0; i < n; i++) {
		ytmp2[i] = y0[i] + h * ytmp[i] / 2;
	}

	for (ui i = 0; i < n; i++) {
		yres[i] += 2 * (ytmp[i] = f.at(i)(x0 + h / 2, ytmp2));
	}

	for (ui i = 0; i < n; i++) {
		ytmp2[i] = y0[i] + h * ytmp[i];
	}

	for (ui i = 0; i < n; i++) {
		yres[i] += f.at(i)(x0 + h, ytmp2);
		yres[i] = yres[i] * h / 6;
	}

}

int DifUrSolver::SolveDiffUr(double x0, double y0, double x1, double h, double(*f)(double x, double y),double &res, Methods method)
{
	double(*stepF)(double x0, double y0, double h, double(*f)(double x, double y));
	
	switch (method)
	{
	case Euler:
		stepF = &EulerStep;
		break;
	case RK2:
		stepF = &RK2Step;
		break;
	case RK4:
		stepF = &RK4Step;
		break;
	default:
		return -1;
		break;
	}
	int i = 0;
	while (x0 + h < x1) {
		y0 += stepF(x0, y0, h, f);
		i++;
		x0 += i * h;
	}
	res = y0 + stepF(x0, y0, x1 - x0, f);
	return 0;
}

int DifUrSolver::SolveDiffUrArr(double x0, double y0, double x1, double h, double(*f)(double x, double y), double* res, Methods method)
{
	double(*stepF)(double x0, double y0, double h, double(*f)(double x, double y));

	switch (method)
	{
	case Euler:
		stepF = &EulerStep;
		break;
	case RK2:
		stepF = &RK2Step;
		break;
	case RK4:
		stepF = &RK4Step;
		break;
	default:
		return -1;
		break;
	}
	int i = 0;
	while (x0 + h < x1) {
		res[i] = y0;
		y0 += stepF(x0, y0, h, f);
		i++;
		x0 = i * h;
		
	}
	res[i] = y0;
	res[i+1] = y0 + stepF(x0, y0, x1 - x0, f);
	return 0;
}

int DifUrSolver::SolveDiffUrSystm(double x0, double* y0, double x1, ui n, double h, vecF f, double* res, Methods method)
{
	double* ystep = new double[n];
	double* ytmp = new double[n];
	double* ytmp2 = new double[n];
	int k = 0;
	switch (method)
	{
	case Euler:

		while (x0 + h < x1) {

			EulerSystemStep(x0, y0, ystep, n, h, f);
			for (ui i = 0; i < n; i++) {
				y0[i] += ystep[i];
			}
			k++;
			x0 = h * k;
		}

		EulerSystemStep(x0, y0, ystep, n, x1 - x0, f);
		for (ui i = 0; i < n; i++) {
			res[i] = y0[i] + ystep[i];
		}
		break;
	case RK2:

		while (x0 + h < x1) {
			RK2SystemStep(x0, y0, ystep, ytmp, n, h, f);
			for (ui i = 0; i < n; i++) {
				y0[i] += ystep[i];
			}
			k++;
			x0 = h * k;
		}

		RK2SystemStep(x0, y0, ystep, ytmp, n, x1 - x0, f);
		for (ui i = 0; i < n; i++) {
			res[i] = y0[i] + ystep[i];
		}
		break;
	case RK4:

		while (x0 + h < x1) {
			RK4SystemStep(x0, y0, ystep, ytmp, ytmp2, n, h, f);
			for (ui i = 0; i < n; i++) {
				y0[i] += ystep[i];
			}
			k++;
			x0 = h * k;
		}

		RK4SystemStep(x0, y0, ystep, ytmp, ytmp2, n, x1 - x0, f);
		for (ui i = 0; i < n; i++) {
			res[i] = y0[i] + ystep[i];
		}
		break;
	default:
		return -1;
		break;
	}

	return 0;
}

int DifUrSolver::SolveDiffUrSystmArr(double x0, double* y0, double x1, ui n, double h, vecF f, double** res, Methods method)
{
	double* ystep = new double[n];
	double* ytmp = new double[n];
	double* ytmp2 = new double[n];

	int k = 0;
	switch (method)
	{
	case Euler:
		while (x0 + h + 0.00000001 < x1) {

			EulerSystemStep(x0, y0, ystep, n, h, f);
			for (ui i = 0; i < n; i++) {
				res[k][i] = y0[i];
				y0[i] += ystep[i];
			}
			
			k++;
			x0 = h * k;
		}
		
		EulerSystemStep(x0, y0, ystep, n, x1 - x0, f);
		for (ui i = 0; i < n; i++) {
			res[k][i] = y0[i];
			res[k+1][i] = y0[i] + ystep[i];
		}
		break;
	case RK2:
		while (x0 + h + 0.00000001 < x1) {
			RK2SystemStep(x0, y0, ystep, ytmp, n, h, f);
			for (ui i = 0; i < n; i++) {
				res[k][i] = y0[i];
				y0[i] += ystep[i];
			}
			k++;
			x0 = h * k;
			
		}

		RK2SystemStep(x0, y0, ystep, ytmp, n, x1 - x0, f);
		for (ui i = 0; i < n; i++) {
			res[k][i] = y0[i];
			res[k + 1][i] = y0[i] + ystep[i];
		}
		break;
	case RK4:
		while (x0 + h + 0.00000001 < x1) {
			RK4SystemStep(x0, y0, ystep, ytmp, ytmp2, n, h, f);
			for (ui i = 0; i < n; i++) {
				res[k][i] = y0[i];
				y0[i] += ystep[i];
			}
			
			k++;
			x0 = h * k;
		}

		RK4SystemStep(x0, y0, ystep, ytmp, ytmp2, n, x1 - x0, f);
		for (ui i = 0; i < n; i++) {
			res[k][i] = y0[i];
			res[k + 1][i] = y0[i] + ystep[i];
		}
		break;
	default:
		return -1;
		break;
	}

	return 0;
}

void DifUrSolver::Test(double x0, double y0, double x1, double h,
	double(*f)(double x, double y), double(*realY)(double x), Methods method, std::string header) {
	
	double(*stepF)(double x0, double y0, double h, double(*f)(double x, double y));
		
	switch (method)
	{
	case Euler:
		stepF = EulerStep;
		break;
	case RK2:
		stepF = RK2Step;
		break;
	case RK4:
		stepF = RK4Step;
		break;
	default:
		return;
		break;
	}

	double localErr[2] = {0,0};
	double GlobalErr[2] = { 0,0};

	int coutCount = 10;

	std::cout << "                             " << header << std::endl;

	for (int i = 0; i < 2; i++) {

		double x = x0, y = y0;
		double Y = y;
		double h2 = h / std::pow(2,i);

		while (x + h2 < x1) {
			localErr[i] = std::max(std::abs(Y + stepF(x, Y, h2, f) - realY(x + h2)), localErr[i]);
			y += stepF(x, y, h2, f);

			x += h2;
			Y = realY(x);

			GlobalErr[i] = std::max(std::abs(Y - y), GlobalErr[i]);
		}


		std::cout << "h: " << std::fixed << std::setprecision(coutCount) << h2 << " localErr: " << std::fixed << std::setprecision(coutCount) << localErr[i]
						<< " GlobalErr: " << std::fixed << std::setprecision(coutCount) << GlobalErr[i] << std::endl;

	}

	std::cout << "                  O(h^" << std::fixed <<  std::setprecision(coutCount) << std::log2(localErr[0]/ localErr[1]) 
					<< ")        O(h^" << std::fixed << std::setprecision(coutCount) << std::log2(GlobalErr[0] / GlobalErr[1]) <<")\n";
}


void DifUrSolver::TestSystem(double x0, double* y0, double x1, ui n, double h,
	vecF f, vecRealF RealY, Methods method, std::string header) {


	double localErr[2] = { 0,0 };
	double GlobalErr[2] = { 0,0 };

	int coutCount = 10;

	std::cout << "                             " << header << std::endl;

	double* y = new double[n];
	double* Y = new double[n];
	double* resY = new double[n];

	double* tmpY = new double[n];
	double* tmpY2 = new double[n];

	for (int i = 0; i < 2; i++) {

		double x = x0;

		for (ui i = 0; i < n; i++) {
			Y[i] = y[i] = y0[i];
		}

		double h2 = h / std::pow(2, i);

		while (x + h2 < x1) {

			switch (method)
			{
			case Euler:
				EulerSystemStep(x, Y, resY, n, h2, f);
				break;
			case RK2:
				RK2SystemStep(x, Y, resY,tmpY, n, h2, f);
				break;
			case RK4:
				RK4SystemStep(x, Y, resY,tmpY,tmpY2, n, h2, f);
				break;
			default:
				return;
				break;
			}

			double lr = 0;

			for (ui i = 0; i < n; i++) {
				lr += pow(Y[i] + resY[i] - RealY[i](x + h2),2);
			}

			lr = sqrt(lr);
			
			localErr[i] = std::max(lr, localErr[i]);

			switch (method)
			{
			case Euler:
				EulerSystemStep(x, y, resY, n, h2, f);
				break;
			case RK2:
				RK2SystemStep(x, y, resY, tmpY, n, h2, f);
				break;
			case RK4:
				RK4SystemStep(x, y, resY, tmpY, tmpY2, n, h2, f);
				break;
			default:
				return;
				break;
			}

			for (ui i = 0; i < n; i++) {
				y[i] += resY[i];
			}

			x += h2;
			for (ui i = 0; i < n; i++) {
				Y[i] = RealY.at(i)(x);
			}

			lr = 0;

			for (ui i = 0; i < n; i++) {
				lr += pow(Y[i] - y[i], 2);
			}

			lr = sqrt(lr);

			GlobalErr[i] = std::max(lr, GlobalErr[i]);
		}


		std::cout << "h: " << std::fixed << std::setprecision(coutCount) << h2 << " localErr: " << std::fixed << std::setprecision(coutCount) << localErr[i]
			<< " GlobalErr: " << std::fixed << std::setprecision(coutCount) << GlobalErr[i] << std::endl;

	}

	std::cout << "                  O(h^" << std::fixed << std::setprecision(coutCount) << std::log2(localErr[0] / localErr[1])
		<< ")        O(h^" << std::fixed << std::setprecision(coutCount) << std::log2(GlobalErr[0] / GlobalErr[1]) << ")\n";
}