#include "pch.h"
#include "DiffUrSolverDLL2.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>



double EulerStep(double x0, double y0, double h, F1 f)
{
	return f(x0, y0) * h;
}

double RK2Step(double x0, double y0, double h, F1 f)
{

	//return h * f(x0 + h/2, y0 + h * f(x0,y0)/2);

	double tmp = f(x0, y0);
	return h * (tmp + f(x0 + h, y0 + h * tmp)) / 2;
}

double RK4Step(double x0, double y0, double h, F1 f)
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

double AutoHStep(double x0, double y0, double& h, F1 f, StepF stepF, double p, double eps)
{
	double res1;
	double res2;
	double R;

	double eps2;
	bool isLowing = true;

	eps2 = eps / 128;

	res1 = stepF(x0, y0, h, f);
	res2 = stepF(x0, y0, h / 2, f);
	res2 += stepF(x0 + h / 2, y0 + res2, h / 2, f);
	R = (res2 - res1) / (pow(2, p) - 1);
	if (abs(R) > eps) {
		h /= 2;
		isLowing = true;
	}
	else if (abs(R) <= eps2) {
		h *= 2;
		isLowing = false;
	}
	else return res1;

	while (true) {
		res1 = stepF(x0, y0, h, f);
		res2 = stepF(x0, y0, h / 2, f);
		res2 += stepF(x0 + h / 2, y0 + res2, h / 2, f);
		R = (res2 - res1) / (pow(2, p) - 1);
		if (abs(R) > eps) {
			if (!isLowing) break;
			h /= 2;
		}
		else if (abs(R) <= eps2) {
			if (isLowing) break;
			h *= 2;
		}
		else break;
	}

	return res1;
}

double AutoMethodStep(double x0, double y0, double h, F1 f, double eps,Methods& method)
{
	double res1;
	double res2;
	bool isLowing = true;

	double k = 512;

	

	switch (method)
	{
		case Euler: {
			res1 = EulerStep(x0, y0, h, f);
			res2 = RK2Step(x0, y0, h, f);

			if (abs(res2 - res1) > eps) { // условие повышения точности
				method = RK2;
				isLowing = false;
			}
			else {
				return res1;
			}

			break;
		}
		case RK2: {
			res1 = RK2Step(x0, y0, h, f);
			res2 = RK4Step(x0, y0, h, f);

			if (abs(res2 - res1) < eps / k) { // условие понижения точности
				method = Euler;
				isLowing = true;
			}
			else if (abs(res2 - res1) / 3 > eps) { // условие повышения точности
				method = RK4;
				isLowing = false;
			}
			else {
				return res1;
			}

			break;
		}
		case RK4: {
			res1 = RK4Step(x0, y0, h, f);
			res2 = RK2Step(x0, y0, h, f);
			if (abs(res2 - res1) < eps / k) { // условие понижения точности
				method = RK2;
				isLowing = true;
			}
			else {
				return res1;
			}
			break;

		}
		default:
			return -1;
			break;
	}

	


	while (true) {

		switch (method)
		{
			case Euler: {
				res1 = EulerStep(x0, y0, h, f);
				res2 = RK2Step(x0, y0, h, f);

				if (abs(res2 - res1) > eps) { // условие повышения точности
					if (isLowing) return res1;
					method = RK2;
					continue;
				}
				else {
					return res1;
				}

				break;
			}
			case RK2: {
				res1 = RK2Step(x0, y0, h, f);
				res2 = RK4Step(x0, y0, h, f);

				if (abs(res2 - res1) < eps / k) { // условие понижения точности
					if (!isLowing) return res1;
					method = Euler;
					continue;
				}
				else if (abs(res2 - res1) / 3 > eps) { // условие повышения точности
					if (isLowing) return res1;
					method = RK4;
					continue;
				}
				else {
					return res1;
				}

				break;
			}
			case RK4: {
				res1 = RK4Step(x0, y0, h, f);
				res2 = RK2Step(x0, y0, h, f);
				if (abs(res2 - res1) < eps / k) { // условие понижения точности
					if (!isLowing) return res1;
					method = RK2;
					continue;
				}
				else {
					return res1;
				}
				break;

			}
		}

	}

	return res1;
}


void EulerSystemStep(double x0, double* y0, double* yres, ui n, double h, F1System* f)
{
	for (ui i = 0; i < n; i++) {
		yres[i] = h * f[i](x0, y0);
	}
}

void RK2SystemStep(double x0, double* y0, double* yres, double* ytmp, ui n, double h, F1System* f)
{

	for (ui i = 0; i < n; i++) {
		ytmp[i] = f[i](x0, y0);
	}

	for (ui i = 0; i < n; i++) {
		y0[i] += h * ytmp[i];
	}

	for (ui i = 0; i < n; i++) {
		yres[i] = h * (ytmp[i] + f[i](x0 + h, y0)) / 2;
	}

	for (ui i = 0; i < n; i++) {
		y0[i] -= h * ytmp[i];
	}
}

void RK4SystemStep(double x0, double* y0, double* yres, double* ytmp, double* ytmp2, ui n, double h, F1System* f) {

	for (ui i = 0; i < n; i++) {
		yres[i] = 0;
	}

	for (ui i = 0; i < n; i++) {
		yres[i] += ytmp[i] = f[i](x0, y0);
	}

	for (ui i = 0; i < n; i++) {
		ytmp2[i] = y0[i] + h * ytmp[i] / 2;
	}

	for (ui i = 0; i < n; i++) {
		yres[i] += 2 * (ytmp[i] = f[i](x0 + h / 2, ytmp2));
	}

	for (ui i = 0; i < n; i++) {
		ytmp2[i] = y0[i] + h * ytmp[i] / 2;
	}

	for (ui i = 0; i < n; i++) {
		yres[i] += 2 * (ytmp[i] = f[i](x0 + h / 2, ytmp2));
	}

	for (ui i = 0; i < n; i++) {
		ytmp2[i] = y0[i] + h * ytmp[i];
	}

	for (ui i = 0; i < n; i++) {
		yres[i] += f[i](x0 + h, ytmp2);
		yres[i] = yres[i] * h / 6;
	}

}





extern "C" {
	
	DIFFURSOLVERDLL_API int SolveDiffUr(double x0, double y0, double x1, double h, F1 f, double& res, Methods method)
	{
		double(*stepF)(double x0, double y0, double h, F1 f);

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
		double realx0 = x0;
		while (x0 + h < x1) {
			y0 += stepF(x0, y0, h, f);
			i++;
			x0 = realx0 + i * h;
		}
		res = y0 + stepF(x0, y0, x1 - x0, f);
		return 0;
	}

	DIFFURSOLVERDLL_API int SolveDiffUrAutoH(double x0, double y0, double x1, double h, double eps, F1 f, double& res, Methods method)
	{
		StepF stepF;
		double p;

		switch (method)
		{
		case Euler:
			stepF =EulerStep;
			p = 1;
			break;
		case RK2:
			stepF = RK2Step;
			p = 2;
			break;
		case RK4:
			stepF = RK4Step;
			p = 4;
			break;
		default:
			return -1;
			break;
		}
		int i = 0;
		while (true) {
			double tmp = AutoHStep(x0, y0, h, f, stepF, p, eps);
			if (x0 + h > x1) {
				if ((x0 + h - x1) < eps) {
					res = y0 + tmp;
					break;
				}

				switch (method)
				{
				case Euler:
					res = y0 + EulerStep(x0, y0, x1 - x0, f);
					break;
				case RK2:
					res = y0 + RK2Step(x0, y0, x1 - x0, f);
					break;
				case RK4:
					res = y0 + RK4Step(x0, y0, x1 - x0, f);
					break;
				default:
					return -1;
					break;
				}

				break;
			}

			y0 += tmp;
			x0 += h;
		}

		return 0;
	}

	DIFFURSOLVERDLL_API int SolveDiffUrAutoHArr(double x0, double y0, double x1, double h, double eps, F1 f, double* resX, double* resY,int& resSize, Methods method)
	{

		StepF stepF;
		double p;

		switch (method)
		{
		case Euler:
			stepF =EulerStep;
			p = 1;
			break;
		case RK2:
			stepF = RK2Step;
			p = 2;
			break;
		case RK4:
			stepF = RK4Step;
			p = 4;
			break;
		default:
			return -1;
			break;
		}
		int i = 0;
		resX[i] =x0;
		resY[i] =y0;
		while (true) {
			double tmp = AutoHStep(x0, y0, h, f, stepF, p, eps);
			if (x0 + h > x1) {
				i++;
				if ((x0 + h - x1) < eps) {
					resX[i] = x0 + h;
					resY[i] =y0 + tmp;
					resSize = i + 1;
					return 0;
				}

				double LastRes = 0;

				switch (method)
				{
				case Euler:
					LastRes = y0 + EulerStep(x0, y0, x1 - x0, f);
					break;
				case RK2:
					LastRes = y0 + RK2Step(x0, y0, x1 - x0, f);
					break;
				case RK4:
					LastRes = y0 + RK4Step(x0, y0, x1 - x0, f);
					break;
				default:
					return -1;
					break;
				}
				resX[i] =x1;
				resY[i] = LastRes;
				resSize = i + 1;

				return 0;
			}

			y0 += tmp;
			x0 += h;
			i++;
			resX[i] =x0;
			resY[i] =y0;

			
		}

		return 0;
	}

	DIFFURSOLVERDLL_API int SolveDiffUrAutoMethod(double x0, double y0, double x1, double h, double eps, F1 f, double& res, Methods method)
	{
		int i = 0;
		double realx0 = x0;
		while (x0 + h < x1) {
			y0 +=AutoMethodStep(x0,y0,h,f,eps,method);
			i++;
			x0 = realx0 + i * h;
		}
		res = y0 + AutoMethodStep(x0, y0, x1-x0, f, eps, method);
		return 0;
	}

	DIFFURSOLVERDLL_API int SolveDiffUrAutoMethodArr(double x0, double y0, double x1, double h, double eps, F1 f, double* res, int& sizeRes, Methods method)
	{

		int i = 0;
		double realx0 = x0;

		res[i] =y0;

		while (x0 + h < x1) {
			y0 += AutoMethodStep(x0, y0, h, f, eps, method);
			
			i++;
			res[i] = y0;
			x0 = realx0 + i * h;
		}
		i++;
		res[i] =y0 + AutoMethodStep(x0, y0, x1 - x0, f, eps, method);
		sizeRes = i + 1;
		return 0;
	}

	DIFFURSOLVERDLL_API int SolveDiffUrAutoMethodArr2(double x0, double y0, double x1, double h, double eps, F1 f, double* res, int& sizeRes, Methods* methodsArr, Methods method)
	{

		int i = 0;
		double realx0 = x0;

		res[i] =y0;
		methodsArr[i]=method;

		while (x0 + h < x1) {
			y0 += AutoMethodStep(x0, y0, h, f, eps, method);
			i++;
			res[i] = y0;
			methodsArr[i] = method;
			
			x0 = realx0 + i * h;
		}
		i++;
		res[i] = y0 + AutoMethodStep(x0, y0, x1 - x0, f, eps, method);
		methodsArr[i] = method;
		sizeRes = i + 1;
		return 0;
	}


	DIFFURSOLVERDLL_API int SolveDiffUrArr(double x0, double y0, double x1, double h, F1 f, double* res, int& sizeRes, Methods method)
	{
		StepF stepF;

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
		double realx0 = x0;
		res[i] = y0;
		while (x0 + h < x1) {
			
			y0 += stepF(x0, y0, h, f);
			i++;
			res[i] = y0;
			x0 = realx0 + i * h;

		}
		i++;
		res[i] = y0 + stepF(x0, y0, x1 - x0, f);

		sizeRes = i + 1;
		return 0;
	}

	DIFFURSOLVERDLL_API int SolveDiffUrSystem(double x0, double* y0, double x1, ui n, double h, F1System* f, double* res, Methods method)
	{
		double* ystep = new double[n];
		double* ytmp = new double[n];
		double* ytmp2 = new double[n];
		int k = 0;
		double realx0 = 0;
		switch (method)
		{
		case Euler:

			while (x0 + h < x1) {

				EulerSystemStep(x0, y0, ystep, n, h, f);
				for (ui i = 0; i < n; i++) {
					y0[i] += ystep[i];
				}
				k++;
				x0 = realx0 + h * k;
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
				x0 = realx0 + h * k;
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
				x0 = realx0 + h * k;
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

	DIFFURSOLVERDLL_API int SolveDiffUrSystemArr(double x0, double* y0, double x1, ui n, double h, F1System* f, double** res, int& resSize, Methods method)
	{
		double* ystep = new double[n];
		double* ytmp = new double[n];
		double* ytmp2 = new double[n];

		int k = 0;
		double realx0 = 0;
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
				x0 = realx0 + h * k;
			}

			EulerSystemStep(x0, y0, ystep, n, x1 - x0, f);
			for (ui i = 0; i < n; i++) {
				res[k][i] = y0[i];
				res[k + 1][i] = y0[i] + ystep[i];
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
				x0 = realx0 + h * k;

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
				x0 = realx0 + h * k;
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

		resSize = k + 1;

		return 0;
	}


}