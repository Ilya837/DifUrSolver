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

double AutoHStep(double x0, double y0,double& h, F1 f, StepF stepF, double p, double eps, bool& more)
{
	double res1;
	double res2;
	double R;

	double eps2 = eps / pow(2, p + 1);


	while (true) {
		res1 = stepF(x0, y0, h, f);
		res2 = stepF(x0, y0, h / 2, f);
		res2 += stepF(x0 + h / 2, y0 + res2, h / 2, f);
		R = (res2 - res1) / (pow(2, p) - 1);
		if (abs(R) > eps) {
			h /= 2;
		}
		else if (abs(R) <= eps2) {
			more = true;
			return res1;
		}
		else break;
	}

	more = false;
	return res1;
}

double AutoMethodStep(double x0, double y0, double h, F1 f, double eps,bool& lessMethod,Methods& method)
{
	double res1;
	double res2;
	double res0;

	while (true) {

		switch (method)
		{
			case Euler: {
				res1 = EulerStep(x0, y0, h, f);
				res2 = RK2Step(x0, y0, h, f);

				if (abs(res2 - res1) > eps) { // условие повышения точности
					method = RK2;
					continue;
				}
				else {
					return res1;
				}

				break;
			}
			case RK2: {
				res0 = EulerStep(x0, y0, h, f);
				res1 = RK2Step(x0, y0, h, f);
				res2 = RK4Step(x0, y0, h, f);

				if (abs(res1 - res0) < eps / pow(2,2 + 1)) { // условие понижения точности
					lessMethod = true;
					return res1;
				}
				else if (abs(res2 - res1) > eps) { // условие повышения точности
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
				if (abs(res2 - res1) < eps / pow(2,4+1)) { // условие понижения точности
					lessMethod = true;
					return res1;
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


void EulerSystemStep(double x0, double* y0, double* yStep, ui n, double h, F1System* f)
{
	for (ui i = 0; i < n; i++) {
		yStep[i] = h * f[i](x0, y0);
	}
}

void RK2SystemStep(double x0,const double* y0, double* yStep, double* ytmp, ui n, double h, F1System* f)
{

	for (ui i = 0; i < n; i++) {
		yStep[i] = f[i](x0, y0);
	}

	for (ui i = 0; i < n; i++) {
		ytmp[i] = y0[i] + h * yStep[i];
		
	}

	for (ui i = 0; i < n; i++) {
		yStep[i] = h * (yStep[i] + f[i](x0 + h, ytmp)) / 2;
	}
}

void RK4SystemStep(double x0,const double* y0, double* yStep, double* ytmp, double* ytmp2, ui n, double h, F1System* f) {

	
	for (ui i = 0; i < n; i++) {
		yStep[i] = ytmp[i] = f[i](x0, y0);
	}

	for (ui i = 0; i < n; i++) {
		ytmp2[i] = y0[i] + h * ytmp[i] / 2;
	}

	for (ui i = 0; i < n; i++) {
		yStep[i] += 2 * (ytmp[i] = f[i](x0 + h / 2, ytmp2));
	}

	for (ui i = 0; i < n; i++) {
		ytmp2[i] = y0[i] + h * ytmp[i] / 2;
	}

	for (ui i = 0; i < n; i++) {
		yStep[i] += 2 * (ytmp[i] = f[i](x0 + h / 2, ytmp2));
	}

	for (ui i = 0; i < n; i++) {
		ytmp2[i] = y0[i] + h * ytmp[i];
	}

	for (ui i = 0; i < n; i++) {
		yStep[i] += f[i](x0 + h, ytmp2);
		yStep[i] = yStep[i] * h / 6;
	}

}





extern "C" {
	
	DIFFURSOLVERDLL_API int SolveDiffUr(const Task& task, double& res)
	{
		double(*stepF)(double x0, double y0, double h, F1 f);

		switch (task.method)
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
		double xNow = task.x0;
		res = task.y0;

		while (task.x0 + task.h < task.x1) {
			res += stepF(xNow, res, task.h, (*task.f));
			i++;
			xNow = task.x0 + i * task.h;
		}
		res += stepF(xNow, res, task.x1 - xNow, (*task.f));
		return 0;
	}

	DIFFURSOLVERDLL_API int SolveDiffUrAutoH(const Task& task, double eps, double& res)
	{
		StepF stepF;
		double p;

		switch (task.method)
		{
		case Euler:
			stepF = EulerStep;
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
		res = task.y0;
		double xNow = task.x0;
		double hNow = task.h;
		while (true) {
			bool more = false;
			double tmp = AutoHStep(xNow, res, hNow, *task.f, stepF, p, eps, more);

			if ((xNow + hNow > task.x1)) {
				if ((xNow + hNow - task.x1) < eps) {
					res += tmp;
					break;
				}

				switch (task.method)
				{
				case Euler:
					res += EulerStep(xNow, res, task.x1 - xNow, *task.f);
					break;
				case RK2:
					res += RK2Step(xNow, res, task.x1 - xNow, *task.f);
					break;
				case RK4:
					res += RK4Step(xNow, res, task.x1 - xNow, *task.f);
					break;
				default:
					return -1;
					break;
				}
				break;
			}

			xNow += hNow;
			if (more) hNow *= 2;
			res += tmp;
		}

		return 0;
	}

	DIFFURSOLVERDLL_API int SolveDiffUrAutoHArr(const Task& task, double eps, double* resX, double* resY, int& resSize)
	{

		StepF stepF;
		double p;

		switch (task.method)
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
		double hNow = task.h;
		resX[i] =task.x0;
		resY[i] =task.y0;
		while (true) {
			bool more = false;
			double tmp = AutoHStep(resX[i], resY[i], hNow, *task.f, stepF, p, eps, more);
			if ((resX[i] + hNow > task.x1)) {
				

				if ((resX[i] + hNow - task.x1) < eps) {
					i++;
					resX[i] = resX[i-1] + hNow;
					resY[i] = resY[i-1] + tmp;
					resSize = i + 1;
					return 0;
				}

				double LastRes = 0;

				switch (task.method)
				{
				case Euler:
					LastRes = resY[i] + EulerStep(resX[i], resY[i], task.x1 - resX[i], *task.f);
					break;
				case RK2:
					LastRes = resY[i] + RK2Step(resX[i], resY[i], task.x1 - resX[i], *task.f);
					break;
				case RK4:
					LastRes = resY[i] + RK4Step(resX[i], resY[i], task.x1 - resX[i], *task.f);
					break;
				default:
					return -1;
					break;
				}

				i++;
				resX[i] = task.x1;
				resY[i] = LastRes;
				resSize = i + 1;

				return 0;
			}

			i++;
			resX[i] =resX[i-1] + hNow;
			resY[i] = resY[i-1] + tmp;

			if (more) hNow *= 2;

			
		}

		return 0;
	}

	DIFFURSOLVERDLL_API int SolveDiffUrAutoMethod(const Task& task, double eps, double& res)
	{
		int i = 0;
		double xNow = task.x0; 
		Methods methodNow = task.method;
		bool lessMethod = false;
		while (xNow + task.h < task.x1) {
			lessMethod = false;
			res +=AutoMethodStep(task.x0,res,task.h,*task.f,eps,lessMethod, methodNow);
			if (lessMethod) {
				switch (methodNow)
				{
				case Euler:
					break;
				case RK2:
					methodNow = Euler;
					break;
				case RK4:
					methodNow = RK2;
					break;
				default:
					break;
				}
			}
			i++;
			xNow = task.x0 + i * task.h;
		}
		res += AutoMethodStep(xNow, res, task.x1-xNow, *task.f, eps, lessMethod, methodNow);
		return 0;
	}

	DIFFURSOLVERDLL_API int SolveDiffUrAutoMethodArr(const Task& task, double eps, double* res)
	{

		int i = 0;
		double xNow = task.x0;
		bool lessMethod = false;
		Methods methodNow = task.method;

		res[i] = task.y0;

		double yAdd = 0;

		while (xNow + task.h < task.x1) {
			lessMethod = false;
			yAdd = AutoMethodStep(xNow, res[i], task.h, (*task.f), eps, lessMethod, methodNow);
			if (lessMethod) {
				switch (methodNow)
				{
				case Euler:
					break;
				case RK2:
					methodNow = Euler;
					break;
				case RK4:
					methodNow = RK2;
					break;
				default:
					break;
				}
			}
			i++;
			res[i] = res[i-1] + yAdd;
			xNow = task.x0 + i * task.h;
		}
		i++;
		res[i] = res[i-1] + AutoMethodStep(xNow, res[i-1], task.x1 - xNow, *task.f, eps, lessMethod, methodNow);
		return 0;
	}

	DIFFURSOLVERDLL_API int SolveDiffUrAutoMethodArr2(const Task& task, double eps, double* res, Methods* methodsArr)
	{

		int i = 0;
		double xNow = task.x0;
		double yAdd = 0;
		Methods methodNow = task.method;

		res[i] = task.y0;
		methodsArr[i] = task.method;
		bool lessMethod = false;
		while (xNow + task.h < task.x1) {
			lessMethod = false;
			yAdd = AutoMethodStep(xNow, res[i], task.h, *task.f, eps, lessMethod, methodNow);
			i++;
			res[i] = res[i-1] + yAdd;
			methodsArr[i] = methodNow;


			if (lessMethod) {
				switch (methodNow)
				{
				case Euler:
					break;
				case RK2:
					methodNow = Euler;
					break;
				case RK4:
					methodNow = RK2;
					break;
				default:
					break;
				}
			}
			
			xNow = task.x0 + i * task.h;
		}
		i++;
		res[i] = res[i-1] + AutoMethodStep(xNow, res[i - 1], task.x1 - xNow, *task.f, eps, lessMethod, methodNow);
		methodsArr[i] = methodNow;
		return 0;
	}

	DIFFURSOLVERDLL_API int SolveDiffUrArr(const Task& task, double* res)
	{
		StepF stepF;

		switch (task.method)
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
		double xNow = task.x0;
		res[i] = task.y0;
		while (xNow + task.h < task.x1) {
			i++;
			res[i] =res[i-1] + stepF(xNow, res[i-1], task.h, *task.f);
			
			xNow = task.x0 + i * task.h;

		}
		i++;
		res[i] = res[i-1] + stepF(xNow, res[i-1], task.x1 - xNow, *task.f);

		return 0;
	}

	DIFFURSOLVERDLL_API int SolveDiffUrSystem(const SystemTask &task, double* res)
	{
		double* yStep = new double[task.n];
		double* yTmp = new double[task.n];
		double* yTmp2 = new double[task.n];
		int k = 0;
		double tNow = task.t0;

		switch (task.method)
		{
		case Euler:

			for (int i = 0; i < task.n; i++) {
				res[i] = task.y0[i];
			}

			while (tNow + task.h < task.t1) {

				EulerSystemStep(tNow,res,yStep,task.n,task.h,task.f);
				for (ui i = 0; i < task.n; i++) {
					res[i] += yStep[i];
				}
				k++;
				tNow = task.t0 + task.h * k;
			}

			EulerSystemStep(tNow, res, yStep, task.n, task.t1 - tNow, task.f);
			for (ui i = 0; i < task.n; i++) {
				res[i] += yStep[i];
			}
			break;
		case RK2:

			for (int i = 0; i < task.n; i++) {
				res[i] = task.y0[i];
			}

			while (tNow + task.h < task.t1) {

				RK2SystemStep(tNow, res, yStep, yTmp, task.n, task.h, task.f);
				for (ui i = 0; i < task.n; i++) {
					res[i] += yStep[i];
				}
				k++;
				tNow = task.t0 + task.h * k;
			}

			RK2SystemStep(tNow, res, yStep,yTmp, task.n, task.t1 - tNow, task.f);
			for (ui i = 0; i < task.n; i++) {
				res[i] += yStep[i];
			}
			break;
		case RK4:

			for (int i = 0; i < task.n; i++) {
				res[i] = task.y0[i];
			}

			while (tNow + task.h < task.t1) {

				RK4SystemStep(tNow, res, yStep,yTmp,yTmp2, task.n, task.h, task.f);
				for (ui i = 0; i < task.n; i++) {
					res[i] += yStep[i];
				}
				k++;
				tNow = task.t0 + task.h * k;
			}

			RK4SystemStep(tNow, res, yStep, yTmp, yTmp2, task.n, task.t1 - tNow, task.f);
			for (ui i = 0; i < task.n; i++) {
				res[i] += yStep[i];
			}
			break;
		default:
			return -1;
			break;
		}

		return 0;
	}

	DIFFURSOLVERDLL_API int SolveDiffUrSystemArr(const SystemTask& task, double* res)
	{
		double* yStep = new double[task.n];
		double* yTmp = new double[task.n];
		double* yTmp2 = new double[task.n];
		int k = 0;
		double tNow = task.t0;


		for (int i = 0; i < task.n; i++) {
			res[k * task.n + i] = task.y0[i];

		}
		k++;

		switch (task.method)
		{
		case Euler:


			while (tNow + task.h < task.t1) {

				EulerSystemStep(tNow, res + (k-1)*task.n, yStep, task.n, task.h, task.f);
				for (ui i = 0; i < task.n; i++) {
					res[k*task.n +i] = res[(k-1) * task.n + i] + yStep[i];
				}
				k++;
				tNow = task.t0 + task.h * k;
			}

			EulerSystemStep(tNow, res + (k - 1) * task.n, yStep, task.n, task.t1 - tNow, task.f);
			for (ui i = 0; i < task.n; i++) {
				res[k * task.n + i] = res[(k - 1) * task.n + i] + yStep[i];
			}
			break;
		case RK2:


			while (tNow + task.h < task.t1) {

				RK2SystemStep(tNow, res, yStep, yTmp, task.n, task.h, task.f);
				for (ui i = 0; i < task.n; i++) {
					res[k * task.n + i] = res[(k - 1) * task.n + i] + yStep[i];
				}
				k++;
				tNow = task.t0 + task.h * k;
			}

			RK2SystemStep(tNow, res, yStep, yTmp, task.n, task.t1 - tNow, task.f);
			for (ui i = 0; i < task.n; i++) {
				res[k * task.n + i] = res[(k - 1) * task.n + i] + yStep[i];
			}
			break;
		case RK4:


			while (tNow + task.h < task.t1) {

				RK4SystemStep(tNow, res, yStep, yTmp, yTmp2, task.n, task.h, task.f);
				for (ui i = 0; i < task.n; i++) {
					res[k * task.n + i] = res[(k - 1) * task.n + i] + yStep[i];
				}
				k++;
				tNow = task.t0 + task.h * k;
			}

			RK4SystemStep(tNow, res, yStep, yTmp, yTmp2, task.n, task.t1 - tNow, task.f);
			for (ui i = 0; i < task.n; i++) {
				res[k * task.n + i] = res[(k - 1) * task.n + i] + yStep[i];
			}
			break;
		default:
			return -1;
			break;
		}

		return 0;
	}


}