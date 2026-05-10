#include "pch.h"
#include "DiffUrSolverDLL2.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <algorithm>
#include <vector>
#include <limits>
#include <functional>
#include <omp.h>

typedef std::vector<std::function<double(double*)>> DiffFSystem;

double* MemoryTmp::stepArr = nullptr;
ui* MemoryTmp::arrUI1 = nullptr;
double* MemoryTmp::arr1 = nullptr;
double* MemoryTmp::arr2 = nullptr;
double* MemoryTmp::arr3 = nullptr;
double* MemoryTmp::matrix1 = nullptr;
double* MemoryTmp::matrix2 = nullptr;

ui MemoryTmp::NstepArr = 0;
ui MemoryTmp::NarrUI1 = 0;
ui MemoryTmp::Narr1 = 0;
ui MemoryTmp::Narr2 = 0;
ui MemoryTmp::Narr3 = 0;
ui MemoryTmp::Nmatrix1 = 0;
ui MemoryTmp::Nmatrix2 = 0;



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

double BackEulerStep(double x0, double y0, double h, F1 f)
{
	double diff_lim = 1e-12;
	int max_iter = 100000;

	double y_next = y0 + h * f(x0, y0);

	double x_next = x0 + h;

	auto F = [&](double yy) {
		return yy - y0 - h * f(x_next, yy);
		};

	for (int i = 0; i < max_iter;i++) {
		double Fy = F(y_next);

		double eps = 1e-8;

		double dFy = (F(y_next + eps) - Fy) / eps;

		if (dFy == 0) {
			break;
		}
		else {
			y_next -= Fy / dFy;
			double diff = std::abs(Fy / dFy);
			//std::cout << diff;

			if (diff < diff_lim) break;
		}

	}

	return y_next - y0;

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


void EulerSystemStep(double x0,const double* y0, double* yStep, ui n, double h, F1System* f)
{
	for (ui i = 0; i < n; i++) {
		yStep[i] = h * f[i](x0, y0);
	}
}

void RK2SystemStep(double x0,const double* y0, double* yStep, ui n, double h, F1System* f)
{

	for (ui i = 0; i < n; i++) {
		yStep[i] = f[i](x0, y0);
	}

	for (ui i = 0; i < n; i++) {
		MemoryTmp::arr1[i] = y0[i] + h * yStep[i];
		
	}

	for (ui i = 0; i < n; i++) {
		yStep[i] = h * (yStep[i] + f[i](x0 + h, MemoryTmp::arr1)) / 2;
	}
}

void RK4SystemStep(double x0,const double* y0, double* yStep, ui n, double h, F1System* f) {

	
	for (ui i = 0; i < n; i++) {
		yStep[i] = MemoryTmp::arr1[i] = f[i](x0, y0);
	}

	for (ui i = 0; i < n; i++) {
		MemoryTmp::arr2[i] = y0[i] + h * MemoryTmp::arr1[i] / 2;
	}

	for (ui i = 0; i < n; i++) {
		yStep[i] += 2 * (MemoryTmp::arr1[i] = f[i](x0 + h / 2, MemoryTmp::arr2));
	}

	for (ui i = 0; i < n; i++) {
		MemoryTmp::arr2[i] = y0[i] + h * MemoryTmp::arr1[i] / 2;
	}

	for (ui i = 0; i < n; i++) {
		yStep[i] += 2 * (MemoryTmp::arr1[i] = f[i](x0 + h / 2, MemoryTmp::arr2));
	}

	for (ui i = 0; i < n; i++) {
		MemoryTmp::arr2[i] = y0[i] + h * MemoryTmp::arr1[i];
	}

	for (ui i = 0; i < n; i++) {
		yStep[i] += f[i](x0 + h, MemoryTmp::arr2);
		yStep[i] = yStep[i] * h / 6;
	}

}

void J( double* y, double* res, ui n, DiffFSystem f) {
	double exp = 1e-8;

	#pragma omp parallel
	{
		double* copy_plus = new double[n];
		std::copy(y, y + n, copy_plus);

		double* copy_minus = new double[n];
		std::copy(y, y + n, copy_minus);

		#pragma omp for
		for (int j = 0; j < n; j++) {
			copy_plus[j] += exp;
			copy_minus[j] -= exp;
			for (int i = 0; i < n; i++) {
				double f_plus = f[i](copy_plus);
				double f_minus = f[i](copy_minus);
				res[i * n + j] = (f_plus - f_minus) / (2 * exp);
			}
			copy_plus[j] = y[j];
			copy_minus[j] = y[j];

		}

		delete[] copy_plus;
		delete[] copy_minus;
	}
}

// will break original matrix
void InverseMatrix(double* matrix, double* tmp, ui n) {
	
	std::memset(tmp, 0.0, n * n * sizeof(double));

	double** matrix2 = new double*[n];
	double** tmp2 = new double* [n];

	for (int i = 0; i < n; i++) {
		tmp[i * n + i] = 1.0;

		matrix2[i] = matrix + i * n;
		tmp2[i] = tmp + i * n;
	}

	for (ui i = 0; i < n; i++) {

		ui j = i;
		for (; j < n; j++) {
			if (matrix2[j][i] != 0) break;
		}

		if (j == n) {
			return;
		}

		if (j != i) {
			std::swap(matrix2[i], matrix2[j]);
			std::swap(tmp2[i], tmp2[j]);
		}

		double diag = matrix2[i][i];

		for (ui j = i; j < n; j++) {
			matrix2[i][j] /= diag;
			tmp2[i][j] /= diag;
		}
		matrix2[i][i] = 1;

		for (int k = i+1; k < n; k++) {

			double factor = matrix2[k][i];
			for (int j = 0; j < n; j++) {
				matrix2[k][j] -= factor * matrix2[i][j];
				tmp2[k][j] -= factor * tmp2[i][j];
			}
			
		}
	}

	for (int i = n-2; i >= 0; i--) {
		for (int j = i+1; j < n; j++) {
			tmp2[i][j] -= matrix2[i][j] * tmp2[j][j];
		}
	}

	for (ui i = 0; i < n; i++) {
		for (ui j = 0; j < n; j++) {
			matrix[i * n + j] = tmp2[i][j];
		}
	}

}

void LU(double* A, double* res, ui n) {

	std::copy(A, A + n * n, res);

#pragma omp parallel
	{

	#pragma omp for schedule(static,1)
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < i; j++) {
				res[i * n + j] /= res[i * n + i];
			}
		}

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < i; j++) {
			#pragma omp for
				for (int k = i; k < i; k++) {
					res[i * n + k] -= res[i * n + j] * res[j * n + k];
				}
			}
		}
	}

}

void countP(ui* P, double* A, ui n) {
	double global_max = 0;
	ui global_idx = 0;


	for (int i = 0; i < n; i++) {
		P[i] = i;
	}

#pragma omp parallel if( n > 5000) 
	{

		for (int i = 0; i < n; i++) {

			double local_max = 0;
			ui local_idx = 0;

#pragma omp for 
			for (int j = i; j < n; j++) {
				if (std::abs(A[P[j] * n + i]) > local_max) {
					local_max = std::abs(A[P[j] * n + i]);
					local_idx = j;
				}
			}

#pragma omp barrier

#pragma omp critical
			{
				if (local_max > global_max) {
					global_max = local_max;
					global_idx = local_idx;
				}
			}

#pragma omp single
			{
				if (i != local_idx) {
					std::swap(P[i], P[local_idx]);
				}
			}

#pragma omp barrier

		}
	}


}

void  countLU(ui start, ui finish, double* LUMatrix, ui n) {

	for (int i = start; i < finish; i++) {

		double* i_row = &LUMatrix[i * n];

		for (int l = start; l < i; l++) {

			double* l_row = &LUMatrix[l * n];

			const double tmp = i_row[l];

			for (int k = i; k < finish; k++) {
				i_row[k] -= l_row[k] * tmp;
			}
		}

		for (int k = i + 1; k < finish; k++) {

			double* k_row = &LUMatrix[k * n];

			for (int l = start; l < i; l++) {
				k_row[i] -= k_row[l] * LUMatrix[l * n + i];
			}

			k_row[i] /= i_row[i];
		}
	}
}

void countU(ui i1, ui i2, ui j1, ui j2, double* LUMatrix, ui n) {

	if (j2 > n)	j2 = n;

	double sum = 0;

	for (int i = i1; i < i2; i++) {

		double* row_i = LUMatrix + i * n;

		for (int k = i1; k < i; k++) {

			double* row_k = LUMatrix + k * n;

			const double tmp = row_i[k];

			for (int j = j1; j < j2; j++) {
				row_i[j] -= tmp * row_k[j];
			}

		}
	}


}

void countL(ui i1, ui i2, ui j1, ui j2, double* LUMatrix, ui n) {

	if (i2 > n) i2 = n;
	double  sum = 0;

	for (int j = j1; j < j2; j++) {

		double* row_j = &LUMatrix[j * n];

		for (int i = i1; i < i2; i++) {

			double* row_i = &LUMatrix[i * n];
			sum = 0;

			for (int k = j1; k < j; k++) {
				sum += row_i[k] * LUMatrix[k * n + j];
			}

			row_i[j] = (row_i[j] - sum) / row_j[j];
		}
	}

}

void countA(ui i1, ui i2, ui j1, ui  j2, ui k1, ui k2, double* LUMatrix, ui n) {

	if (i2 > n) i2 = n;
	if (j2 > n) j2 = n;
	if (k2 > n) k2 = n;

	for (int i = i1; i < i2; i++) {

		double* row_i = &LUMatrix[i * n];

		for (int k = k1; k < k2; k++) {

			double* row_k = &LUMatrix[k * n];

			const double tmp = row_i[k];

			for (int j = j1; j < j2; j++) {
				row_i[j] -= tmp * row_k[j];
			}
		}
	}

}

void PLU(double* A, ui* P, double* LUMatrix, ui n, ui block_size) {

	countP(P, A, n);

#pragma omp parallel for num_threads(2)
	for (int i = 0; i < n; i++) {
		std::copy(A + P[i] * n, A + P[i] * n + n, LUMatrix + i * n);
	}

	int stage = 0;

	for (stage = 0; (stage + 1) * block_size < n; stage++) {

		int start = stage * block_size;

		countLU(start, start + block_size, LUMatrix, n);

#pragma omp parallel
		{
#pragma omp for nowait schedule(dynamic, 1)
			for (int j1 = start + block_size; j1 < n; j1 += block_size) {
				countU(start, start + block_size, j1, j1 + block_size, LUMatrix, n);
			}

#pragma omp for schedule(dynamic, 1)
			for (int i1 = start + block_size; i1 < n; i1 += block_size) {

				countL(i1, i1 + block_size, start, start + block_size, LUMatrix, n);
			}


#pragma omp for collapse(2)  schedule(dynamic, 1)
			for (int i1 = start + block_size; i1 < n; i1 += block_size) {
				for (int j1 = start + block_size; j1 < n; j1 += block_size) {
					countA(i1, i1 + block_size, j1, j1 + block_size, start, start + block_size, LUMatrix, n);
				}
			}

		}

	}

	int start = stage * block_size;

	countLU(start, n, LUMatrix, n);
	
}

void Axb(double* A,double* b, double* res, ui n) {

	double* LUmatrix = MemoryTmp::matrix2;

	ui* P = MemoryTmp::arrUI1;

	PLU(A, P, LUmatrix, n,80);

	double* y = MemoryTmp::arr2;

	for (int i = 0; i < n; i++) {
		y[i] = b[P[i]];

		double sum = 0;

		double* row_i = &LUmatrix[i * n];

		for (int j = 0; j < i; j++) {
			sum += y[j] * row_i[j];
		}

		y[i] -= sum;
	}

	for (int i = n - 1; i >= 0; i--) {
		res[i] = y[i];

		double sum = 0;

		double* row_i = &LUmatrix[i * n];

		for (int j = i + 1; j < n; j++) {
			sum += row_i[j] * res[j];
		}

		res[i] = (res[i] - sum) / row_i[i];
	}

}

void BackEulerSystemStep(double x0, const double* y0, double* yStep, ui n, double h, F1System* f) {
	double diff_lim = 1e-12;
	int max_iter = 10;
	double x_next = x0 + h;

	double* A = MemoryTmp::matrix1;

	double sum = 0;
	double eps = 1e-6;

	DiffFSystem system(n);

	for (int i = 0; i < n; i++) {
		yStep[i] = y0[i] + h * f[i](x0,y0);
	}

	for (int i = 0; i < n; i++) {
		system[i] = [&,i](double* yy) {
			return yy[i] - y0[i] - h * f[i](x_next, yy);
			};
	}

	double* b = MemoryTmp::arr1;

	double* stepdel = MemoryTmp::arr3;

	for (int i = 0; i < max_iter; i++) {

		J(yStep, A, n, system);

		for (int i = 0; i < n; i++) {
			b[i] = -system[i](yStep);
		}

		Axb(A, b, stepdel, n);

		for (int j = 0; j < n; j++) {
			yStep[j] += stepdel[j];
		}

		bool flag = false;

		for (ui j = 0; j < n; j++) {
			double diff = std::abs(yStep[j] - h * f[j](x_next, yStep) - y0[j]);
			if (diff > diff_lim) {
				flag = true;
				break;
			}
		}

		if (!flag) break;

	}

	for (int i = 0; i < n; i++) {
		yStep[i] -= y0[i];
	}
	
}

extern "C" {

	DIFFURSOLVERDLL_API void InitMemory(ui n) {
		MemoryTmp::initStepArr(n);
		MemoryTmp::initArrUI1(n);
		MemoryTmp::initArr1(n);
		MemoryTmp::initArr2(n);
		MemoryTmp::initArr3(n);
		MemoryTmp::initMatrix1(n);
		MemoryTmp::initMatrix2(n);
	}

	DIFFURSOLVERDLL_API void SetThreadCount(ui n) {
		omp_set_num_threads(n);
	}
	
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
		case BackEuler:
			stepF = &BackEulerStep;
			break;
		default:
			return -1;
			break;
		}
		int i = 0;
		double xNow = task.x0;
		res = task.y0;

		while (xNow + task.h < task.x1) {
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
		case BackEuler:
			stepF = &BackEulerStep;
			p = 1;
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
		case BackEuler:
			stepF = &BackEulerStep;
			p = 1;
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
		case BackEuler:
			stepF = &BackEulerStep;
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
		MemoryTmp::initStepArr(task.n);
		double* yStep = MemoryTmp::stepArr;
		int k = 0;
		double tNow = task.t0;

		for (int i = 0; i < task.n; i++) {
			res[i] = task.y0[i];
		}
		k++;

		StepSystemF stepF;

		switch (task.method)
		{
		case Euler: {
			stepF = &EulerSystemStep;
			break;
		}
		case RK2: {
			MemoryTmp::initArr1(task.n);
			stepF = &RK2SystemStep;
			break;
		}
		case RK4: {

			MemoryTmp::initArr1(task.n);
			MemoryTmp::initArr2(task.n);
			stepF = &RK4SystemStep;
			break;
		}
		case BackEuler: {

			MemoryTmp::initArr1(task.n);
			MemoryTmp::initArr2(task.n);
			MemoryTmp::initArr3(task.n);
			MemoryTmp::initMatrix1(task.n);
			MemoryTmp::initMatrix2(task.n);
			stepF = &BackEulerSystemStep;
			break;
		}
		default: {
			return -1;
			break;
		}
		}

		while (tNow + task.h < task.t1) {

			stepF(tNow, res, yStep, task.n, task.h, task.f);
			for (ui i = 0; i < task.n; i++) {
				res[i] += yStep[i];
			}

			tNow = task.t0 + task.h * k;
			k++;
		}

		stepF(tNow, res, yStep, task.n, task.t1 - tNow, task.f);
		for (ui i = 0; i < task.n; i++) {
			res[i] += yStep[i];
		}

		return 0;
	}

	DIFFURSOLVERDLL_API int SolveDiffUrSystemArr(const SystemTask& task, double* res)
	{
		MemoryTmp::initStepArr(task.n);
		double* yStep = MemoryTmp::stepArr;
		
		int k = 0;
		double tNow = task.t0;


		for (int i = 0; i < task.n; i++) {
			res[k * task.n + i] = task.y0[i];
		}
		k++;

		StepSystemF stepF;

		switch (task.method)
		{
		case Euler: {

			stepF = &EulerSystemStep;
			break;
		}
		case RK2: {

			MemoryTmp::initArr1(task.n);
			stepF = &RK2SystemStep;
			break;
		}
		case RK4: {

			MemoryTmp::initArr1(task.n);
			MemoryTmp::initArr2(task.n);
			stepF = &RK4SystemStep;
			break;
		}
		case BackEuler: {

			MemoryTmp::initArr1(task.n);
			MemoryTmp::initArr2(task.n);
			MemoryTmp::initArr3(task.n);
			MemoryTmp::initArrUI1(task.n);
			MemoryTmp::initMatrix1(task.n);
			MemoryTmp::initMatrix2(task.n);
			stepF = &BackEulerSystemStep;
			break;
		}
		default:
			return -1;
			break;
		}


		while (tNow + task.h < task.t1) {

			stepF(tNow, res + (k - 1) * task.n, yStep, task.n, task.h, task.f);
			for (ui i = 0; i < task.n; i++) {
				res[k * task.n + i] = res[(k - 1) * task.n + i] + yStep[i];
			}

			tNow = task.t0 + task.h * k;
			k++;
		}

		stepF(tNow, res + (k - 1) * task.n, yStep, task.n, task.t1 - tNow, task.f);
		for (ui i = 0; i < task.n; i++) {
			res[k * task.n + i] = res[(k - 1) * task.n + i] + yStep[i];
		}


		return 0;
	}


}