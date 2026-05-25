#include "pch.h"
#include <vector>
#include <new>
#include <exprtk.hpp>

#ifdef DIFFURSOLVERDLL_EXPORTS
#define DIFFURSOLVERDLL_API __declspec(dllexport)
#else
#define DIFFURSOLVERDLL_API __declspec(dllimport)
#endif

typedef unsigned int ui;

enum Methods
{
    Euler,
    RK2,
    RK4,
    BackEuler,
	RadoIIA3,
	RadoIA3,
	RadoIIA5
};

static struct MemoryTmp {

	static double* stepArr;
	static ui NstepArr;
	static void initStepArr(ui n) {
		if (n > NstepArr) {
			delete[] stepArr;
			stepArr = new double[n];
			NstepArr = n;
		}
	}

	static ui* arrUI1;
	static ui NarrUI1;
	static void initArrUI1(ui n) {
		if (n > NarrUI1) {
			delete[] arrUI1;
			arrUI1 = new ui[n];
			NarrUI1 = n;
		}
	}


	static double* arr1;
	static ui Narr1;
	static void initArr1(ui n) {
		if (n > Narr1) {
			delete[] arr1;
			arr1 = new double[n];
			Narr1 = n;
		}
	}

	static double* arr2;
	static ui Narr2;
	static void initArr2(ui n) {
		if (n > Narr2) {
			delete[] arr2;
			arr2 = new double[n];
			Narr2 = n;
		}
	}

	static double* arr3;
	static ui Narr3;
	static void initArr3(ui n) {
		if (n > Narr3) {
			delete[] arr3;
			arr3 = new double[n];
			Narr3 = n;
		}
	}

	static double* arr4;
	static ui Narr4;
	static void initArr4(ui n) {
		if (n > Narr4) {
			delete[] arr4;
			arr4 = new double[n];
			Narr4 = n;
		}
	}


	static double* matrix1;
	static ui Nmatrix1;
	static void initMatrix1(ui n) {
		if (n > Nmatrix1) {
			delete[] matrix1;
			matrix1 = new double[n * n];
			Nmatrix1 = n;
		}
	}

	static double* matrix2;
	static ui Nmatrix2;
	static void initMatrix2(ui n) {
		if (n > Nmatrix2) {
			delete[] matrix2;
			matrix2 = new double[n * n];
			Nmatrix2 = n;
		}
	}
};

static struct Config {
	static ui block_size;

	static void SetBlockSize(ui n) {
		if (n > 0) {
			block_size = n;
		}
	}
};


class CompiledFunction
{
public:

	using expression_t = exprtk::expression<double>;
	using parser_t = exprtk::parser<double>;
	using symbol_table_t = exprtk::symbol_table<double>;

private:

	double x_ = 0.0;

	exprtk::vector_view<double> y_;

	std::vector<exprtk::vector_view<double>> vecs_;

	symbol_table_t symbols_;
	expression_t expression_;
	parser_t parser_;

public:

	CompiledFunction(
		const char* expr,
		double* yBase,
		ui ySize,
		double* scal_consts,
		const char** scal_names,
		const ui scal_count,
		double* vec_consts,
		const ui* vec_lengts,
		const char** vec_names,
		const ui vec_count) : y_(yBase,ySize)
	{

		symbols_.add_variable("x", x_);

		symbols_.add_vector("y",y_);

		for (ui i = 0; i < scal_count; i++) {
			symbols_.add_variable(scal_names[i], scal_consts[i]);
		}

		
		

		ui ind = 0;

		vecs_.reserve(vec_count);

		for (ui i = 0; i < vec_count; i++) {

			vecs_.emplace_back(exprtk::vector_view<double>(vec_consts + ind, vec_lengts[i]));
			symbols_.add_vector(vec_names[i], vecs_[i]);
			ind += vec_lengts[i];
		}




		expression_.register_symbol_table(symbols_);

		if (!parser_.compile(std::string(expr), expression_))
		{
			throw std::runtime_error(parser_.error().c_str());
		}
	}

	double eval(
		double xv,
		double* y)
	{
		x_ = xv;

		y_.rebase(y);

		return expression_.value();
	}
};

typedef std::vector<std::unique_ptr<CompiledFunction>> F1System;



extern "C" {

    typedef double(__cdecl* F1)(double x, double y);

    //typedef double(__cdecl* F1System)(double x,const double* y);

    typedef double(__cdecl* F)(double x);

    typedef double(__cdecl* StepF)(double x0, double y0, double h, F1 f);

    typedef void(__cdecl* StepSystemF)(double x0, double* y0, double* yStep, ui n, double h, F1System& f);


    struct Task {
        double x0;
        double x1;
        double y0;
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
        const char** f;
		double** scal_consts;
		const char*** scal_names;
		const ui* scal_count;
		double** vec_consts;
		const ui** vec_lengts;
		const char*** vec_names;
		const ui* vec_count;
        Methods method;

    };

	DIFFURSOLVERDLL_API void InitMemory(ui n);

	DIFFURSOLVERDLL_API void SetThreadCount(ui n);

	DIFFURSOLVERDLL_API void SetBlockSizeCount(ui n);

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

