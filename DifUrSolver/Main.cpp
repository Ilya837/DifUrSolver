#include <iostream>
#include <vector>
#include <iomanip>
#include <float.h>
#include "DifUrSolver.h"
#include "GraphDrawer.h"

double func(double x, double y) {
    return cos(x);
}

double realFunc(double x) {
    return sin(x);
}

double f1(double x, double* y) {
    return y[0] - y[1];
}

double f2(double x, double* y) {
    return - 4 * y[0] + y[1];
}

double f1Real(double x) {
    return ( - exp(3 * x) + exp(-x)) /4;
}

double f2Real(double x) {
    return (exp(3 * x) + exp(-x)) / 2;
}

void getArrFromF(double x0, double x1, double h, double* resx, double* resy, double (*f)(double x)) {
    int i = 0;
    while (x0 + h + 0.0000001 < x1) {
        resy[i] = f(x0);
        resx[i] = x0;
        i++;
        x0 = i * h;
    }
    resy[i] = f(x0);
    resx[i] = x0;

    resy[i+1] = f(x1);
    resx[i+1] = x1;
}

void getMaxMin(double* arr,ui n,double& max, double& min) {
    max = DBL_MIN;
    min = DBL_MAX;
    for (ui i = 0; i < n; i++) {
        if (arr[i] > max) {
            max = arr[i];
        }
        if (arr[i] < min) {
            min = arr[i];
        }
    }
}

void ex1(Methods m) {
    DifUrSolver* difSolv = new DifUrSolver();
    GraphDrawer* drawer = new GraphDrawer();

    double x0 = 0;
    double x1 = 7;
    double h = 0.1;

    ui size = ui(ceil((x1 - x0) / h) + 1);

    double* res = new double[size];

    difSolv->SolveDiffUrArr(x0, 0, x1, h, func, res, m);
    /*for (int i = 0; i <= double(1 - 0) / 0.05; i++) {
        std::cout <<0+ 0.05*i <<"  " << res[i] << "  " << sin(0 + 0.05 * i) << std::endl;
    }*/

    double* resx = new double[size * 3 - 1];
    double* resy = new double[size * 3 - 1];

    int sum = 0;

    for (ui i = 0; i < size; i++) {
        resx[i] = i * h;
        resy[i] = res[i];
    }
    sum += size;

    getArrFromF(x0, x1, h / 2, resx + sum, resy + sum, realFunc);

    double xmin = DBL_MAX;
    double xmax = DBL_MIN;
    double ymin = DBL_MAX;
    double ymax = DBL_MIN;

    getMaxMin(resx, size * 3 - 2, xmax, xmin);
    getMaxMin(resy, size * 3 - 2, ymax, ymin);


    ui ns[2] = { size,size * 2 - 2 };

    std::string colors[2] = { "red","green" };
    std::string legend[2] = { "red","green" };

    drawer->DrawGraph(resx, resy, ns, 2, xmin, xmax, 1, ymin, ymax, 1, colors, legend);
}

void ex2() {
    DifUrSolver* difSolv = new DifUrSolver();

    difSolv->Test(0, 0, 6, 0.1, func, realFunc, Euler, "Euler");
    difSolv->Test(0, 0, 6, 0.1, func, realFunc, RK2, "RK2");
    difSolv->Test(0, 0, 6, 0.1, func, realFunc, RK4, "RK4");
}

void ex3(Methods m) {
    /*
            x' = x - y
            y' = y - 4x

            x(0) = 0
            y(0) = 1

            ans:

            x = 1/4 (-e^(3t) + e^(-t))
            y = 1/2 (e^(3t) + e^(-t))

            x(1) ~ -4.92941
            y(0) ~ 10.2267
    */

    DifUrSolver* difSolv = new DifUrSolver();
    GraphDrawer* drawer = new GraphDrawer();


    double x0 = 0;
    double x1 = 1;
    double h = 0.1;
    ui size = ui(ceil((x1 - x0) / h) + 1);

    double** res = new double* [size];
    for (ui i = 0; i < size; i++) {
        res[i] = new double[2];
    }

    std::vector<double (*)(double x, double* y)> f;
    f.push_back(f1);
    f.push_back(f2);

    std::vector<double (*)(double x)> realF;
    realF.push_back(f1Real);
    realF.push_back(f2Real);

    double* y = new double[2];
    y[0] = 0;
    y[1] = 1;

    difSolv->SolveDiffUrSystmArr(x0, y, x1, 2, h, f, res, m);

    double* resx = new double[size * 3 - 1];
    double* resy = new double[size * 3 - 1];

    int sum = 0;

    for (ui i = 0; i < size; i++) {
        resx[i] = res[i][0];
        resy[i] = res[i][1];
    }
    sum += size;

    for (ui i = 0; i < size * 2 - 2; i++) {
        resx[i + sum] = realF[0](i * 0.05);
        resy[i + sum] = realF[1](i * 0.05);
    }

    double xmin = DBL_MAX;
    double xmax = DBL_MIN;
    double ymin = DBL_MAX;
    double ymax = DBL_MIN;

    getMaxMin(resx, size * 3 - 2, xmax, xmin);
    getMaxMin(resy, size * 3 - 2, ymax, ymin);


    ui ns[2] = { size,size * 2 - 2 };

    std::string colors[2] = { "red","green" };
    std::string legend[2] = { "red","green" };

    drawer->DrawGraph(resx, resy, ns, 2, xmin, xmax, 1, ymin, ymax, 1, colors, legend);

    /*std::cout << "t          x      realx     y     realy" << std::endl;
    for (int i = 0; i <= double(1 - 0) / 0.05; i++) {
        std::cout << std::fixed << std::setprecision(4) << 0 + 0.05 * i << "  ";
        for (int k = 0; k < 2; k++) {
            std::cout << std::fixed << std::setprecision(4) << res[i][k] << "  ";
            std::cout << std::fixed << std::setprecision(4) << realF[k](0 + 0.05*i) << "  ";
        }

        std::cout<< std::endl;
    }*/

}

void ex4() {
    DifUrSolver* difSolv = new DifUrSolver();

    std::vector<double (*)(double x, double* y)> f;
    f.push_back(f1);
    f.push_back(f2);

    std::vector<double (*)(double x)> realF;
    realF.push_back(f1Real);
    realF.push_back(f2Real);

    double* y = new double[2];
    y[0] = 0;
    y[1] = 1;

    difSolv->TestSystem(0, y, 1, 2, 0.001, f, realF, Euler, "Euler");
    difSolv->TestSystem(0, y, 1, 2, 0.001, f, realF, RK2, "RK2");
    difSolv->TestSystem(0, y, 1, 2, 0.001, f, realF, RK4, "RK4");
}

int main()
{
    ex1(Euler);
    ex2();
    //ex3(Euler);
    //ex4();

}