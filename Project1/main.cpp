#include <iostream>
#include <vector>
#include <iomanip>
#include <float.h>
#include "../DiffUrSolverDll2/DiffUrSolverDLL2.h"
#include "omp.h"

double func(double x, double y) {
    return cos(x);
}

double realFunc(double x) {
    return sin(x);
}

double func2(double x, double y) {
    return (2 * x * cos(x) - sin(x))/ (2 * x * sqrt(x));
}

double realFunc2(double x) {
    return sin(x) / sqrt(x);
}

double f1(double x, double* y) {
    return y[0] - y[1];
}

double f2(double x, double* y) {
    return -4 * y[0] + y[1];
}

double f1Real(double x) {
    return (-exp(3 * x) + exp(-x)) / 4;
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

    resy[i + 1] = f(x1);
    resx[i + 1] = x1;
}

void getArrFromF(double x0, double x1, double h, double* resx, double* resy, double (*xf)(double x), double (*yf)(double x, double y)) {
    int i = 0;
    while (x0 + h + 0.0000001 < x1) {
        resy[i] = yf(x0, 0);
        resx[i] = xf(x0);
        i++;
        x0 = i * h;
    }
    resy[i] = yf(x0, 0);
    resx[i] = xf(x0);

    resy[i + 1] = yf(x1, 0);
    resx[i + 1] = xf(x0);
}

void getMaxMin(double* arr, ui n, double& max, double& min) {
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

    double x0 = 0;
    double x1 = 7;
    double h = 0.01;

    ui size = ui(ceil((x1 - x0) / h) + 1);

    std::vector<double> res = std::vector<double>();
    res.clear();

    SolveDiffUrArr(x0, 0, x1, h, func, res, m);
    for (int i = 0; i < res.size(); i++) {
        std::cout <<x0+ h*i <<"  " << res[i] << "  " << realFunc(x0 + h*i) << std::endl;
    }

    /*double* resx = new double[size * 3 - 1];
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

    drawer->DrawGraph(resx, resy, ns, 2, xmin, xmax, 1, ymin, ymax, 1, colors, legend);*/
}

//void ex1_5(Methods m) {
//    DifUrSolver* difSolv = new DifUrSolver();
//    GraphDrawer* drawer = new GraphDrawer();
//
//    double x0 = 0;
//    double x1 = 7;
//    double h = 0.1;
//
//    ui size = ui(ceil((x1 - x0) / h) + 1);
//
//    std::vector<double> res = std::vector<double>();
//
//    SolveDiffUrArr(x0, 0, x1, h, func, res, m);
//    for (int i = 0; i <= double(1 - 0) / 0.05; i++) {
//        std::cout <<0+ 0.05*i <<"  " << res[i] << "  " << sin(0 + 0.05 * i) << std::endl;
//    }
//
//    double* resx = new double[size * 3 - 1];
//    double* resy = new double[size * 3 - 1];
//
//    int sum = 0;
//
//    for (ui i = 0; i < size; i++) {
//        resx[i] = res[i];
//        resy[i] = func(i * h, 0);
//    }
//    sum += size;
//
//    getArrFromF(x0, x1, h / 2, resx + sum, resy + sum, realFunc, func);
//
//
//
//    double xmin = DBL_MAX;
//    double xmax = DBL_MIN;
//    double ymin = DBL_MAX;
//    double ymax = DBL_MIN;
//
//    getMaxMin(resx, size * 3 - 2, xmax, xmin);
//    getMaxMin(resy, size * 3 - 2, ymax, ymin);
//
//
//    ui ns[2] = { size,size * 2 - 2 };
//
//    std::string colors[2] = { "red","green" };
//    std::string legend[2] = { "red","green" };
//
//    drawer->DrawGraph(resx, resy, ns, 2, xmin, xmax, 1, ymin, ymax, 1, colors, legend);
//}
//
//void ex2() {
//    DifUrSolver* difSolv = new DifUrSolver();
//
//    difSolv->Test(0, 0, 6, 0.1, func, realFunc, Euler, "Euler");
//    difSolv->Test(0, 0, 6, 0.1, func, realFunc, RK2, "RK2");
//    difSolv->Test(0, 0, 6, 0.1, func, realFunc, RK4, "RK4");
//}
//
//void ex3(Methods m) {
//    /*
//            x' = x - y
//            y' = y - 4x
//
//            x(0) = 0
//            y(0) = 1
//
//            ans:
//
//            x = 1/4 (-e^(3t) + e^(-t))
//            y = 1/2 (e^(3t) + e^(-t))
//
//            x(1) ~ -4.92941
//            y(0) ~ 10.2267
//    */
//
//    DifUrSolver* difSolv = new DifUrSolver();
//    GraphDrawer* drawer = new GraphDrawer();
//
//
//    double x0 = 0;
//    double x1 = 1;
//    double h = 0.1;
//    ui size = ui(ceil((x1 - x0) / h) + 1);
//
//    double** res = new double* [size];
//    for (ui i = 0; i < size; i++) {
//        res[i] = new double[2];
//    }
//
//    std::vector<double (*)(double x, double* y)> f;
//    f.push_back(f1);
//    f.push_back(f2);
//
//    std::vector<double (*)(double x)> realF;
//    realF.push_back(f1Real);
//    realF.push_back(f2Real);
//
//    double* y = new double[2];
//    y[0] = 0;
//    y[1] = 1;
//
//    difSolv->SolveDiffUrSystemArr(x0, y, x1, 2, h, f, res, m);
//
//    double* resx = new double[size * 3 - 1];
//    double* resy = new double[size * 3 - 1];
//
//    int sum = 0;
//
//    for (ui i = 0; i < size; i++) {
//        resx[i] = res[i][0];
//        resy[i] = res[i][1];
//    }
//    sum += size;
//
//    for (ui i = 0; i < size * 2 - 2; i++) {
//        resx[i + sum] = realF[0](i * 0.05);
//        resy[i + sum] = realF[1](i * 0.05);
//    }
//
//    double xmin = DBL_MAX;
//    double xmax = DBL_MIN;
//    double ymin = DBL_MAX;
//    double ymax = DBL_MIN;
//
//    getMaxMin(resx, size * 3 - 2, xmax, xmin);
//    getMaxMin(resy, size * 3 - 2, ymax, ymin);
//
//
//    ui ns[2] = { size,size * 2 - 2 };
//
//    std::string colors[2] = { "red","green" };
//    std::string legend[2] = { "red","green" };
//
//    drawer->DrawGraph(resx, resy, ns, 2, xmin, xmax, 1, ymin, ymax, 1, colors, legend);
//
//    /*std::cout << "t          x      realx     y     realy" << std::endl;
//    for (int i = 0; i <= double(1 - 0) / 0.05; i++) {
//        std::cout << std::fixed << std::setprecision(4) << 0 + 0.05 * i << "  ";
//        for (int k = 0; k < 2; k++) {
//            std::cout << std::fixed << std::setprecision(4) << res[i][k] << "  ";
//            std::cout << std::fixed << std::setprecision(4) << realF[k](0 + 0.05*i) << "  ";
//        }
//
//        std::cout<< std::endl;
//    }*/
//
//}
//
//void ex4() {
//    DifUrSolver* difSolv = new DifUrSolver();
//
//    std::vector<double (*)(double x, double* y)> f;
//    f.push_back(f1);
//    f.push_back(f2);
//
//    std::vector<double (*)(double x)> realF;
//    realF.push_back(f1Real);
//    realF.push_back(f2Real);
//
//    double* y = new double[2];
//    y[0] = 0;
//    y[1] = 1;
//
//    difSolv->TestSystem(0, y, 1, 2, 0.001, f, realF, Euler, "Euler");
//    difSolv->TestSystem(0, y, 1, 2, 0.001, f, realF, RK2, "RK2");
//    difSolv->TestSystem(0, y, 1, 2, 0.001, f, realF, RK4, "RK4");
//}
//
//void ex5(Methods m) {
//    DifUrSolver* difSolv = new DifUrSolver();
//    double res = -50;
//    double t0 = 0, t1 = 0;
//    double x0 = 0, x1 = 6, y0 = 0;
//
//    t0 = omp_get_wtime();
//    difSolv->SolveDiffUrAutoH(x0, y0, x1, 0.3, 0.0001, func, res, m);
//    t1 = omp_get_wtime();
//    std::cout << "time : " << t1 - t0 << " sec" << std::endl;
//    std::cout << "result: " << res << " real result " << realFunc(x1) << std::endl;
//}
//
//void ex5full() {
//    std::cout << "Euler: " << std::endl;
//    ex5(Euler);
//    std::cout << "RK2: " << std::endl;
//    ex5(RK2);
//    std::cout << "RK4: " << std::endl;
//    ex5(RK4);
//}
//
void ex6(Methods m) {
    std::vector<double> resX;
    std::vector<double> resY;
    double t0 = 0, t1 = 0;
    double x0 = 0, x1 = 6, y0 = realFunc(x0);

    t0 = omp_get_wtime();
    SolveDiffUrAutoHArr(x0, y0, x1, 0.3, 0.000001, func, resX, resY, m);
    t1 = omp_get_wtime();
    std::cout << "time : " << t1 - t0 << " sec" << std::endl;
    for (int i = 0; i < resX.size(); i++) {
        std::cout << "X: " << std::fixed << std::setprecision(6) << resX[i] << " result: " << std::fixed << std::setprecision(6) <<
            resY[i] << " real result " << std::fixed << std::setprecision(6) << realFunc(resX[i]) << std::endl;
    }

}

void ex6full() {
    std::cout << "Euler: " << std::endl;
    ex6(Euler);
    std::cout << "RK2: " << std::endl;
    ex6(RK2);
    std::cout << "RK4: " << std::endl;
    ex6(RK4);
}

void ex7() {
    std::vector<double> res;
    double t0 = 0, t1 = 0;
    double x0 = 1, x1 = 20, y0 = realFunc(x0), h = 0.01;

    t0 = omp_get_wtime();
    SolveDiffUrAutoMethodArr(x0, y0, x1, h, 0.0001, func, res, Euler);
    t1 = omp_get_wtime();
    std::cout << "time : " << t1 - t0 << " sec" << std::endl;
    for (int i = 0; i < res.size(); i++) {
        std::cout << "X: " << std::fixed << std::setprecision(6) << x0 + i*h << " result: " << std::fixed << std::setprecision(6) <<
            res[i] << " real result " << std::fixed << std::setprecision(6) << realFunc(x0 + i * h) << std::endl;
    }

}

void ex8() {
    std::vector<double> res;
    std::vector<Methods> m;
    double t0 = 0, t1 = 0;
    double x0 = 1, x1 = 2, y0 = realFunc2(x0), h = 0.1;

    t0 = omp_get_wtime();
    SolveDiffUrAutoMethodArr2(x0, y0, x1, h, 0.00001, func2, res,m, Euler);
    t1 = omp_get_wtime();
    std::cout << "time : " << t1 - t0 << " sec" << std::endl;
    for (int i = 0; i < res.size(); i++) {
        std::cout << "X: " << std::fixed << std::setprecision(6) << x0 + i * h << " result: " << std::fixed << std::setprecision(6) <<
            res[i] << " real result " << std::fixed << std::setprecision(6) << realFunc2(x0 + i * h) <<" Method: " << m[i] << std::endl;
    }

}



int main()
{
    //ex1(Euler);
    //ex1_5(RK4);
    //ex2();
    //ex3(Euler);
    //ex4();
    //ex6full();
    ex8();







}