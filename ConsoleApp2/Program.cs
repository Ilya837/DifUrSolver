using System;
using System.Runtime.InteropServices;
using DiffUrSolver;


public enum Methods
{
    Euler,
    RK2,
    RK4
}

class Program
{
    

    // Пример функции для решения
    public static double SampleFunction(double x, double y)
    {
        return Math.Cos(x);
    }

    public static double RealSampleFunction(double x)
    {
        return Math.Sin(x);
    }

    public static double SampleFunction2(double x, double[] y)
    {
        return x * x + y[0];
    }

    public static double SampleFunction3(double x0, double y0, double h, DiffUrSolver.DiffUrSolver.F1 f)
    {
        return 0;
    }

    public static double SampleFunction4(double x)
    {
        return x * x;
    }


    static void experiment1(Methods method)
    {
        double x0 = 0.0;
        double y0 = 0.0;
        double x1 = 1.0;
        double h = 0.1;
        double result = 0;

        DiffUrSolver.DiffUrSolver.F1 func = SampleFunction;


        int status = DiffUrSolver.DiffUrSolver.SolveDiffUr(x0, y0, x1, h, func, ref result,method);

        Console.WriteLine(result.ToString());
        Console.WriteLine(RealSampleFunction(x1));
    }

    static void experiment2(Methods method)
    {
        double x0 = 0.0;
        double y0 = 0.0;
        double x1 = 6.0;
        double h = 0.3;
        double eps = 0.0001;
        double[] resultX = new double[1000];
        double[] resultY = new double[1000];
        int resultSize = 0;

        // Решение ОДУ методом Эйлера
        DiffUrSolver.DiffUrSolver.F1 func = SampleFunction;
        DiffUrSolver.DiffUrSolver.F1System f1System = SampleFunction2;
        DiffUrSolver.DiffUrSolver.StepF stepF = SampleFunction3;
        DiffUrSolver.DiffUrSolver.F f = SampleFunction4;


        int status = DiffUrSolver.DiffUrSolver.SolveDiffUrAutoHArr(x0,y0,x1,h,eps,func,resultX,resultY,ref resultSize,method);

        //List<double> finalResult = result.TakeWhile(x => x != 0).ToList();

        for(int i =0; i< resultSize; i++)
        {
            Console.WriteLine("{0:F7}|{1:F7}|{2:F7}|{3:F7}", resultX[i], resultY[i], RealSampleFunction(resultX[i]), RealSampleFunction(resultX[i]) - resultY[i]);
                
        }

    }

    static void Main(string[] args)
    {
        //experiment1(Methods.RK4);
        experiment2(Methods.RK4);

    }
}
