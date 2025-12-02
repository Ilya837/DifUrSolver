using System;
using System.Runtime.InteropServices;
using DiffUrSolver;
using OxyPlot;
using OxyPlot.Series;
using OxyPlot.Axes;
using OxyPlot.Core.Drawing;
using OxyPlot.WindowsForms;




class Program
{
    

    // Пример функции для решения
    public static double func1(double x, double y)
    {
        return Math.Cos(x);
    }

    public static double realFunc1(double x)
    {
        return Math.Sin(x);
    }

    public static double SampleFunction2(double x, double[] y)
    {
        return x * x + y[0];
    }

    public static double SampleFunction3(double x0, double y0, double h, DiffUrSolver.F1 f)
    {
        return 0;
    }

    public static double SampleFunction4(double x)
    {
        return x * x;
    }


    static void experiment1(DiffUrSolver.Methods method)
    {

        DiffUrSolver.Task task = new DiffUrSolver.Task();

        task.x0 = 0.0;
        task.y0 = 0.0;
        task.x1 = 1.0;
        task.h = 0.1;
        double result = 0;

        task.f = DiffUrSolver.DiffUrSolver.F1ToIntPtr( func1);
        task.method = method;

        int status = DiffUrSolver.DiffUrSolver.SolveDiffUr(task, ref result);

        Console.WriteLine(result.ToString());
        Console.WriteLine(realFunc1(task.x1));
    }

    static void experiment2(DiffUrSolver.Methods method)
    {

        DiffUrSolver.Task task = new DiffUrSolver.Task();

        task.x0 = 0.0;
        task.y0 = 0.0;
        task.x1 = 6.0;
        task.h = 0.3;
        double eps = 0.0001;
        double[] resultX = new double[1000];
        double[] resultY = new double[1000];
        int resultSize = 0;

        // Решение ОДУ методом Эйлера
        task.f = DiffUrSolver.DiffUrSolver.F1ToIntPtr(func1);
        task.method = method;

        int status = DiffUrSolver.DiffUrSolver.SolveDiffUrAutoHArr(task,eps,resultX,resultY,ref resultSize);

        //List<double> finalResult = result.TakeWhile(x => x != 0).ToList();

        


        var pm = new PlotModel();
        var ls = new LineSeries() { Title = "1", MarkerStroke = OxyColors.Black };
        var ls2 = new LineSeries() { Title = "2", MarkerStroke = OxyColors.Black };

        for (int i =0; i< resultSize; i++)
        {
            //Console.WriteLine("{0:F7}|{1:F7}|{2:F7}|{3:F7}", resultX[i], resultY[i], RealSampleFunction(resultX[i]), RealSampleFunction(resultX[i]) - resultY[i]);

            ls.Points.Add(new DataPoint(resultX[i], resultY[i]));
            ls2.Points.Add(new DataPoint(resultX[i], realFunc1(resultX[i])));


        }

        pm.Series.Add(ls);
        pm.Series.Add(ls2);

    }
     static void TestMethodsEps(DiffUrSolver.F1 func, DiffUrSolver.F realFunc, DiffUrSolver.Methods method)
    {

        DiffUrSolver.Task task = new DiffUrSolver.Task();

        task.x0 = 0.0;
        task.y0 = 0.0;
        task.x1 = 4.0;
        task.h = 0.1;

        task.f = DiffUrSolver.DiffUrSolver.F1ToIntPtr(func);
        task.method = method;
        
        double[] result = new double[10000];
        double[] realResult = new double[10000];
        double localRes = 0;

        for(int i=0; task.x0+ i*task.h <= task.x1; i++)
        {
            realResult[i] = realFunc(task.x0+ i*task.h);
        }

        DiffUrSolver.DiffUrSolver.SolveDiffUrArr(task, result);

        //for (int i = 0; task.x0 + i * task.h < task.x1; i++)
        //{
        //    Console.Write(result[i] + " ");
        //}

        //Console.WriteLine();

        //for (int i = 0; task.x0 + i * task.h < task.x1; i++)
        //{
        //    Console.Write(realResult[i] + " ");
        //}

        //Console.WriteLine();

        double globalError1 = -10000000;
        double localError1 = -10000000;

        double startx0 = task.x0;
        double starty0 = task.y0;
        double startx1 = task.x1;

        for (int i = 0; task.x0 + task.h < startx1; i++)
        {

            task.x0 = startx0 + task.h * i;
            task.y0 = realFunc(task.x0);
            task.x1 = task.x0 + task.h;

            if (Math.Abs( result[i] - realResult[i]) > globalError1)
            {
                globalError1 = Math.Abs(result[i] - realResult[i]);
            }

            DiffUrSolver.DiffUrSolver.SolveDiffUr(task,ref localRes);

            if (Math.Abs(localRes - realResult[i+1]) > localError1)
            {
                localError1 = Math.Abs(localRes - realResult[i + 1]);
            }
        }

        task.x0 = startx0;
        task.y0 = starty0;
        task.x1 = startx1;
        task.h /= 2;

        for (int i = 0; task.x0 + i * task.h <= task.x1; i++)
        {
            realResult[i] = realFunc(task.x0 + i * task.h);
        }

        DiffUrSolver.DiffUrSolver.SolveDiffUrArr(task, result);

        //for (int i = 0; task.x0 + i * task.h < task.x1; i++)
        //{
        //    Console.Write(result[i] + " ");
        //}

        //Console.WriteLine();

        //for (int i = 0; task.x0 + i * task.h < task.x1; i++)
        //{
        //    Console.Write(realResult[i] + " ");
        //}

        //Console.WriteLine();

        double globalError2 = -10000000;
        double localError2 = -10000000;

        startx0 = task.x0;
        startx1 = task.x1;
        starty0 = task.y0;

        for (int i = 0; task.x0 + task.h < startx1; i++)
        {

            task.x0 = startx0 + task.h * i;
            task.y0 = realFunc(task.x0);
            task.x1 = task.x0 + task.h;

            if (Math.Abs(result[i] - realResult[i]) > globalError2)
            {
                globalError2 = Math.Abs(result[i] - realResult[i]);
            }

            DiffUrSolver.DiffUrSolver.SolveDiffUr(task, ref localRes);

            if (Math.Abs(localRes - realResult[i+1]) > localError2)
            {
                localError2 = Math.Abs(localRes - realResult[i + 1]);
            }
        }

        task.x0 = startx0;
        task.y0 = starty0;
        task.x1 = startx1;
        task.h *= 2;

        Console.WriteLine("          |    Global error     |   Local error");

        Console.WriteLine(" h = {0:F2} |     {1:F10}    |   {2:F10}",task.h, globalError1, localError1);

        Console.WriteLine(" h = {0:F2} |     {1:F10}    |   {2:F10}", task.h/2, globalError2, localError2);

        Console.WriteLine("          |     h^{0:F10}  | h^{1:F10}", Math.Log2( globalError1/globalError2) ,Math.Log2(localError1 / localError2) );


    }

    static void FullTestMethodsEps(DiffUrSolver.F1 func, DiffUrSolver.F realFunc)
    {
        Console.WriteLine("                      Euler");
        TestMethodsEps(func,realFunc,DiffUrSolver.Methods.Euler);
        Console.WriteLine();
        Console.WriteLine("                       RK2");
        TestMethodsEps(func, realFunc, DiffUrSolver.Methods.RK2);
        Console.WriteLine();
        Console.WriteLine("                       RK4");
        TestMethodsEps(func, realFunc, DiffUrSolver.Methods.RK4);
        Console.WriteLine();
    }

    public static double SystemFunctionSingle(double x, IntPtr y)
    {

        double[] array = new double[1];
        Marshal.Copy(y, array, 0, 1);

        return Math.Cos(array[0]);
    }
    public static double SystemFunction1(double x, IntPtr y)
    {

        double[] array = new double[2];
        Marshal.Copy(y, array, 0, 2);


        return array[0] - 5 * array[1];
    }

    public static double SystemFunction2(double x, IntPtr y)
    {
        double[] array = new double[2];
        Marshal.Copy(y, array, 0, 2);

        return 5 * array[0] + array[1];
    }

    public static double RealSystemFunction1(double x)
    {
        return Math.Pow(Math.E, x) * (Math.Sin(-5 * x) + Math.Cos(-5 * x));
    }

    public static double RealSystemFunction2(double x)
    {
        return Math.Pow(Math.E, x) * (Math.Sin(5 * x) + Math.Cos(5 * x));
    }

    static void experiment3(DiffUrSolver.Methods method)
    {
        DiffUrSolver.SystemTask task = new DiffUrSolver.SystemTask();

        task.t0 = 0;
        task.t1 = 4.0;
        task.h = 0.001;
        task.method = method;
        task.n = 2;

        double[] y0 = new double[6];
        y0[0] = 1;
        y0[1] = 1;

        int size = Marshal.SizeOf(typeof(double)) * y0.Length;
        IntPtr y0Ptr = Marshal.AllocHGlobal(size);
        Marshal.Copy(y0, 0, y0Ptr, y0.Length);

        task.y0 = y0Ptr;


        double[] resultY = new double[task.n];

        // Создаём массив IntPtr для хранения указателей на делегаты
        IntPtr[] functionPointers = new IntPtr[2];

        // Фиксируем делегаты в памяти, чтобы GC не перемещал их
        GCHandle[] handles = new GCHandle[2];

        F1System[] system = new F1System[2];

        system[0] = SystemFunction1;
        system[1] = SystemFunction2;

        for (int i = 0; i < system.Length; i++)
        {
            handles[i] = GCHandle.Alloc(system[i]);
            functionPointers[i] = Marshal.GetFunctionPointerForDelegate(system[i]);
        }


        GCHandle functionsArrayHandle = GCHandle.Alloc(functionPointers, GCHandleType.Pinned);

        task.f = functionsArrayHandle.AddrOfPinnedObject();


        DiffUrSolver.DiffUrSolver.SolveDiffUrSystem(task, resultY);

        for (int i = 0; i < task.n; i++)
        {
            Console.WriteLine(resultY[i]);
        }

        Console.WriteLine(RealSystemFunction1(task.t1));
        Console.WriteLine(RealSystemFunction2(task.t1));

    }
    static void experiment4(DiffUrSolver.Methods method)
    {
        DiffUrSolver.SystemTask task = new DiffUrSolver.SystemTask();

        task.t0 = 0;
        task.t1 = 4.0;
        task.h = 0.001;
        task.method = method;
        task.n = 1;

        double[] y0 = new double[1];
        y0[0] = 0;
        //y0[1] = 1;

        task.y0 = DiffUrSolver.DiffUrSolver.DoubleArrayToIntPtr(y0);

        double[] resultY = new double[task.n];

        F1System[] system = new F1System[task.n];

        system[0] = SystemFunctionSingle;

        task.f = DiffUrSolver.DiffUrSolver.F1SystemToIntPtr(system);

        DiffUrSolver.DiffUrSolver.SolveDiffUrSystem(task, resultY);

        for (int i = 0; i < task.n; i++)
        {
            Console.WriteLine(resultY[i]);
        }

        Console.WriteLine(realFunc1(task.t1));

    }

    static void experiment5(DiffUrSolver.Methods method)
    {
        double c = 1.0;
        double m = 1.0;
        double k = 1.0;

        DiffUrSolver.SecondTask task = new DiffUrSolver.SecondTask();
        task.t0 = 0;
        task.t1 = 10.0;
        task.y_0 = 2;
        task.y1_0 = 0;
        task.h = 0.1;
        task.A = 1;
        task.B = c/m;
        task.C = k/m;
        task.D = 0;
        task.method = method;
        double result = 0;
        DiffUrSolver.DiffUrSolver.SolveSecondDiffUr(task, ref result);
        Console.WriteLine(result.ToString());
    }
    static void experiment5Arr(DiffUrSolver.Methods method)
    {
        double c = 4.47;
        double m = 10.0;
        double k = 50.0;

        DiffUrSolver.SecondTask task = new DiffUrSolver.SecondTask();
        task.t0 = 0;
        task.t1 = 10;
        task.y_0 = 0.01;
        task.y1_0 = 0;
        task.h = 0.1;
        task.A = 1;
        task.B = c / m;
        task.C = k / m;
        task.D = 0;
        task.method = method;
        double[] result = new double[10000];
        DiffUrSolver.DiffUrSolver.SolveSecondDiffUrArr(task, ref result);
        

        double start = task.t0;
        int j = 0;

        while (start < task.t1)
        {
            Console.WriteLine(result[j]);
            j++;
            start = task.t0 + j * task.h;
        }
    }
    static void Main(string[] args)
    {
        //experiment1(Methods.Euler);
        //experiment1(Methods.BackEuler);

        experiment5Arr(Methods.BackEuler);

        //experiment2(DiffUrSolver.Methods.RK4);

        //FullTestMethodsEps(func1,realFunc1);

        //FullTestMethodsEps(func1,realFunc1);

        //experiment4(Methods.Euler);

        
    }
}
